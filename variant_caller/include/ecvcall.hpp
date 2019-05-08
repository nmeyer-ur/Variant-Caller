/*
    Allele-specific variant caller
    Copyright (C) 2019  Nils Meyer, University of Regensburg, Germany
    <nils.meyer@ur.de>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/  // end legal

#ifndef _ECVCALL_H_
#define _ECVCALL_H_

#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <memory>
#include <cctype>
#include <cstdio>
#include "ecdefs.hpp"


class Variant_caller {

public:
    std::string region;
    int allele_count;
    position_t start, end;
    length_t length;
    std::string master_str;

    Variant_caller(const Parameters& p, const base_vector_t& ref,
        const cluster_map_t& cmap, const std::string& _region,
        const position_t& _start, const position_t& _end):
        region(_region), allele_count(0), start(_start),
        end(_end), params(p), cluster_map(cmap), reference(ref) {

        length = _end - _start;
        allele_count = cluster_map.size();
        init_alleles();
    }

    void variant_call(void) {

        // gather pool of variant positions
        auto& vpool = variant_position_pool;
        vpool.reserve(allele_count * _init_variant_position_pool);
        for(auto& allele : allele_vector) {
            allele.gather_variant_positions(vpool);
        }

        // unify positions
        std::sort(vpool.begin(), vpool.end());
        vpool.erase(std::unique(vpool.begin(), vpool.end()), vpool.end());

        // use a vector of vectors to store consensuses
        // Variant 1: [ Consensus Allele 1, Consensus Allele 2, ... ]
        // ...
        // Variant N: [ Consensus Allele 1, Consensus Allele 2, ... ]
        const auto position_count = vpool.size();
        consensus_pool_vector.resize(position_count);

        // reserve memory
        for(auto& pool : consensus_pool_vector) {
            pool.reserve(position_count);
        }

        // determine consensuses
        for(auto& allele : allele_vector) {
            allele.consensus(variant_position_pool, consensus_pool_vector);
        }


        // insertions
        consensus_pool_vector_t consensus_insertions_pool_vector;
        for(auto& allele : allele_vector) {
//            allele.consensus_insertions(consensus_insertions_pool_vector);
            allele.consensus_insertions(consensus_pool_vector);
        }
        // sort the insertions pool wrt positions
        std::sort(consensus_pool_vector.begin(), \
            consensus_pool_vector.end(), \
            [](const std::vector<Consensus_allele>& lhs, \
                const std::vector<Consensus_allele>& rhs) \
                { return (lhs[0].position < rhs[0].position); });
/*
        std::sort(consensus_insertions_pool_vector.begin(), \
            consensus_insertions_pool_vector.end(), \
            [](const std::vector<Consensus_allele>& lhs, \
                const std::vector<Consensus_allele>& rhs) \
                { return (lhs[0].position < rhs[0].position); });
*/
        _allele_coverage();
    }

    void tmpvcf(void) {
        _build_tmpvcf();
    }

    void display_tmpvcf(void) {
        _display_tmpvcf();
    }

    void gather_consensus_stats(void) {

        for(const auto& consensus_pool : consensus_pool_vector) {
            statistics.total_consensuses += (counter_t)consensus_pool.size();

            for(const auto& cs : consensus_pool) {
                if (cs.is_conclusive == true)
                    statistics.conclusive_consensuses += 1;
                if (cs.is_ambiguous == true)
                    statistics.ambiguous_consensuses += 1;

                if (cs.is_insertion == true)
                    statistics.total_consensuses_insertions += 1;
                if ((cs.is_insertion == true) && (cs.is_conclusive == true))
                    statistics.conclusive_consensuses_insertions += 1;
                if ((cs.is_insertion == true) && (cs.is_ambiguous == true))
                    statistics.ambiguous_consensuses_insertions += 1;

                if (cs.base == _char_del_base)
                    statistics.total_consensuses_deletions += 1;
                if ((cs.base == _char_del_base) && (cs.is_conclusive == true))
                    statistics.conclusive_consensuses_deletions += 1;
                if ((cs.base == _char_del_base) && (cs.is_conclusive == false))
                    statistics.inconclusive_consensuses_deletions += 1;
            }
        }
    }

    void update_statistics(Statistics& ext_statistics) {

        ext_statistics.positions_tested += statistics.positions_tested;
        ext_statistics.positions_with_variants += variant_position_pool.size();

        ext_statistics.total_consensuses += statistics.total_consensuses;
        ext_statistics.conclusive_consensuses += statistics.conclusive_consensuses;
        ext_statistics.ambiguous_consensuses += statistics.ambiguous_consensuses;

        ext_statistics.total_consensuses_insertions += statistics.total_consensuses_insertions;
        ext_statistics.conclusive_consensuses_insertions += statistics.conclusive_consensuses_insertions;
        ext_statistics.ambiguous_consensuses_insertions += statistics.ambiguous_consensuses_insertions;

        ext_statistics.total_consensuses_deletions += statistics.total_consensuses_deletions;
        ext_statistics.conclusive_consensuses_deletions += statistics.conclusive_consensuses_deletions;
        ext_statistics.inconclusive_consensuses_deletions += statistics.inconclusive_consensuses_deletions;

        const auto& ctr1 = statistics.allele_variant_counter_vector;
        for(size_t i=0; i != ctr1.size(); ++i) {
            ext_statistics.allele_variant_counter_vector[i] += ctr1[i];
        }

        const auto& ctr2 = statistics.allele_position_counter_vector;
        for(size_t i=0; i != ctr2.size(); ++i) {
            ext_statistics.allele_position_counter_vector[i] += ctr2[i];
        }

        ext_statistics.valid_block_length_after_clustering_vector.push_back((counter_t)statistics.positions_tested);
        ext_statistics.temp_variable += statistics.temp_variable;
        ext_statistics.block_region_map[region] += 1;
    }

private:
    Statistics statistics;
    const Parameters& params;
    const cluster_map_t& cluster_map;
    const base_vector_t& reference;
    allele_vector_t allele_vector;
    variant_position_vector_t variant_position_pool;
    consensus_pool_vector_t consensus_pool_vector;
    coverage_vector_t allele_coverage;

    void init_alleles(void) {

        allele_vector.reserve(allele_count);
        auto allele_id = 1;
        for(auto& mapping : cluster_map) {
            auto& read_vector = mapping.second;
            auto& tag = read_vector[0]->tag_encoded;
            allele_vector.emplace_back( \
                Allele(params, read_vector, tag, allele_id, reference, start));

            //allele.display_allele();
            ++allele_id;
        }
    }

    void _allele_coverage(void) {

        allele_coverage.resize(length);
        const auto min_reads = params.min_reads;
        for(const auto& allele : allele_vector) {
            const auto allele_offset = allele.start - start;
            const auto allele_length = allele.length;

            for(size_t i=0; i != allele_length; ++i) {
                if (allele.coverage[i] >= min_reads)
                    allele_coverage[allele_offset+i] += 1;
            }
        }

        for(size_t i=0; i != length; ++i) {
            const auto cov = allele_coverage[i];
            if (cov > 0) {
                statistics.positions_tested += 1;
                if (cov > _max_allele_position_counter)
                    statistics.allele_position_counter_vector[0] += 1;
                else
                    statistics.allele_position_counter_vector[cov] += 1;
            }
        }
    }

    void _display_tmpvcf(void) {

        fwrite(master_str.data(), 1, master_str.size(), stdout);
    }

    void _build_tmpvcf(void) {

        // snprintf limits the performance
        // if performance is critical, snprintf should be
        // coded by hand
        constexpr int str_reserve = 1000;

        master_str.reserve(str_reserve * str_reserve);

        std::string ref_str; ref_str.reserve(str_reserve);
        std::string alt_str; alt_str.reserve(str_reserve);
        std::string bases_str; bases_str.reserve(str_reserve * str_reserve);
        std::string phred_str; phred_str.reserve(str_reserve * str_reserve);
        std::string mapq_str; mapq_str.reserve(str_reserve);

        std::string bases; bases.reserve(str_reserve);
        std::string phred; phred.reserve(str_reserve);
        std::string mapq; mapq.reserve(str_reserve);

        constexpr int dummy_size = 100;
        char dummy[dummy_size+1];

        for(const auto& consensus_pool : consensus_pool_vector) {

            if (consensus_pool.size() == 0)
                continue;

            auto pos = consensus_pool[0].position;
            std::string position = std::to_string(pos);

            ref_str.clear();
            alt_str.clear();
            bases_str.clear();
            phred_str.clear();
            mapq_str.clear();

            bases.clear();
            phred.clear();
            mapq.clear();

            auto allele_ctr = 0;
            const bool is_insertion = consensus_pool[0].is_insertion;

            for(const auto& cs : consensus_pool) {

                ++allele_ctr;
                const int allele_id = (int)cs.allele_id;

#ifndef FAST_INT_TO_STRING
                auto n = snprintf(dummy, dummy_size, "(%d,%d,%d,%d,%d,%d;%d),",
                    (int)cs.base_count[_bin_A_base], \
                    (int)cs.base_count[_bin_T_base], \
                    (int)cs.base_count[_bin_C_base], \
                    (int)cs.base_count[_bin_G_base], \
                    (int)cs.base_count[_bin_N_base], \
                    (int)cs.base_count[_bin_del_base], \
                    allele_id);

                if (n < 0)
                    assert(!"VCF encoding error. Exiting.");
#else
                tuple_7(dummy,
                    (int)cs.base_count[_bin_A_base], \
                    (int)cs.base_count[_bin_T_base], \
                    (int)cs.base_count[_bin_C_base], \
                    (int)cs.base_count[_bin_G_base], \
                    (int)cs.base_count[_bin_N_base], \
                    (int)cs.base_count[_bin_del_base], \
                    allele_id);
#endif
                bases_str += dummy;

#ifndef FAST_INT_TO_STRING
                n = snprintf(dummy, dummy_size, "(%d,%d,%d,%d,%d,%d;%d),",
                    (int)cs.avg_phred[_bin_A_base], \
                    (int)cs.avg_phred[_bin_T_base], \
                    (int)cs.avg_phred[_bin_C_base], \
                    (int)cs.avg_phred[_bin_G_base], \
                    (int)cs.avg_phred[_bin_N_base], \
                    (int)cs.avg_phred[_bin_del_base], \
                    allele_id);

                if (n < 0)
                    assert(!"VCF encoding error. Exiting.");
#else
                tuple_7(dummy,
                    (int)cs.avg_phred[_bin_A_base], \
                    (int)cs.avg_phred[_bin_T_base], \
                    (int)cs.avg_phred[_bin_C_base], \
                    (int)cs.avg_phred[_bin_G_base], \
                    (int)cs.avg_phred[_bin_N_base], \
                    (int)cs.avg_phred[_bin_del_base], \
                    allele_id);
#endif
                phred_str += dummy;

#ifndef FAST_INT_TO_STRING
                n = snprintf(dummy, dummy_size, "(%d;%d),", \
                    (int)cs.avg_mapq, \
                    allele_id);

                if (n < 0)
                    assert(!"VCF encoding error. Exiting.");
#else
                tuple_2(dummy,
                    (int)cs.avg_mapq, \
                    allele_id);
#endif
                mapq_str += dummy;

                // consensus base
                char b = (char)cs.base;
                if (!cs.is_conclusive)
                    b = (char)tolower(b);

                if (cs.base == cs.reference) {
                    ref_str += b;
                    ref_str += ',';
                } else {
                    alt_str += b;
                    alt_str += ',';
                }
            }

            // yet empty strings
            if (ref_str.empty()) ref_str = ".,";
            if (alt_str.empty()) alt_str = ".,";

            // remove trailing ','
            ref_str.pop_back();
            alt_str.pop_back();
            bases_str.pop_back();
            phred_str.pop_back();
            mapq_str.pop_back();

            master_str += "TMPVARCALL ";
            master_str += region;
            master_str += '\t';
            master_str += position;
            master_str += '\t';
            master_str += (char)(consensus_pool[0].reference);
            master_str += '\t';
            master_str += ref_str;
            master_str += '\t';
            master_str += alt_str;
            master_str += '\t';
            master_str += std::to_string(allele_ctr);
            master_str += '\t';
            master_str += bases_str;
            master_str += '\t';
            master_str += phred_str;
            master_str += '\t';
            master_str += mapq_str;
            master_str += '\n';

            //cout << master_str;
            //master_str.clear();

            // update allele counter
            if (is_insertion == false) {
                auto& ctr = statistics.allele_variant_counter_vector;
                if (allele_ctr > _max_allele_variant_counter)
                    ctr[0] += 1;
                else
                    ctr[allele_ctr] += 1;
            }
        }
        //cout << master_str;
        statistics.temp_variable += 1;
        //fputs(master_str.c_str(), stdout);
        //fwrite(master_str.data(), 1, master_str.size(), stdout);
    }



    // Ultra fast integer to string conversion
    /*
        https://stackoverflow.com/questions/3982320/convert-integer-to-string-without-access-to-libraries

        buffer   pointer to char array
        index    array index
        n        integer number to be converted
        returns buffer index increment
    */
    int int_to_string(char* buffer, int index, int n) {

        int i = 0;

        // zero
        if (n == 0) {
            buffer[index + i++] = '0';
            return i;
        }

        // pile up, order is reversed
        while(n != 0) {
            buffer[index + i++] = (n%10) + '0';
            n = n / 10;
        }

        // reverse without additional buffer space
        for(int t = 0; t < i / 2; t++) {
            const int index1 = index + t;
            const int index2 = index + i - t - 1;
            std::swap(buffer[index1], buffer[index2]);
        }

        return i;
    }

    int tuple_2(char* buffer, int i1, int i2) {

        int index = 0;
        buffer[index++] = '(';
        index += int_to_string(buffer, index, i1);
        buffer[index++] = ';';
        index += int_to_string(buffer, index, i2);
        buffer[index++] = ')';
        buffer[index++] = ',';
        buffer[index++] = '\0';

        return index;
    }

    int tuple_7(char* buffer, int i1, int i2, int i3, int i4, int i5, int i6, int i7) {

        int index = 0;
        buffer[index++] = '(';
        index += int_to_string(buffer, index, i1);
        buffer[index++] = ',';
        index += int_to_string(buffer, index, i2);
        buffer[index++] = ',';
        index += int_to_string(buffer, index, i3);
        buffer[index++] = ',';
        index += int_to_string(buffer, index, i4);
        buffer[index++] = ',';
        index += int_to_string(buffer, index, i5);
        buffer[index++] = ',';
        index += int_to_string(buffer, index, i6);
        buffer[index++] = ';';
        index += int_to_string(buffer, index, i7);
        buffer[index++] = ')';
        buffer[index++] = ',';
        buffer[index++] = '\0';

        return index;
    }
};


#endif // _ECVCALL_H_
