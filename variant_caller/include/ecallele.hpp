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

#ifndef _ECALLELE_H_
#define _ECALLELE_H_

#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <array>
#include "ecdefs.hpp"


class Consensus_allele {

public:
    const id_t allele_id;
    const position_t position;
    const int insertion_index;
    const base_t reference;
    int avg_mapq;
    coverage_t coverage;
    base_t base;
    std::array<int, _num_bases> base_count;
    std::array<int, _num_bases> avg_phred;
    bool is_conclusive;
    bool is_ambiguous;
    bool is_insertion;

    Consensus_allele(const Parameters& p, const id_t _allele_id, const position_t _position,
        coverage_t _coverage, const base_t _reference, const bool _is_insertion,
        const int _insertion_index):
        allele_id(_allele_id), position(_position), insertion_index(_insertion_index),
        reference(_reference), avg_mapq(0), coverage(_coverage), base(_char_del_base),
        is_insertion(_is_insertion), params(p) {

        for(auto i=0; i<_num_bases; ++i) {
            base_count[i] = 0;
            avg_phred[i] = 0;
        }
    }

    void update(const int base_index, const phred_t phred_score, const mapq_t mapq) {
        base_count[base_index] += 1;
        avg_phred[base_index] += (int)phred_score;
        avg_mapq += (int)mapq;
    }

    void finalize(void) {
        const float cov = (float)coverage;
        avg_mapq = (int)((float)avg_mapq / cov);

        // determine average phred score, consensus base and ambiguity
        bool conclusive = false;
        bool ambiguous = false;
        int cs_base_count = 0;
        for(auto i=0; i<_num_bases; ++i) {
            if (base_count[i] > 0)
                avg_phred[i] = \
                    (phred_t)((float)avg_phred[i] / (float)base_count[i]);

            if (cs_base_count < base_count[i]) {
                cs_base_count = base_count[i];
                base = bin_base_to_char_lookup(i);
                ambiguous = false;
            } else if (cs_base_count == base_count[i]) {
                base = bin_base_to_char_lookup(_bin_N_base);
                ambiguous = true;
            }
        }

        // conclusiveness of the consensus base
        if (ambiguous == false) {
            const float thres = params.base_acceptance_threshold;
            if ((float)cs_base_count >= thres * cov)
                conclusive = true;
        }

        is_conclusive = conclusive;
        is_ambiguous = ambiguous;
    }


private:
    const Parameters& params;

};

class Allele {

public:
    id_t id;
    position_t start, end;
    length_t length;
    size_t read_count;
    tag_encoded_t tag;
    coverage_vector_t coverage;
    int max_depth;
    variant_position_vector_t variant_positions;

    Allele(const Parameters& p, const read_vector_t& rv,
        const tag_encoded_t& cluster_tag, const int& allele_id,
        const base_vector_t& _reference, const int& ref_start):
        id(allele_id), start(0), end(0), length(0), tag(cluster_tag),
        max_depth(0), params(p), read_vector(rv), reference(_reference) {

        read_count = read_vector.size();
        init_allele();

        ref_offset = start - ref_start;

        if (params.verbosity >= 4)
            display_allele();
    }

    void gather_variant_positions(variant_position_vector_t& variant_positions) {
        _gather_variant_positions(variant_positions);
    }

    void consensus(const variant_position_vector_t& variant_positions,
        consensus_pool_vector_t& consensus_pool_vector) {
        _consensus(variant_positions, consensus_pool_vector);
    }

    void consensus_insertions(consensus_pool_vector_t& consensus_insertions_pool_vector) {
        _consensus_insertions(consensus_insertions_pool_vector);
    }

    void display_allele(void) {

        std::cout << "Allele " << id << " # reads "<< read_vector.size() << '\n';
        std::cout << tag_encoded_to_string(tag) << "  -> start " << start << "   end " << end << '\n';

        std::string ref = base_vector_to_string(reference);
        std::cout << start - ref_offset << " REF -- - " << tag_encoded_to_string(0) << ' ' << ref << '\n';

        // reorder reads for improved visualization
        // we do it indirectly
        std::vector<std::pair<int, position_t>> reordering(read_vector.size());
        for(size_t i=0; i != reordering.size(); ++i) {
            reordering[i].first = i;
            reordering[i].second = read_vector[i]->start;
        }
        std::sort(reordering.begin(), reordering.end(), \
            [](std::pair<int, position_t>& a, std::pair<int, position_t>& b)
                { return a.second < b.second; });

        for(size_t i=0; i != reordering.size(); ++i) {
            const auto index = reordering[i].first;
            const auto delta = read_vector[index]->start - start + ref_offset;
            const auto tag   = read_vector[index]->tag_encoded;
			const char direction = read_vector[index]->is_forward ? 'F' : 'R';
			std::string orientation;
			if (read_vector[index]->orientation.test(_orientation_index_FR))
				orientation = "FR";
			else if (read_vector[index]->orientation.test(_orientation_index_RF))
				orientation = "RF";
			else if (read_vector[index]->orientation.test(_orientation_index_FF))
				orientation = "FF";
			else if (read_vector[index]->orientation.test(_orientation_index_RR))
				orientation = "RR";
			else
				orientation = "--";
            std::string aln = base_vector_to_color_string(read_vector[index]);
            aln.insert(0, delta, ' ');
            std::cout << reordering[i].second << " ALN " << orientation << ' ' << direction << ' ' \
				<< tag_encoded_to_string(tag) << ' ' << aln << '\n';
        }

        /*
        for(const auto& read : read_vector) {
            auto delta = read->start - start + ref_offset;
            //string aln = base_vector_to_string(read->base);
            string aln = base_vector_to_color_string(read);
            // TODO there might be a better way to insert spaces
            aln.insert(0, delta, ' ');
            //cout << read->start << ' ' << read->id << ' ' << aln << '\n';
            cout << read->start << " ALN " << aln << '\n';
        }
        */

        std::cout << "Max depth of coverage " << max_depth << '\n';
    }

    std::string base_vector_to_color_string(const read_ptr_t read) {

        const auto read_offset = read->start - start;
        const auto& read_base  = read->base;

        std::string s;
        for(size_t i=0; i != read_base.size(); ++i) {
            const auto offset = ref_offset+read_offset+i;

            // FIXME
            if (coverage[read_offset+i] < params.min_reads) {
                s += "\x1B[1;31m-\x1B[0m";
                continue;
            }

            const char c = read_base[i];
            if (c == reference[offset]) {
                s += '-';
            } else {
                s += "\x1B[1;31m";
                s += c;
                s += "\x1B[0m";
            }
        }

        return s;
    }

    void display_variants(void) {

        const auto num_variants = variant_positions.size();
        std::cout << "Allele " << id << " variants " << num_variants << '\n';
        if (num_variants > 0) {
            for (const auto& variant : variant_positions)
                std::cout << variant << ' ';
            std::cout << '\n';
        }
    }

private:
    const Parameters& params;
    const read_vector_t& read_vector;
    const base_vector_t& reference;
    position_t ref_offset;
    std::vector<int> deviations;

    void init_allele(void) {

        // evaluate allele boundaries
        start = read_vector[0]->start;
        for(const auto& read : read_vector) {
            auto s = read->start;
            auto e = read->end;
            start = std::min(s, start);
            end   = std::max(e, end);
        }

        length = end - start;

        // evaluate the depth of coverage
        // first create step counter
        coverage_vector_t delta(length+1, 0);
        for(const auto& read : read_vector) {
            auto s = read->start - start;
            auto e = read->end - start;
            delta[s] += 1;
            delta[e] -= 1;
        }

        // next evaluate coverage
        coverage.reserve(length);
        auto cum_sum = 0;
        for(auto it=delta.begin(); it != --delta.end(); ++it) {
            cum_sum += *it;
            max_depth = std::max(max_depth, cum_sum);
            coverage.push_back(cum_sum);
        }

        // this should never happen
        if (max_depth < params.min_reads) {
            std::cout << read_count << '\n';
            std::cout << max_depth << '\n';
            assert(!"Coverage in allele below threshold. Exiting.");
        }
    }

    void _gather_variant_positions(variant_position_vector_t& variant_positions) {

        deviations.resize(length);

        // count deviations from reference
        for(const auto& read : read_vector) {

// late read data decode
#ifdef LATE_READ_DATA_DECODING
            read->process_read2();
#endif

            const auto read_offset = read->start - start;
            const auto read_length = read->length;
            const auto& read_base  = read->base;

            for(size_t i=0; i != read_length; ++i) {
                if (read_base[i] != reference[ref_offset+read_offset+i])
                    deviations[read_offset+i] += 1;
            }
        }

        // check deviations
        // store absolute position in variant position pool if
        // threshold exceeded
        const auto min_reads = params.min_reads;
        for(size_t i=0; i != length; ++i) {
            const auto dev = deviations[i];
            if (dev == 0)
                continue;

            const auto cov = coverage[i];
            if (cov < min_reads)
                continue;

            const float thres = params.base_acceptance_threshold;
            const float cov_float = (float)cov;
            const float n_ref = (float)(cov - dev);
            if (n_ref < (thres * cov_float)) {
                position_t position = start + i;
                variant_positions.push_back(position);
            }
        }
    }

    void _consensus(const variant_position_vector_t& variant_positions,
        consensus_pool_vector_t& consensus_pool_vector) {

        // generate position index pool
        std::vector<int> variant_positions_index_pool;
        variant_positions_index_pool.reserve(variant_positions.size());

        const auto min_reads = params.min_reads;
        for(size_t index=0; index != variant_positions.size(); ++index) {

            const auto position = variant_positions[index];

            // boundary check
            if (position < start)
                continue;
            if (position >= end)
                break;

            // coverage check
            const auto allele_index = position - start;
            const auto cov = coverage[allele_index];
            if (cov < min_reads)
                continue;

            variant_positions_index_pool.push_back(index);
            const auto ref = reference[ref_offset+allele_index];
            consensus_pool_vector[index].emplace_back(\
                Consensus_allele(params, id, position, cov, ref, false, -1));
        }

        // return if no variants in allele
        if (variant_positions_index_pool.size() == 0)
            return;

        // loop over reads
        for(const auto& read : read_vector) {

            const auto read_start = read->start;
            const auto read_end = read->end;
            const auto mapq = read->mapq;

            // loop over variant positions
            for(const auto index : variant_positions_index_pool) {

                const auto position = variant_positions[index];
                if (position < read_start)
                    continue;
                if (position >= read_end)
                    break;

                const auto read_index = position - read_start;
                const auto base_index = char_base_to_bin_lookup(read->base[read_index]);
                const auto phred = read->phred[read_index];

                auto& cs = consensus_pool_vector[index].back();
                cs.update(base_index, phred, mapq);
            }
        }

        // finalize consensuses
        for(const auto index : variant_positions_index_pool) {
            auto& cs = consensus_pool_vector[index].back();
            cs.finalize();
        }

    }

    void _consensus_insertions(consensus_pool_vector_t& consensus_insertions_pool_vector) {

        insertions_vector_map_t insertions_map;

        // create a look-up table of insertions if coverage is sufficient
        // Position 1: [ Insertions read 1, Insertions read 2, ... ]
        // ...
        // Position N: [ Insertions read 1, Insertions read 2, ... ]
        const auto min_reads = params.min_reads;
        for(const auto& read : read_vector) {
            const auto& ins = read->insertions;
            for(const auto& entry : ins) {
                const auto position = entry.first;
                const auto allele_index = position - start;
                const auto cov = coverage[allele_index];

                if (cov < min_reads)
                    continue;

                insertions_map[position].push_back(entry.second);
            }
        }

        if (insertions_map.size() == 0)
            return;

        // consensus
        for(const auto& entry : insertions_map) {
            const auto& ins_vector = entry.second;
            const auto position = entry.first;

            // threshold check, treat insertions as deviations
            const auto allele_index = position - start;
            const auto cov = coverage[allele_index];
            const auto dev = ins_vector.size();
            const float thres = params.base_acceptance_threshold;
            const float cov_float = (float)cov;
            const float n_ref = (float)(cov - dev);

            // remove noise
            if (n_ref >= (thres * cov_float))
                continue;

            //std::cout << position << ' ';

            int mapq = 0;
            int base_index = 0;
            std::vector<Consensus_allele> insertions_vector;
            while(true) {
                Consensus_allele cs(params, id, position, cov, _char_del_base, true, base_index);
                int ctr = 0;
                for(auto ins_index = 0; ins_index != (int)dev; ++ins_index) {
                    const auto& insertion = ins_vector[ins_index];
                    auto base = _bin_del_base;
                    auto phred = _phred_del_base;
                    if (base_index < (int)insertion.base.size()) {
                        base = char_base_to_bin_lookup(insertion.base[base_index]);
                        phred = insertion.phred[base_index];
                        ++ctr;
                    }
                    cs.update(base, phred, 0);
                }

                // break if insufficient insertion coverage
/*
                if (n_ref >= (thres * (float)ctr)) {
                    std::cout << " noise";
                    break;
                }
*/
                // complete missing deletions in consensus
                const int reads_without_insertions = cov - dev;
                cs.base_count[_bin_del_base] += reads_without_insertions;
                cs.avg_phred[_bin_del_base] += reads_without_insertions * _phred_del_base;

                // determine mapq score, needs to be done only once for each position
                if (base_index == 0) {
                    for(const auto& read : read_vector) {
                        if (position < read->start)
                            continue;
                        if (position >= read->end)
                            continue;
                        mapq += (int)read->mapq;
                    }
                }

                cs.avg_mapq = mapq;
                cs.finalize();

                // break if consensus base is a deletion
                if (cs.base == _char_del_base)
                    break;

                insertions_vector.push_back(cs);

                //std::cout << cs.base;
                ++base_index;
            }

            if (insertions_vector.size() > 0)
                consensus_insertions_pool_vector.push_back(insertions_vector);

            //std::cout << '\n';
        }
    }

};

#endif // _ECALLELE_H_
