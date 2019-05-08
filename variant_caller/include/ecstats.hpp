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

#ifndef _ECSTATS_H_
#define _ECSTATS_H_


#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <zlib.h> // for crc32

class S {

public:
    counter_t n, min, max, q1, q2, q3, sum;
    double mean, stddev, stderr, var;

    S(counter_vector_t& v):
        min(0), max(0), q1(0), q2(0), q3(0), sum(0),
        mean(0.), stddev(0.), stderr(0.), var(0.) {

        n = v.size();
        if (n == 0)
            return;

        sort(v.begin(), v.end());
        min = v[0];
        max = v[n-1];
        q1 = v[(int)(n / 4.)];
        q2 = v[(int)(2*n / 4.)];
        q3 = v[(int)(3*n / 4.)];

        sum = 0;
        for(const auto& el : v)
            sum += el;
        mean = (double)sum / n;

        var = 0.;
        for(const auto& el : v)
            var += ((double)el - mean) * ((double)el - mean);
        var /= n;

        stddev = sqrt(var);
        stderr = stddev / sqrt((double)n);
    }
};

class Statistics {

public:
    counter_t regions;
    counter_t total_blocks;
    counter_t valid_blocks;
    counter_t total_reads_in_alignment_file;
    counter_t valid_reads_in_alignment_file;
    counter_t valid_reads_in_valid_blocks;
    counter_t total_reads_in_blocks;
    counter_t tag_before_clustering;
    counter_t tag_after_clustering;
    counter_t reads_before_clustering;
    counter_t reads_after_clustering;
    counter_t reference_bases_loaded_from_file;
    counter_t positions_tested;
    counter_t positions_with_variants;
    counter_t leven_complete_evaluations;
    counter_t leven_total_evaluations;
    counter_t total_consensuses;
    counter_t conclusive_consensuses;
    counter_t ambiguous_consensuses;
    counter_t total_consensuses_insertions;
    counter_t conclusive_consensuses_insertions;
    counter_t ambiguous_consensuses_insertions;
    counter_t total_consensuses_deletions;
    counter_t conclusive_consensuses_deletions;
    counter_t inconclusive_consensuses_deletions;
    counter_t temp_variable;  // useful for debugging

    counter_vector_t valid_block_length_before_clustering_vector;
    counter_vector_t valid_block_length_after_clustering_vector;
    counter_vector_t allele_variant_counter_vector;
    counter_vector_t allele_position_counter_vector;
    counter_vector_t tag_before_clustering_vector;
    counter_vector_t tag_after_clustering_vector;
    tag_counter_map_t tag_counter_map;
    std::map<std::string, counter_t> block_region_map;


    Statistics():
        regions(0),
        total_blocks(0),
        valid_blocks(0),
        total_reads_in_alignment_file(0),
        valid_reads_in_alignment_file(0),
        valid_reads_in_valid_blocks(0),
        total_reads_in_blocks(0),
        tag_before_clustering(0),
        tag_after_clustering(0),
        reads_before_clustering(0),
        reads_after_clustering(0),
        reference_bases_loaded_from_file(0),
        positions_tested(0),
        positions_with_variants(0),
        leven_complete_evaluations(0),
        leven_total_evaluations(0),
        total_consensuses(0),
        conclusive_consensuses(0),
        ambiguous_consensuses(0),
        total_consensuses_insertions(0),
        conclusive_consensuses_insertions(0),
        ambiguous_consensuses_insertions(0),
        total_consensuses_deletions(0),
        conclusive_consensuses_deletions(0),
        inconclusive_consensuses_deletions(0),
        temp_variable(0)
    {
        allele_variant_counter_vector.resize(_max_allele_variant_counter+1);
        allele_position_counter_vector.resize(_max_allele_position_counter+1);
    }

    void display_stats(metadata_t metadata) {

        std::sort(metadata.begin(), metadata.end(), \
            [](Alignment_metadata a, Alignment_metadata b){ return a.index < b.index; });
        _display_stats(metadata);
    }

    void display_crc32() { _crc32(); }

    Statistics& operator += (const Statistics& rhs) {

        this->total_blocks += rhs.total_blocks;
        this->valid_blocks += rhs.valid_blocks;
        this->total_reads_in_alignment_file += rhs.total_reads_in_alignment_file;
        this->valid_reads_in_alignment_file += rhs.valid_reads_in_alignment_file;
        this->valid_reads_in_valid_blocks += rhs.valid_reads_in_valid_blocks;
        this->total_reads_in_blocks += rhs.total_reads_in_blocks;
        this->tag_before_clustering += rhs.tag_before_clustering;
        this->tag_after_clustering += rhs.tag_after_clustering;
        this->reads_before_clustering += rhs.reads_before_clustering;
        this->reads_after_clustering += rhs.reads_after_clustering;
        this->reference_bases_loaded_from_file += rhs.reference_bases_loaded_from_file;
        this->positions_tested += rhs.positions_tested;
        this->positions_with_variants += rhs.positions_with_variants;
        this->leven_complete_evaluations += rhs.leven_complete_evaluations;
        this->leven_total_evaluations += rhs.leven_total_evaluations;
        this->total_consensuses += rhs.total_consensuses;
        this->conclusive_consensuses += rhs.conclusive_consensuses;
        this->ambiguous_consensuses += rhs.ambiguous_consensuses;
        this->total_consensuses_insertions += rhs.total_consensuses_insertions;
        this->conclusive_consensuses_insertions += rhs.conclusive_consensuses_insertions;
        this->ambiguous_consensuses_insertions += rhs.ambiguous_consensuses_insertions;
        this->total_consensuses_deletions += rhs.total_consensuses_deletions;
        this->conclusive_consensuses_deletions += rhs.conclusive_consensuses_deletions;
        this->inconclusive_consensuses_deletions += rhs.inconclusive_consensuses_deletions;
        this->temp_variable += rhs.temp_variable;

        for(size_t i=0; i != this->allele_variant_counter_vector.size(); ++i) {
            this->allele_variant_counter_vector[i] += rhs.allele_variant_counter_vector[i];
        }

        for(size_t i=0; i != this->allele_position_counter_vector.size(); ++i) {
            this->allele_position_counter_vector[i] += rhs.allele_position_counter_vector[i];
        }

        for(const auto& mapping : rhs.block_region_map)
            this->block_region_map[mapping.first] += mapping.second;

        this->valid_block_length_before_clustering_vector.insert(\
            this->valid_block_length_before_clustering_vector.end(),
            rhs.valid_block_length_before_clustering_vector.begin(),
            rhs.valid_block_length_before_clustering_vector.end());

        this->valid_block_length_after_clustering_vector.insert(\
            this->valid_block_length_after_clustering_vector.end(),
            rhs.valid_block_length_after_clustering_vector.begin(),
            rhs.valid_block_length_after_clustering_vector.end());

        this->tag_before_clustering_vector.insert(\
            this->tag_before_clustering_vector.end(),
            rhs.tag_before_clustering_vector.begin(),
            rhs.tag_before_clustering_vector.end());

        this->tag_after_clustering_vector.insert(\
            this->tag_after_clustering_vector.end(),
            rhs.tag_after_clustering_vector.begin(),
            rhs.tag_after_clustering_vector.end());

        for(const auto& mapping : rhs.tag_counter_map) {
            this->tag_counter_map[mapping.first] += mapping.second;
        }

        return *this;
    }

private:
    void _crc32(void) {

        counter_vector_t scratch(1000);
        scratch.push_back((counter_t)regions);
        scratch.push_back((counter_t)total_blocks);
        scratch.push_back((counter_t)valid_blocks);
        scratch.push_back((counter_t)total_reads_in_alignment_file);
        scratch.push_back((counter_t)valid_reads_in_alignment_file);
        scratch.push_back((counter_t)valid_reads_in_valid_blocks);
        scratch.push_back((counter_t)total_reads_in_blocks);
        scratch.push_back((counter_t)tag_before_clustering);
        scratch.push_back((counter_t)tag_after_clustering);
        scratch.push_back((counter_t)reference_bases_loaded_from_file);
        scratch.push_back((counter_t)positions_tested);
        scratch.push_back((counter_t)positions_with_variants);
        scratch.push_back((counter_t)leven_complete_evaluations);
        scratch.push_back((counter_t)leven_total_evaluations);
        scratch.push_back((counter_t)total_consensuses);
        scratch.push_back((counter_t)conclusive_consensuses);
        scratch.push_back((counter_t)ambiguous_consensuses);
        scratch.push_back((counter_t)temp_variable);
        scratch.push_back((counter_t)total_consensuses_insertions);
        scratch.push_back((counter_t)conclusive_consensuses_insertions);
        scratch.push_back((counter_t)ambiguous_consensuses_insertions);
        scratch.push_back((counter_t)total_consensuses_deletions);
        scratch.push_back((counter_t)conclusive_consensuses_deletions);
        scratch.push_back((counter_t)inconclusive_consensuses_deletions);

        unsigned int result = crc32(0x12345678, (const unsigned char*)scratch.data(), scratch.size() * sizeof(counter_t));
        std::cout << "CRC32  0x" << std::setfill('0') << std::setw(8) << std::hex << result << '\n';
    }

    void _display_stats(const metadata_t& metadata) {

        constexpr int w = 12;

        std::cout << '\n' << header("Summary") << "\n\n";

        std::cout << "Selected regions                            " \
            << std::setw(w) << block_region_map.size() << '\n';

        std::cout << "Total reads in regions                      " \
            << std::setw(w) << total_reads_in_alignment_file << '\n';

        std::cout << "Reads after filtering                       " \
            << std::setw(w) << valid_reads_in_valid_blocks \
            << frac(valid_reads_in_valid_blocks, total_reads_in_alignment_file, "of total reads") << '\n';

        std::cout << "Block count before clustering               " \
            << std::setw(w) << valid_block_length_before_clustering_vector.size() << '\n';

        std::cout << "Block count after clustering                " \
            << std::setw(w) << valid_block_length_after_clustering_vector.size() << '\n';

        std::cout << "Length    | blocks\tmin\tq1\tmedian\tq3\tmax\tsum\tmean\tstddev\tstderr\tvar\n";

        auto b1 = S(valid_block_length_before_clustering_vector);
        std::cout << "  before  | " << b1.n << '\t' << b1.min << '\t' << b1.q1 << '\t' << b1.q2 << '\t' << b1.q3 << '\t' << b1.max \
            << '\t' << b1.sum << '\t' << std::setprecision(4) << b1.mean << '\t' << b1.stddev << '\t' << b1.stderr<< '\t' << b1.var << '\n';

        b1 = S(valid_block_length_after_clustering_vector);
        std::cout << "  after   | " << b1.n << '\t' << b1.min << '\t' << b1.q1 << '\t' << b1.q2 << '\t' << b1.q3 << '\t' << b1.max \
            << '\t' << b1.sum << '\t' << std::setprecision(4) << b1.mean << '\t' << b1.stddev << '\t' << b1.stderr<< '\t' << b1.var << '\n';

        std::cout << '\n';
        std::cout << "  Region          block count after clustering\n";

        for(size_t i = 0; i != metadata.size(); ++i) {
            const std::string& region = metadata[i].name;
            const counter_t count = block_region_map[region];
            std::cout << std::setw(w) << region << "    " << std::setw(w) << count \
                << frac(count, valid_block_length_after_clustering_vector.size(), "of total blocks") << '\n';
        }
        std::cout << '\n';


        std::cout << header("Clustering statistics") << "\n\n";

        const counter_t leven_et = leven_total_evaluations - leven_complete_evaluations;

        std::cout << "Total Levenshtein distances evaluated       " \
                << std::setw(w) << leven_total_evaluations << '\n';

        std::cout << "Levenshtein distances evaluated (NW)        " \
                << std::setw(w) << leven_complete_evaluations \
                << frac(leven_complete_evaluations, leven_total_evaluations, "of total distances") << '\n';

        std::cout << "Levenshtein distances exceeding lower bound " \
                << std::setw(w) << leven_et \
                << frac(leven_et, leven_total_evaluations, "of total distances") << '\n';

        std::cout << "Alleles before TAG clustering               " \
            << std::setw(w) << tag_before_clustering << '\n';

        const counter_t ared = tag_before_clustering - tag_after_clustering;

        std::cout << "Alleles after TAG clustering and filtering  " \
            << std::setw(w) << tag_after_clustering \
            << frac(ared, tag_before_clustering, "reduction") << '\n';

        std::cout << "Alleles   | blocks\tmin\tq1\tmedian\tq3\tmax\tsum\tmean\tstddev\tstderr\tvar\n";

        auto s1 = S(tag_before_clustering_vector);
        std::cout << "  before  | " << s1.n << '\t' << s1.min << '\t' << s1.q1 << '\t' << s1.q2 << '\t' << s1.q3 << '\t' << s1.max \
            << '\t' << s1.sum << '\t' << std::setprecision(4) << s1.mean << '\t' << s1.stddev << '\t' << s1.stderr<< '\t' << s1.var << '\n';
        s1 = S(tag_after_clustering_vector);
        std::cout << "  after   | " << s1.n << '\t' << s1.min << '\t' << s1.q1 << '\t' << s1.q2 << '\t' << s1.q3 << '\t' << s1.max \
            << '\t' << s1.sum << '\t' << std::setprecision(4) << s1.mean << '\t' << s1.stddev << '\t' << s1.stderr<< '\t' << s1.var << '\n';

        std::cout << "Reads before TAG clustering                 " \
            << std::setw(w) << reads_before_clustering << '\n';

        const counter_t rred = reads_before_clustering - reads_after_clustering;

        std::cout << "Reads after TAG clustering and filtering    " \
            << std::setw(w) << reads_after_clustering \
            << frac(rred, reads_before_clustering, "reduction") << '\n';

        // TAG statistics
        std::vector<std::pair<tag_encoded_t, counter_t> > tag_counter_vector(tag_counter_map.begin(), tag_counter_map.end());

        sort(tag_counter_vector.begin(), tag_counter_vector.end(), \
            [](const std::pair<tag_encoded_t, counter_t>& a, const std::pair<tag_encoded_t, counter_t>& b) \
            { return a.second > b.second; } );

        std::cout << "TAG sequences after clustering and filtering" \
            << std::setw(w) << tag_counter_vector.size() << '\n';

        std::cout << '\n';
        std::cout << "  TAG sequence      block count             " << '\n';

        auto c = 0;
        for(const auto& p : tag_counter_vector) {
            std::cout << "  " << tag_encoded_to_string(p.first) << "  " << std::setw(w) << p.second << '\n';
            c++;
            if (c > 10) break;
        }

        std::cout << '\n' << header("Reference statistics") << "\n\n";

        std::cout << "Reference bases loaded from file            " \
            << std::setw(w) << reference_bases_loaded_from_file << '\n';

        std::cout << "Reference bases tested                      " \
            << std::setw(w) << positions_tested \
            << frac(positions_tested, reference_bases_loaded_from_file, "of loaded reference bases") << '\n';


        std::cout << '\n' << header("Global allele statistics (w/o insertions)") << "\n\n";

        counter_t one_allele  = 0;
        counter_t two_alleles = 0;
        counter_t at_least_three_alleles = 0;

        if (allele_position_counter_vector.size() > 0) {
            one_allele = allele_position_counter_vector[1];
        }

        if (allele_position_counter_vector.size() > 1) {
            two_alleles = allele_position_counter_vector[2];
        }

        if (allele_position_counter_vector.size() > 2) {
            const auto& ctr = allele_position_counter_vector;
            at_least_three_alleles = ctr[0];
            for(size_t i=3; i != ctr.size(); ++i)
                at_least_three_alleles += ctr[i];
        }

        const counter_t tot_alleles = one_allele + two_alleles + at_least_three_alleles;

        std::cout << "Positions comprising one allele             " \
            << std::setw(w) << one_allele \
            << frac(one_allele, tot_alleles, "of tested positions") << '\n';

        std::cout << "Positions comprising two alleles            " \
            << std::setw(w) << two_alleles \
            << frac(two_alleles, tot_alleles, "of tested positions") << '\n';

        std::cout << "Positions comprising at least three alleles " \
            << std::setw(w) << at_least_three_alleles \
            << frac(at_least_three_alleles, tot_alleles, "of tested positions") << '\n';

        std::cout << "Positions comprising allelic dropout        " \
            << std::setw(w) << one_allele \
            << frac(one_allele, tot_alleles, "of tested positions") << '\n';


        if (allele_position_counter_vector.size() > 0) {
            const auto digits = (_max_allele_position_counter > 0) ? \
                (int)log10((double)_max_allele_position_counter) + 1 : 1;
            const auto& ctr = allele_position_counter_vector;
            counter_t total = 0;
            for(size_t i=0; i != ctr.size(); ++i)
                total += ctr[i];

            std::cout << '\n';
            std::cout << "  Allele count  Positions\n";

            for(size_t i=1; i != ctr.size(); ++i) {
                const auto entry = ctr[i];

                std::cout << "    " << std::setw(digits) << i << "    " << std::setw(w) << entry \
                    << frac(entry, total, "of tested positions") << '\n';
            }
            const auto entry = ctr[0];
            std::cout << " >= " << std::setw(digits) << _max_allele_position_counter+1 \
                << "    " << std::setw(w) << entry \
                << frac(entry, total, "of tested positions") << '\n';
        }


        std::cout << '\n' << header("Variants statistics") << "\n\n";

        std::cout << "Positions comprising variants               " \
            << std::setw(w) << positions_with_variants \
            << frac(positions_with_variants, positions_tested, "of tested positions") << '\n';


        std::cout << "Total consensuses count                     " \
            << std::setw(w) << total_consensuses << '\n';

        std::cout << "  Insertions count                          " \
            << std::setw(w) << total_consensuses_insertions \
            << frac(total_consensuses_insertions, total_consensuses, "of total consensuses") << '\n';

        std::cout << "  Deletions count                           " \
            << std::setw(w) << total_consensuses_deletions \
            << frac(total_consensuses_deletions, total_consensuses, "of total consensuses") << '\n';

        const auto mutations = total_consensuses - total_consensuses_insertions - total_consensuses_deletions;
        std::cout << "  SNPs count                                " \
            << std::setw(w) << mutations \
            << frac(mutations, total_consensuses, "of total consensuses") << '\n';


        std::cout << "Total conclusive consensuses count          " \
            << std::setw(w) << conclusive_consensuses \
            << frac(conclusive_consensuses, total_consensuses, "of total consensuses") << '\n';

        std::cout << "  Insertions count                          " \
            << std::setw(w) << conclusive_consensuses_insertions \
            << frac(conclusive_consensuses_insertions, total_consensuses_insertions, "of total insertions consensuses") << '\n';

        std::cout << "  Deletions count                           " \
            << std::setw(w) << conclusive_consensuses_deletions \
            << frac(conclusive_consensuses_deletions, total_consensuses_deletions, "of total deletions consensuses") << '\n';

        const auto conclusive_mutations = \
            conclusive_consensuses - conclusive_consensuses_insertions - conclusive_consensuses_deletions;
        std::cout << "  SNPs count                                " \
            << std::setw(w) << conclusive_mutations \
            << frac(conclusive_mutations, mutations, "of total SNPs consensuses") << '\n';


        const auto inconclusive_consensuses = \
            total_consensuses - conclusive_consensuses - ambiguous_consensuses;
        std::cout << "Total inconclusive consensuses count        " \
            << std::setw(w) << inconclusive_consensuses \
            << frac(inconclusive_consensuses, total_consensuses, "of total consensuses") << '\n';

        const auto inconclusive_consensuses_insertions = \
            total_consensuses_insertions - conclusive_consensuses_insertions - ambiguous_consensuses_insertions;
        std::cout << "  Insertions count                          " \
            << std::setw(w) << inconclusive_consensuses_insertions \
            << frac(inconclusive_consensuses_insertions, total_consensuses_insertions, "of total insertions consensuses") << '\n';

        const auto inconclusive_consensuses_deletions = \
            total_consensuses_deletions - conclusive_consensuses_deletions;
        std::cout << "  Deletions count                           " \
            << std::setw(w) << inconclusive_consensuses_deletions \
            << frac(inconclusive_consensuses_deletions, total_consensuses_deletions, "of total deletions consensuses") << '\n';

        const auto inconclusive_mutations = \
            inconclusive_consensuses - inconclusive_consensuses_insertions - inconclusive_consensuses_deletions;
        std::cout << "  SNPs count                                " \
            << std::setw(w) << inconclusive_mutations \
            << frac(inconclusive_mutations, mutations, "of total SNPs consensuses") << '\n';


        std::cout << "Total ambiguous consensuses count           " \
            << std::setw(w) << ambiguous_consensuses \
            << frac(ambiguous_consensuses, total_consensuses, "of total consensuses") << '\n';

        std::cout << "  Insertions count                          " \
            << std::setw(w) << ambiguous_consensuses_insertions \
            << frac(ambiguous_consensuses_insertions, total_consensuses_insertions, "of total insertions consensuses") << '\n';

        const auto ambiguous_mutations = \
            ambiguous_consensuses - ambiguous_consensuses_insertions;
        std::cout << "  SNPs count                                " \
            << std::setw(w) << ambiguous_mutations \
            << frac(ambiguous_mutations, mutations, "of total SNPs consensuses") << '\n';


        if (allele_variant_counter_vector.size() > 0) {
            const auto digits = (_max_allele_variant_counter > 0) ? \
                (int)log10((double)_max_allele_variant_counter) + 1 : 1;
            const auto& ctr = allele_variant_counter_vector;
            counter_t total = 0;
            for(size_t i=0; i != ctr.size(); ++i)
                total += ctr[i];

            const auto ado = ctr[1];
            std::cout << "Allelic dropout                             " \
                << std::setw(w) << ado \
                << frac(ado, total, "of total variants") << '\n';

            std::cout << '\n';
            std::cout << "  Allele count  Positions comprising variants\n";

            for(size_t i=1; i != ctr.size(); ++i) {
                const auto entry = ctr[i];

                std::cout << "    " << std::setw(digits) << i << "    " << std::setw(w) << entry \
                    << frac(entry, total, "of total variants") << '\n';
            }
            const auto entry = ctr[0];
            std::cout << " >= " << std::setw(digits) << _max_allele_variant_counter+1 \
                << "    " << std::setw(w) << entry \
                << frac(entry, total, "of total variants") << '\n';
        } else {
            std::cout << "Allelic dropout                             " \
                << std::setw(w) << 0 \
                << frac(0, 0, "of total variants") << '\n';
        }

/*
        if (temp_variable) {
            std::cout << "\nTemporary variable                    " \
            << temp_variable << '\n';
        }
*/
        std::cout << '\n';
    }

    std::string frac(const counter_t nom, const counter_t denom, const std::string& str) {

        constexpr int str_len = 100;
        char s[str_len+1];
        if (denom) {
            double f = 100. * (double)nom / (double)denom;
            snprintf(s, str_len, "  (%3.2f %% %s)", f, str.c_str());
        } else
            snprintf(s, str_len, "  ( - %% %s)", str.c_str());
        return std::string(s);
    }

    std::string header(const std::string& str) {
        auto len = str.size();
        std::string s1, s2;
        s1 = std::string(len+6, '*');
        s2 = "\n*  " + str + "  *\n";
        return std::string(s1 + s2 + s1);
    }

};

#endif // _ECSTATS_H_
