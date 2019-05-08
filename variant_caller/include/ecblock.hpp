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

#ifndef _ECBLOCK_H_
#define _ECBLOCK_H_

#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include "ecdefs.hpp"


class Block {

public:
    identity_t block_id_ctr;
    std::string region;
    position_t start, end;
    length_t length;
    read_vector_t read_vector;
    int max_depth;

    Block(const Parameters& p, const std::string _region):
        region(_region), start(0), end(0), length(0), max_depth(0), params(p) {

        block_id_ctr = 1;
        read_vector.reserve(_init_block_read_vector_size);
    }

    bool append_read(const read_ptr_t& read) {

        if (read_vector.size() == 0) {
            start  = read->start;
            end    = read->end;
            length = end - start;
        } else if (read->start >= end + _section_gap) {
            return false;
        } else {
            start  = std::min(start, read->start);
            end    = std::max(end, read->end);
            length = end - start;
        }

        read_vector.push_back(std::move(read));
        return true;
    }

    void display_stats(void) {

        std::cout << "Block " << std::setw(5)  << block_id_ctr << " : " << region << ", " \
            << "0-base pos "  << std::setw(12) << start << " - " << std::setw(12) << end << ", " \
            << "bases count " << std::setw(5)  << length << ", " \
            << "read count "  << std::setw(5)  << read_vector.size() << ", " \
            << "max. depth "  << std::setw(5)  << max_depth << '\n';
    }

    // clear block
    void clear(void) {

        start  = 0;
        end    = 0;
        length = 0;
        read_vector.clear();
        //region = "none";
        max_depth = 0;
        ++block_id_ctr;
    }

    // early check of acceptance
    bool read_count_is_sufficient(void) {

        return (read_vector.size() >= (size_t)params.min_reads);
    }

    // check if block meets requirements for acceptance
    bool is_valid(void) {

        // first check the number of reads
        if (read_vector.size() < (size_t)params.min_reads)
            return false;

        // next evaluate the maximum coverage
        // build step counter
        coverage_vector_t delta(length+1, 0);

        for(const auto& read : read_vector) {
            auto s = read->start - start;
            auto e = read->end - start;
            delta[s] += 1;
            delta[e] -= 1;
        }

        // evaluate maximum coverage
        auto cum_sum = 0;
        for(const auto& step : delta) {
            cum_sum += step;
            max_depth = std::max(max_depth, cum_sum);
        }

        if (max_depth < params.min_reads)
            return false;

        // block is valid
        return true;
    }

    void display_reads(void) {

        std::cout << "Block reads" << '\n';
        for(const auto& read : read_vector) {
            const auto delta = read->start - start;
            auto aln = base_vector_to_string(read->base);
            // TODO there might be a better way to insert spaces
            aln.insert(0, delta, ' ');
            std::cout << read->start << ' ' << read->end << ' ' \
                << read->id << ' ' << aln << '\n';
        }
    }

    void update_statistics(Statistics& ext_statistics) {

        ext_statistics.valid_blocks += 1;
        ext_statistics.valid_reads_in_valid_blocks += read_vector.size();
        ext_statistics.valid_block_length_before_clustering_vector.push_back((counter_t)length);
    }

private:
    Statistics statistics;
    const Parameters& params;
};

// check if depth of coverage fulfills minimum requirements
// true if yes, false if no
bool check_coverage(const read_vector_t& read_vector, const int min_coverage) {

    bool result = false;

    // have to check coverage within allele at this step
    if (read_vector.size() < (size_t)min_coverage)
        return result;

    // evaluate allele boundaries
    auto start = read_vector[0]->start;
    position_t end = 0;
    for(const auto& read : read_vector) {
        auto s = read->start;
        auto e = read->end;
        start = std::min(s, start);
        end   = std::max(e, end);
    }

    auto length = end - start;

    // evaluate coverage
    coverage_vector_t delta(length+1, 0);
    for(const auto& read : read_vector) {
        auto s = read->start - start;
        auto e = read->end - start;
        delta[s] += 1;
        delta[e] -= 1;
    }

    auto cum_sum = 0;
    for(auto it=delta.begin(); it != --delta.end(); ++it) {
        cum_sum += *it;
        if (cum_sum >= min_coverage)
            return true;
    }
    return false;
}

#endif // _ECBLOCK_H_
