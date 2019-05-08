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

#ifndef _ECSAM_H_
#define _ECSAM_H_

//refined

#include <iostream>
#include <cassert>
#include <sstream>
#include "sam.h"
#include "ecdefs.hpp"


class Alignment_file {

public:
    Alignment_file(const Parameters& p): params(p), aln_file_handle(nullptr),
        _aln_header(nullptr), _aln(nullptr), _idx(nullptr), _itr(nullptr),
        _eof(false), aln_pending(false) {

        aln_file_handle = hts_open(params.alignment_file.c_str(), "r");

        if (params.hts_threads > 0) {
            const auto result = hts_set_threads(aln_file_handle, params.hts_threads);
            if (result)
                assert(!"HTSlib threading error. Exiting.");
        }

        _aln_header = sam_hdr_read(aln_file_handle);
        if (_aln_header == nullptr)
            assert(!"Unable to open alignment file. Exiting.");

        _idx = sam_index_load(aln_file_handle, params.alignment_file.c_str());
        if (_aln_header == nullptr)
            assert(!"Unable to open alignment index file. Exiting.");
    }

    ~Alignment_file() {

        if (_itr != nullptr)
            hts_itr_destroy(_itr);
        if (_idx != nullptr)
            hts_idx_destroy(_idx);
        if (_aln != nullptr)
            bam_destroy1(_aln);
        if (_aln_header != nullptr)
            bam_hdr_destroy(_aln_header);
        if (aln_file_handle != nullptr)
            hts_close(aln_file_handle);
    }

    // define iterator
    void select_region(const std::string& region) {

        _itr = sam_itr_querys(_idx, _aln_header, region.c_str());
        aln_pending = false;
        _eof = false;
    }

    // fetch reads until parameter settings are valid or eof
    bool fetch(void) {

        if ((aln_pending == false) && (_eof == false)) {
            _aln = bam_init1();
            //if (sam_read1(aln_file_handle, _aln_header, _aln) >= 0) {
            if (sam_itr_next(aln_file_handle, _itr, _aln) >= 0) {
                statistics.total_reads_in_alignment_file += 1;
                aln_pending = true;
            } else {
                aln_pending = false;
                _eof = true;
                bam_destroy1(_aln);
                _aln = nullptr;
            }
        }
        return aln_pending;
    }

    metadata_t get_metadata() {

        auto entries = _aln_header->n_targets;
        int output_index = 0;

        for(decltype(entries) i = 0; i != entries; ++i) {
            unsigned long long mapped, unmapped;
            hts_idx_get_stat(_idx, i, (uint64_t*)&mapped, (uint64_t*)&unmapped);
            const auto num_reads = mapped;
            if (num_reads >= (unsigned long long)params.min_reads) {
                Alignment_metadata m;
                m.index = output_index;
                m.name = _aln_header->target_name[i];
                m.num_reads = num_reads;

                // region check
                if (!params.region_check) {
                    _metadata.push_back(m);
                    ++output_index;
                } else {
                    const auto& rv = params.region_vector;
                    if (std::find(rv.begin(), rv.end(), m.name) != rv.end()) {
                        _metadata.push_back(m);
                        ++output_index;
                    }
                }

/* select all regions, verified & working
                _metadata.push_back(m);
                ++output_index;
*/
            }
        }

        // sort by number of reads in region
        std::sort(_metadata.begin(), _metadata.end(), \
            [](Alignment_metadata a, Alignment_metadata b){ return a.num_reads > b.num_reads; });
        return _metadata;
    }

    void update_statistics(Statistics& ext_statistics) {

        ext_statistics.total_reads_in_alignment_file += statistics.total_reads_in_alignment_file;
    }

    std::string get_region_name_by_index(const int index) {

        return std::string(_aln_header->target_name[index]);
    }

    int get_region_len_by_index(const int index) {

        return _aln_header->target_len[index];
    }

    unsigned long long get_region_read_count_by_index(const int index) {

        unsigned long long mapped, unmapped;
        hts_idx_get_stat(_idx, index, (uint64_t*)&mapped, (uint64_t*)&unmapped);
        return (unsigned long long)(mapped + unmapped);
    }

    bam_hdr_t*& aln_header(void) { return _aln_header; }
    bam1_t*& aln(void) { return _aln; }
    void accept_aln(void) { aln_pending = false; }
    bool eof(void) { return _eof; }
    position_t read_start(void) { return (position_t)_aln->core.pos; }
    position_t read_end(void) { return (position_t)_aln->core.pos + (position_t)_aln->core.l_qseq;}
    int get_num_regions(void) { return _aln_header->n_targets; }

private:
    Statistics statistics;
    const Parameters& params;
    metadata_t _metadata;
    samFile* aln_file_handle;
    bam_hdr_t* _aln_header;
    bam1_t* _aln;
    hts_idx_t* _idx;
    hts_itr_t* _itr;
    bool _eof;
    bool aln_pending;
};


#endif // _ECSAM_H_
