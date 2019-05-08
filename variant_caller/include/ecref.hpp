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

#ifndef _ECREF_H_
#define _ECREF_H_


#include <iostream>
#include <cassert>
#include "sam.h"
#include "faidx.h"
#include "ecdefs.hpp"


class Reference_file {

public:
    Reference_file(const Parameters& p): params(p), ref_file_handle(nullptr) {

        ref_file_handle = fai_load(params.ref_fasta_file.c_str());
        if (!ref_file_handle)
            assert(!"Error opening reference. Exiting.");

    }

    ~Reference_file() { fai_destroy(ref_file_handle); }

    void fetch(base_vector_t& reference, const std::string& region,
        const position_t& start, const position_t& end) {

        int length;
        char* ref_ptr = \
            faidx_fetch_seq(ref_file_handle, region.c_str(), start, end, &length);

        if (!ref_ptr)
            assert(!"Error loading reference. Exiting.");

        reference.reserve(end - start);
        for(auto index=0; index<length-1; ++index)
            reference.push_back(toupper(ref_ptr[index]));

        // spec says free the reference sequence manually
        free(ref_ptr);

        statistics.reference_bases_loaded_from_file += length-1;
    }

    void update_statistics(Statistics& ext_statistics) {

        ext_statistics.reference_bases_loaded_from_file += statistics.reference_bases_loaded_from_file;
    }

private:
    Statistics statistics;
    const Parameters& params;
    faidx_t* ref_file_handle;

};


#endif // _ECREF_H_
