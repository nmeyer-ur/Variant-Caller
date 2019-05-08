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

#ifndef _ECREADS_H_
#define _ECREADS_H_

#include <iomanip>
#include <iostream>
#include <bitset>
#include <cassert>
#include <vector>
#include "sam.h"
#include "ecdefs.hpp"


class Insertion {

public:
    base_vector_t base;
    phred_vector_t phred;
};


class Read {

public:
    static identity_t read_id_ctr;
    identity_t id;
    position_t start, end;
    length_t length;
    mapq_t mapq;
    base_vector_t base;
    phred_vector_t phred;
    const tag_encoded_t tag_encoded;
    insertions_map_t insertions;
	orientation_t orientation;
	bool is_forward;

    Read(const Parameters& p, bam1_t* aln_ptr,
        const bam_hdr_t* aln_header_ptr, const tag_encoded_t& tag,
		const orientation_t _orientation, const bool _is_forward):
        length(0), tag_encoded(tag),
		orientation(_orientation), is_forward(_is_forward),
		params(p), aln_header(aln_header_ptr), aln(aln_ptr) {

        ++read_id_ctr;
        id = read_id_ctr;
        _process_read1();
#ifndef LATE_READ_DATA_DECODING
        process_read2();
#endif
    }

    ~Read() {
        if (aln != nullptr)
            bam_destroy1(aln);
    }

    void process_read2(void) {
        _process_read2();
    }

    void display_stats(void) {

        std::cout << "Read id:        " << id << '\n';
        std::cout << "TAG:            " << tag_encoded_to_string(tag_encoded) << '\n';
        std::cout << "Start position: " << start << '\n';
        std::cout << "End position:   " << end << '\n';
        std::cout << "Length:         " << length << '\n';
        std::cout << "MAPQ:           " << mapq << '\n';
    }

    void display_read(void) {

        const auto b = base_vector_to_string(base);
        std::cout << start << ' ' << id << " 0x" << \
            ' ' << tag_encoded_to_string(tag_encoded) << \
            ' ' << std::dec << b << '\n';
    }

private:
    const Parameters& params;
    const bam_hdr_t* aln_header;
    bam1_t* aln;

    void _process_read1(void) {

        mapq  = (mapq_t)aln->core.qual;
        start = (position_t)aln->core.pos;

        // calculate the length of the read sequence
        // removing clipping and insertions, but taking
        // into account deletions
        auto n_cigar = aln->core.n_cigar;
        auto* cigar  = bam_get_cigar(aln);

        for (decltype(n_cigar) i=0; i<n_cigar; ++i) {

            auto op = cigar[i] & 0xf; // get operation
            auto l  = cigar[i] >> 4;  // get operation count

            switch(op) {
                // mis/match
                case 0: length += l; break;
                // deletion
                case 2: length += l; break;
                default: break;
            }
        }
        end = start + length;
    }

    void _process_read2(void) {

        base.resize(length);
        phred.resize(length);

        // Base / phred encoding
        auto index_aln = 0; // reencoded sequence used for variant calling
        auto index_seq = 0; // original sequence from file

        auto seq_phred = bam_get_qual(aln);
        auto* seq_base = bam_get_seq(aln);

        const auto phred_thres = params.phred_score_threshold;

        auto n_cigar = aln->core.n_cigar;
        const auto* cigar  = bam_get_cigar(aln);

        for(decltype(n_cigar) i=0; i!=n_cigar; ++i) {

            const auto op = cigar[i] & 0xf;  // get operation
            auto l  = cigar[i] >> 4;    // get operation count

            switch(op) {
                case 0: // mis/match
                    for(decltype(l) j=0; j<l; ++j) {
                        phred_t p = (phred_t)seq_phred[index_seq+j];
                        base_t c = _char_N_base;
                        if (p >= phred_thres)
                            c = (base_t)seq_nt16_str[bam_seqi(seq_base, index_seq+j)];
                        base[index_aln+j] = c;
                        phred[index_aln+j] = p;
                    }
                    index_seq += l;
                    index_aln += l;
                    break;
                case 1: // insertion
                    {
                        Insertion ins;
                        ins.base.reserve(l);
                        ins.phred.reserve(l);
                        for(decltype(l) j=0; j<l; ++j) {
                            phred_t p = (phred_t)seq_phred[index_seq+j];
                            base_t c = _char_N_base;
                            if (p >= phred_thres)
                                c = (base_t)seq_nt16_str[bam_seqi(seq_base, index_seq+j)];
                            ins.base.push_back(c);
                            ins.phred.push_back(p);
/*
                            std::cout << "Read id " << read_id_ctr \
                                << " insertion " << j+1 \
                                << " at " << start + index_seq \
                                << " base " << c << '\n';
*/
                        }
                        insertions[start + index_aln] = ins;

                    }
                    index_seq += l;
                    break;
                case 2: // deletion
                    for(decltype(l) j=0; j<l; ++j) {
                        phred_t p = _phred_del_base;
                        phred[index_aln+j] = p;
                        base_t c = _char_del_base;
                        base[index_aln+j] = c;
                    }
                    index_aln += l;
                    break;
                case 4: // soft clip
                    index_seq += l;
                    break;
                case 5: // hard clip
                    index_seq += l;
                    break;
                default:
                    assert(!"Unknown cigar encoding type. Exiting.");
            }
        }

        // dealloc aln
        bam_destroy1(aln);
        aln = nullptr;
    }
};

identity_t Read::read_id_ctr = 0;



class Read_preprocessing {

public:
    std::string region;
    tag_encoded_t tag_encoded;
	orientation_t orientation;
	bool is_forward;

    Read_preprocessing(const Parameters& p, const bam_hdr_t* aln_header_ptr):
        params(p), _aln_header(aln_header_ptr) {}

    bool is_valid(bam1_t* aln) {

        _aln = aln;
        valid = _preprocess();
        if ((valid == false) && (_aln != nullptr))
            bam_destroy1(_aln);
        return valid;
    }

private:
    Statistics statistics;
    const Parameters& params;
    const bam_hdr_t* _aln_header;
    bam1_t* _aln;
    bool valid;

    bool _preprocess(void) {

        bool valid = true;

        if (_aln == nullptr)
            return false;

/*
        if (params.region_check) {
            // segfault on target name if not checked!
            if (_aln->core.tid < 0)
                return false;

            // get region
            auto* chr = _aln_header->target_name[_aln->core.tid];
            std::string region = chr;
        }
*/

        // check MAPQ
        const mapq_t mapq = (mapq_t)_aln->core.qual;
        if (mapq < params.min_mapq)
            valid = false;

        // check samflags
        const auto flags = (flag_t)_aln->core.flag;
        const std::bitset<12> bitflags(flags);
        if (flags & params.flag_mask)
            valid = false;

        // forward / reverse filter
		is_forward = !bitflags.test(_samflag_index_read_reverse_strand);

        if (!params.forward && is_forward)
            valid = false;

        if (!params.reverse && !is_forward)
            valid = false;

        // check pair state
        if (!(bitflags.test(_samflag_index_read_paired) && \
                bitflags.test(_samflag_index_read_mapped_in_proper_pair))) {
            //std::cout << "not a pair" << '\n';
            valid = false;
        }

        // terminate
        if (params.verbosity >= 5) {
            std::cout << '\n';
            std::cout << "QUERY " << bam_get_qname(_aln) << '\n';
            std::cout << "FLAGS 0b"<< bitflags << "  0x" << std::hex << flags \
                << std::dec << "  " << flags << '\n';
        }

        if (!valid) {
            if (params.verbosity >= 5)
                std::cout << "Invalid\n";
            return false;
        }

        // read TAGs
        std::string tag[4];
        auto *s = bam_get_qname(_aln);
        auto i = 0;
        auto token_index = 0;
        auto tag_index = 0;
        while(*(s+i) != '\0') {
            if (*(s+i) == _tag_delim) {
                ++token_index;
                tag_index = token_index / 2 - 1;
                if (tag_index >= 4) break;
//                std::cout << "token index = " << token_index << "   tag index = " << tag_index << '\n';
//                if (token_index > params.tag_index)
//                    break;
            } else if (!(token_index & 1) && (token_index > 1)) {
                tag[tag_index] += *(s+i);
            }
            ++i;
        }

        if (params.verbosity >= 5) {
            std::cout << "TAG F1  " << tag[_tag_f1_index] << '\n';
            std::cout << "TAG F2  " << tag[_tag_f2_index] << '\n';
            std::cout << "TAG R1  " << tag[_tag_r1_index] << '\n';
            std::cout << "TAG R2  " << tag[_tag_r2_index] << '\n';
        }

		orientation.reset();

        // case: FR
        if ((bitflags.test(_samflag_index_first_in_pair) && \
            !bitflags.test(_samflag_index_read_reverse_strand) && \
            bitflags.test(_samflag_index_mate_reverse_strand)) || \
            (bitflags.test(_samflag_index_second_in_pair) && \
            bitflags.test(_samflag_index_read_reverse_strand) && \
            !bitflags.test(_samflag_index_mate_reverse_strand)))
        {
			orientation.set(_orientation_index_FR);
        }
        // case: RF
        else if ((bitflags.test(_samflag_index_first_in_pair) && \
            bitflags.test(_samflag_index_read_reverse_strand) && \
            !bitflags.test(_samflag_index_mate_reverse_strand)) || \
            (bitflags.test(_samflag_index_second_in_pair) && \
            !bitflags.test(_samflag_index_read_reverse_strand) && \
            bitflags.test(_samflag_index_mate_reverse_strand)))
        {
			orientation.set(_orientation_index_RF);
            std::swap(tag[_tag_f1_index], tag[_tag_r1_index]);
            std::swap(tag[_tag_f2_index], tag[_tag_r2_index]);
        }
        // case FF
        else if (!bitflags.test(_samflag_index_read_reverse_strand) && \
            !bitflags.test(_samflag_index_mate_reverse_strand)) \
        {
			orientation.set(_orientation_index_FF);
        }
        // case RR
        else if (bitflags.test(_samflag_index_read_reverse_strand) && \
            bitflags.test(_samflag_index_mate_reverse_strand)) \
        {
			orientation.set(_orientation_index_RR);
        }
        else
        {
            valid = false;
        }

		if ((orientation & params.orientation) == 0)
			valid = false;

        tag_mask_t tag_mask(0);
        tag_mask.set(_tag_f1_index, tag[_tag_f1_index].size() >= _tag_min_length);
        tag_mask.set(_tag_f2_index, tag[_tag_f2_index].size() >= _tag_min_length);
        tag_mask.set(_tag_r1_index, tag[_tag_r1_index].size() >= _tag_min_length);
        tag_mask.set(_tag_r2_index, tag[_tag_r2_index].size() >= _tag_min_length);

/*
        if (tag_mask == 0b0001 || tag_mask == 0b0100)
        {}
        else    valid = false;
*/

        if (tag_mask == 0b0000)
            valid = false;
        else if (tag_mask == 0b0010)
            valid = false;
        else if (tag_mask == 0b0011)
            valid = false;
        else if (tag_mask == 0b0111)
            valid = false;
        else if (tag_mask == 0b1000)
            valid = false;
        else if (tag_mask == 0b1010)
            valid = false;
        else if (tag_mask == 0b1011)
            valid = false;
        else if (tag_mask == 0b1100)
            valid = false;
        else if (tag_mask == 0b1101)
            valid = false;
        else if (tag_mask == 0b1110)
            valid = false;

        auto tag_selected_index = -1;
        for(size_t i=0; i != params.tag_mask.size(); ++i) {
            if (params.tag_mask.test(i) && tag_mask.test(i)) {
                tag_selected_index = i;
                break;
            }
        }

        if (tag_selected_index == -1)
            valid = false;

        if (params.verbosity >= 5) {
            std::cout << "Orientation ";
            if (orientation.test(_orientation_index_FR)) std::cout << "FR\n";
            else if (orientation.test(_orientation_index_RF)) std::cout << "RF\n";
            else if (orientation.test(_orientation_index_FF)) std::cout << "FF\n";
            else if (orientation.test(_orientation_index_RR)) std::cout << "RR\n";
            else if (orientation == 0) std::cout << "Error\n";
            std::cout << "- rearrange -\n";
            std::cout << "TAG F1  " << tag[_tag_f1_index] << '\n';
            std::cout << "TAG F2  " << tag[_tag_f2_index] << '\n';
            std::cout << "TAG R1  " << tag[_tag_r1_index] << '\n';
            std::cout << "TAG R2  " << tag[_tag_r2_index] << '\n';
            std::cout << "TAG Mask 0b" << tag_mask << '\n';
            std::cout << "TAG Index " << tag_selected_index << '\n';
            if (valid) {
                std::cout << "TAG " << tag[tag_selected_index] << '\n';
                std::cout << "Valid\n";
            } else
                std::cout << "Invalid\n";
        }

        if (!valid)
            return false;

        tag_encoded = 0;
        std::string& tag_selected = tag[tag_selected_index];
        for(size_t i=0; i != tag_selected.size(); ++i) {
            tag_encoded = tag_encoded << 4;
            tag_encoded |= (tag_encoded_t)char_base_to_bin_lookup(tag_selected[i]) & 0xf;
        }

        statistics.valid_reads_in_alignment_file += 1;
        return true;
    }
};


#endif // _ECREADS_H_
