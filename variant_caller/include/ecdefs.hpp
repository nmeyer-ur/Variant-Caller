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

#ifndef _ECDEFS_H_
#define _ECDEFS_H_

#include <memory>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <bitset>


// ec version
#define _version "v0.12"

// use explicit SIMD for evaluation of Levenshtein distance
#define LEVEN_SIMD

// Levenshtein distance, evaluation of lower bound
#define LEVEN_LOWER_BOUND

// Levenshtein distance, early termination
#define LEVEN_EARLY_TERMINATION

// no threading in Eigen
#define EIGEN_DONT_PARALLELIZE

// super fast int to string conversion in tuple generation
#define FAST_INT_TO_STRING

// use hash table for TAG mapping
//#define TAG_CLUSTER_USE_HASH_TABLE

// late read data decoding after TAG clustering, potentially faster
//#define LATE_READ_DATA_DECODING


// SIMD
#define SIMD_WIDTH 1
#define SIMD_TYPE int

#ifdef LEVEN_SIMD

#ifdef __SSSE3__
#ifdef __SSE4_1__

#pragma message("SSE support")
#define USE_SSE

#include <tmmintrin.h>
#include <smmintrin.h>

#undef SIMD_WIDTH
#define SIMD_WIDTH 16

#undef SIMD_TYPE
#define SIMD_TYPE __m128i

#endif // SSE 4.1
#endif // SSSE 3

#ifdef __AVX__

#pragma message("AVX support")
#define USE_AVX

#include <tmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>

#undef SIMD_WIDTH
#define SIMD_WIDTH 32

#undef SIMD_TYPE
#define SIMD_TYPE __m256i

#endif // AVX

#else // LEVEN_SIMD

#pragma message("No explicit SIMD support for Levenshtein distance")

#endif // LEVEN_SIMD

typedef SIMD_TYPE simd_t;


// Statistics
class Statistics;

// Parameters
class Parameters;

// Allele
class Reference_base;
class Consensus_allele;
class Allele;

// Levenshtein distance
class Levenshtein;

// Cluster
class Cluster;

// Read
class Base_phred_pair;
class Insertion;
class Read;
class Read_preprocessing;

// Block
class Block;

// Variant calling
class Variant_caller;

// Alignment file
class Alignment_metadata;
class Alignment_file;

// Reference file
class Reference_file;

// DBSCAN
template<typename Vector, typename Matrix, typename Epsilon>
class DBSCAN;

// init values
constexpr size_t _init_block_vector_size = 200;
constexpr size_t _init_block_read_vector_size = 100;
constexpr size_t _init_variant_position_pool = 100;
constexpr size_t _num_prefetch_reads = 1;
constexpr size_t _init_prefetch_vector_size = 100;
constexpr int _section_gap = 10; //10;
constexpr size_t _init_output_string_size = 100000;
constexpr int _max_allele_variant_counter = 10;
constexpr int _max_allele_position_counter = 10;

// type definitions
typedef std::vector<Alignment_metadata> metadata_t;
typedef long long counter_t;
typedef int identity_t;
typedef size_t length_t;
typedef unsigned int position_t; // htslib: uint32_t:8
typedef unsigned char mapq_t; // htslib: uint32_t:8
typedef unsigned short flag_t; // htslib: uint32_t:16
typedef int coverage_t;
typedef std::vector<coverage_t> coverage_vector_t;
typedef char base_phred_elementary_t;
typedef base_phred_elementary_t base_t;
typedef base_phred_elementary_t phred_t;
typedef std::vector<base_t> base_vector_t;
typedef std::vector<phred_t> phred_vector_t;
typedef std::vector<Block> block_vector_t;
typedef std::shared_ptr<Read> read_ptr_t;
typedef std::vector<read_ptr_t> read_vector_t;
typedef long long tag_encoded_t;
typedef std::vector<tag_encoded_t> tag_encoded_vector_t;
typedef std::vector<counter_t> counter_vector_t;
typedef std::vector<Allele> allele_vector_t;
typedef std::vector<position_t> variant_position_vector_t;
typedef std::vector<Consensus_allele> consensus_pool_t;
typedef std::vector<consensus_pool_t> consensus_pool_vector_t;
typedef std::unordered_map<position_t, Insertion> insertions_map_t;
typedef std::unordered_map<position_t, std::vector<Insertion>> insertions_vector_map_t;
typedef std::array<int, SIMD_WIDTH> indices_array_t;
typedef std::array<char, SIMD_WIDTH> ld_array_t;
typedef std::map<tag_encoded_t, counter_t> tag_counter_map_t;
typedef std::bitset<4> orientation_t;
typedef std::bitset<4> tag_mask_t;

#ifdef TAG_CLUSTER_USE_HASH_TABLE
typedef std::unordered_map<tag_encoded_t, read_vector_t> cluster_map_t;
#else
typedef std::map<tag_encoded_t, read_vector_t> cluster_map_t;
#endif

// lookup definitions for bases conversion
constexpr int  _num_bases      = 6;
constexpr char _char_del_base  = 'D';
constexpr int  _bin_del_base   = 0x5;
constexpr char _char_A_base    = 'A';
constexpr int  _bin_A_base     = 0x1;
constexpr char _char_T_base    = 'T';
constexpr int  _bin_T_base     = 0x2;
constexpr char _char_C_base    = 'C';
constexpr int  _bin_C_base     = 0x3;
constexpr char _char_G_base    = 'G';
constexpr int  _bin_G_base     = 0x4;
constexpr char _char_N_base    = 'N';
constexpr int  _bin_N_base     = 0x0;
constexpr phred_t _phred_del_base = 50;

// TAG
constexpr int _tag_max_length = 16;
constexpr size_t _tag_min_length = 3;
constexpr char _tag_delim = '_';
#define _tag_na "NA"
constexpr int _tag_f1_index = 0;
constexpr int _tag_f2_index = 1;
constexpr int _tag_r1_index = 2;
constexpr int _tag_r2_index = 3;
constexpr int _tag_mask_options = 16;


// Levenshtein
typedef char leven_encoding_t;
typedef char leven_dist_t;
typedef std::array<char, _num_bases> tag_base_freq_array_t;

// DBSCAN / Eigen
typedef Eigen::Matrix<int, 1, Eigen::Dynamic> eigen_vector_t;
typedef Eigen::Matrix<leven_dist_t, Eigen::Dynamic, Eigen::Dynamic> distance_matrix_t;

// SAM Flags
constexpr int _samflag_index_read_paired = 0;
constexpr int _samflag_index_read_mapped_in_proper_pair = 1;
constexpr int _samflag_index_read_unmapped = 2;
constexpr int _samflag_index_mate_unmapped = 3;
constexpr int _samflag_index_read_reverse_strand = 4;
constexpr int _samflag_index_mate_reverse_strand = 5;
constexpr int _samflag_index_first_in_pair = 6;
constexpr int _samflag_index_second_in_pair = 7;
constexpr int _samflag_index_not_primary_alignment = 8;
constexpr int _samflag_index_read_fails_quality_check = 9;
constexpr int _samflag_index_read_is_pcr_or_optical_duplicate = 10;
constexpr int _samflag_index_supplementary_alignment = 11;

constexpr flag_t _samflag_mask_read_paired = 0x1;
constexpr flag_t _samflag_mask_read_mapped_in_proper_pair = 0x2;
constexpr flag_t _samflag_mask_read_unmapped = 0x4;
constexpr flag_t _samflag_mask_mate_unmapped = 0x8;
constexpr flag_t _samflag_mask_read_reverse_strand = 0x10;
constexpr flag_t _samflag_mask_mate_reverse_strand = 0x20;
constexpr flag_t _samflag_mask_first_in_pair = 0x40;
constexpr flag_t _samflag_mask_second_in_pair = 0x80;
constexpr flag_t _samflag_mask_not_primary_alignment = 0x100;
constexpr flag_t _samflag_mask_read_fails_quality_check = 0x200;
constexpr flag_t _samflag_mask_read_is_pcr_or_optical_duplicate = 0x400;
constexpr flag_t _samflag_mask_supplementary_alignment = 0x800;

// orientation
constexpr int _orientation_index_FR = 0;
constexpr int _orientation_index_RF = 1;
constexpr int _orientation_index_FF = 2;
constexpr int _orientation_index_RR = 3;

bool check_coverage(const read_vector_t&, const int);


static inline int char_base_to_bin_lookup(const char base) {

    int encoding = 0x77;
    switch(base) {
        case _char_A_base: encoding = _bin_A_base; break;
        case _char_T_base: encoding = _bin_T_base; break;
        case _char_C_base: encoding = _bin_C_base; break;
        case _char_G_base: encoding = _bin_G_base; break;
        case _char_N_base: encoding = _bin_N_base; break;
        case _char_del_base: encoding = _bin_del_base; break;
        default: break;
    }

    if (encoding == 0x77) {
        std::cout << "Unknown char base type: " << base << '\n';
        //assert(!"Unknown base type. Exiting.");
        encoding = _bin_N_base;
    }

    return encoding;
}

static inline base_t bin_base_to_char_lookup(const int base) {

    char encoding = 0x77;

    if      (base == _bin_A_base) encoding = _char_A_base;
    else if (base == _bin_T_base) encoding = _char_T_base;
    else if (base == _bin_C_base) encoding = _char_C_base;
    else if (base == _bin_G_base) encoding = _char_G_base;
    else if (base == _bin_N_base) encoding = _char_N_base;
    else if (base == _bin_del_base) encoding = _char_del_base;

    if (encoding == 0x77) {
        std::cout << "Unknown binary base type: " << base << '\n';
        //assert(!"Unknown base type. Exiting.");
        encoding = _char_N_base;
    }

    return encoding;
}

static inline std::string base_vector_to_string(const base_vector_t& v) {

    std::string s(v.begin(), v.end());
    return s;
}

static inline std::string tag_encoded_to_string(tag_encoded_t tag) {

    std::string s;
    for(auto i=0; i<_tag_max_length; ++i) {
        int base = (int)(tag & 0xf);
        s += (char)bin_base_to_char_lookup(base);
        tag = tag >> 4;
    }
    return s;
}

tag_encoded_t _tag_encoded_na = \
    ((((tag_encoded_t)char_base_to_bin_lookup(_char_N_base)) & 0xf) << 4) | \
    (((tag_encoded_t)char_base_to_bin_lookup(_char_A_base)) & 0xf);
//_tag_encoded_na = _tag_encoded_na << 4;
//_tag_encoded_na |= ((tag_encoded_t)char_base_to_bin_lookup(_char_A_base)) & 0xf;

class Alignment_metadata {

public:
    int index;
    unsigned long long num_reads;
    std::string name;
};

#endif // _ECDEFS_H_
