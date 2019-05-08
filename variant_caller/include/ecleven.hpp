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

#ifndef _ECLEVEN_H_
#define _ECLEVEN_H_

#include "ecdefs.hpp"

//counter_t tag_base_freq(tag_encoded_t);
//leven_dist_t leven_lower_bound_by_base_freq(counter_t, counter_t);
//int leven_gen(tag_encoded_t, tag_encoded_t, int);



#ifdef BLA // moved into class
/*  TAG bases frequencies encoding

    Microbenchmarks (g++ 5.4/6.x -O3 -march=native)

    Intel Core2 Duo P8700 @ 2.53GHz             36 cycles
    Intel Xeon E5645 Westmere @ 2.40GHz         30 cycles
    Intel Xeon E5-2620 v4 Broadwell @ 2.10GHz   17 cycles

    faster than any explicit SSE implementation
*/
inline counter_t tag_base_freq(tag_encoded_t u) {

    counter_t result = 0;

    // fast scalar implementation, 16 * 5 = 80 elementary Operations (w/o loop)
    const counter_t adder = 1;
    for(int i=0; i != _tag_max_length; i++) {
        counter_t shift = (u & 0xF);
        shift = shift << 3;
        result += (adder << shift);
        u = u >> 4;
    }

    return result;
}

/*  Levenshtein distance lower bound by bases frequencies evaluation

    Microbenchmarks (g++ 5.4/6.x -O3 -march=native)

    Intel Core2 Duo P8700 @ 2.53GHz             10 cycles
    Intel Xeon E5645 Westmere @ 2.40GHz          9 cycles
    Intel Xeon E5-2620 v4 Broadwell @ 2.10GHz    7 cycles

    faster than any explicit SSE implementation
*/
inline leven_dist_t leven_lower_bound_by_base_freq(counter_t u, counter_t v) {

    counter_t result = 0;
    const counter_t mask = 0xFF;

    for(int i=0; i != _num_bases-1; ++i) {
        result += abs((u & mask) - (v & mask));
	    u = u >> 8;
	    v = v >> 8;
    }

    return (leven_dist_t)(result >> 1);
}

/*  Hamming distance between two TAGs
    corresponding to minimum number of insertions/deletions

    not recommended
*/
inline leven_dist_t tag_hamming(tag_encoded_t u, tag_encoded_t v) {

    counter_t result = 0;
    counter_t mask = 0xF;
	const counter_t x = u ^ v;

    for(int i=0; i != _tag_max_length; ++i) {
		result += ((x & mask) != 0);
		mask = mask << 4;
    }
    return (leven_dist_t)result;
}
#endif // BLA

/*  Levenshtein distance, generic version
    slow, but allows for early termination and thus is fast

    Microbenchmarks (g++ 5.4/6.x -O3 -march=native)

    Intel Core2 Duo P8700 @ 2.53GHz             1930 cycles [1 distance; SSE]
    Intel Xeon E5645 Westmere @ 2.40GHz         2190 cycles [1 distance; SSE]
    Intel Xeon E5-2620 v4 Broadwell @ 2.10GHz    860 cycles [1 distance; SSE]
*/
inline int leven_gen(tag_encoded_t a, tag_encoded_t b, int et) {

    std::array<int, _tag_max_length> v, s, t;
    int ret = 0;

    int length_a = 0, length_b = 0;

    for(int i=0; i<_tag_max_length; ++i) {
        v[i] = i+1;
        s[i] = (int)(a & 0xf);
        if (s[i] != _bin_N_base) ++length_a;
        a = a >> 4;
        t[i] = (int)(b & 0xf);
        b = b >> 4;
        if (t[i] != _bin_N_base) ++length_b;
    }

    int length = std::max(length_a, length_b);

    for(int i=0; i != length; ++i) {

        int tmp = i;
        ret = i + 1;

        for (int j=0; j != length; ++j) {

            int tmp2 = (t[i] == s[j]) ? tmp : tmp + 1;
            tmp = v[j];
            // amazingly ugly and confusing
            ret = v[j] = tmp > ret ? tmp2 > ret ? ret + 1 : tmp2 : tmp2 > tmp ? tmp + 1 : tmp2;
        }

#ifdef LEVEN_EARLY_TERMINATION // early termination
        if (v[i] > et) return et+1; // v[i];
#endif
    }
    return ret;
}

#ifdef USE_SSE
/*  transpose 1x8 byte (1x tag_encoded_t) and
    decode/decompress into 16x16 byte array

    Microbenchmarks (g++ 5.4/6.x -O3 -march=native)

    Intel Core2 Duo P8700 @ 2.53GHz               - cycles
    Intel Xeon E5645 Westmere @ 2.40GHz           - cycles
    Intel Xeon E5-2620 v4 Broadwell @ 2.10GHz     - cycles
*/
inline void transpose_decode_1x8_to_16x16_sse( \
    const tag_encoded_t& tag, \
    __m128i* B) {

    tag_encoded_t u[2] = {tag >> 4, tag};
    __m128i F = _mm_load_si128(reinterpret_cast<__m128i*>(&u));
    __m128i mask = _mm_set1_epi8(0xF);
    F = _mm_and_si128(F, mask);

    for(int i=0; i != 16; ++i) {
        __m128i decodemask = _mm_set1_epi8(i);
        B[i] = _mm_shuffle_epi8(F, decodemask);
    }
}


/*  transpose 16x8 byte (16x tag_encoded_t) into 8x16 byte array,
    then decode/decompress into 16x16 byte array

    Microbenchmarks (g++ 5.4/6.x -O3 -march=native)

    Intel Core2 Duo P8700 @ 2.53GHz              98 cycles
    Intel Xeon E5645 Westmere @ 2.40GHz          44 cycles
    Intel Xeon E5-2620 v4 Broadwell @ 2.10GHz    28 cycles
*/
inline void transpose_decode_16x8_to_16x16_sse( \
    const tag_encoded_vector_t& tag, \
    const indices_array_t& index_array, \
    const int offset, \
    __m128i* B) {

    __m128i E[2], F[4], T[8]; // temporaries
    const __m128i shuffle4x4 = _mm_set_epi8(15,11,7,3, 14,10,6,2, 13,9,5,1, 12,8,4,0);
    const __m128i decodemask = _mm_set1_epi8(0xF);

#define load(i) \
    F[0] = _mm_set_epi64x(tag[index_array[(i)+1+offset]], tag[index_array[(i)+0+offset]]); \
    F[1] = _mm_set_epi64x(tag[index_array[(i)+3+offset]], tag[index_array[(i)+2+offset]]); \
    F[2] = _mm_set_epi64x(tag[index_array[(i)+5+offset]], tag[index_array[(i)+4+offset]]); \
    F[3] = _mm_set_epi64x(tag[index_array[(i)+7+offset]], tag[index_array[(i)+6+offset]])

#define cast(dst) \
    T[(dst)+0] = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(F[0]), _mm_castsi128_ps(F[1]), 0b10001000)); \
    T[(dst)+1] = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(F[2]), _mm_castsi128_ps(F[3]), 0b10001000)); \
    T[(dst)+2] = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(F[0]), _mm_castsi128_ps(F[1]), 0b11011101)); \
    T[(dst)+3] = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(F[2]), _mm_castsi128_ps(F[3]), 0b11011101))

#define shuffle(src) \
    F[0] = _mm_shuffle_epi8(T[(src)+0], shuffle4x4); \
    F[1] = _mm_shuffle_epi8(T[(src)+1], shuffle4x4); \
    F[2] = _mm_shuffle_epi8(T[(src)+2], shuffle4x4); \
    F[3] = _mm_shuffle_epi8(T[(src)+3], shuffle4x4)

#define unpack(dst) \
    T[(dst)+0] = _mm_unpacklo_epi32(F[0], F[1]); \
    T[(dst)+1] = _mm_unpackhi_epi32(F[0], F[1]); \
    T[(dst)+2] = _mm_unpacklo_epi32(F[2], F[3]); \
    T[(dst)+3] = _mm_unpackhi_epi32(F[2], F[3])

#define decode(src1, src2, dst) \
    E[0] = _mm_unpacklo_epi64(T[(src1)], T[(src2)]); \
    B[(dst)+0] = _mm_and_si128(E[0], decodemask); \
    B[(dst)+1] = _mm_srli_epi64(E[0], 4); \
    B[(dst)+1] = _mm_and_si128(B[(dst)+1], decodemask); \
    E[1] = _mm_unpackhi_epi64(T[(src1)], T[(src2)]); \
    B[(dst)+2] = _mm_and_si128(E[1], decodemask); \
    B[(dst)+3] = _mm_srli_epi64(E[1], 4); \
    B[(dst)+3] = _mm_and_si128(B[(dst)+3], decodemask)

    // first 8x8

    load(0);
    cast(0);
    shuffle(0);
    unpack(0);

    // second 8x8

    load(8);
    cast(4);
    shuffle(4);
    unpack(4);

    // unpack and decode

    decode(0, 4,  0);
    decode(1, 5,  4);
    decode(2, 6,  8);
    decode(3, 7, 12);

#undef load
#undef cast
#undef shuffle
#undef unpack
#undef decode
}

/*  Levenshtein distance SSE/AVX implementation for 16/32 distances in parallel

    Microbenchmarks (g++ 5.4/6.x -O3 -march=native)

    Intel Core2 Duo P8700 @ 2.53GHz             3850 cycles [16 distances; SSE]
    Intel Xeon E5645 Westmere @ 2.40GHz         3440 cycles [16 distances; SSE]
    Intel Xeon E5-2620 v4 Broadwell @ 2.10GHz   1525 cycles [16 distances; SSE]
*/
inline ld_array_t leven_sse( \
    const tag_encoded_vector_t& tag_encoded_vector, \
    const indices_array_t& row_indices, \
    const indices_array_t& column_indices, \
    int et) {

// SIMD intrinsics definitions
#ifndef USE_AVX
#ifdef USE_SSE
#define simd_set1(a)        _mm_set1_epi8((a))
#define simd_setzero()      _mm_setzero_si128()
#define simd_add(a, b)      _mm_add_epi8((a), (b))
#define simd_cmpeq(a, b)    _mm_cmpeq_epi8((a), (b))
#define simd_min(a, b)      _mm_min_epu8((a), (b))
#define simd_cmpgt(a, b)    _mm_cmpgt_epi8((a), (b))
#define simd_movemask(a)    _mm_movemask_epi8((a))
#define simd_store(a, b)    _mm_store_si128((a), (b))
#endif
#endif
#ifdef USE_AVX
#define simd_set1(a)        _mm256_set1_epi8((a))
#define simd_setzero()      _mm256_setzero_si256()
#define simd_add(a, b)      _mm256_add_epi8((a), (b))
#define simd_cmpeq(a, b)    _mm256_cmpeq_epi8((a), (b))
#define simd_min(a, b)      _mm256_min_epu8((a), (b))
#define simd_cmpgt(a, b)    _mm256_cmpgt_epi8((a), (b))
#define simd_movemask(a)    _mm256_movemask_epi8((a))
#define simd_store(a, b)    _mm256_store_si256((a), (b))
#endif

    ld_array_t result;

	simd_t a[_tag_max_length];
	simd_t b[_tag_max_length];

#ifndef USE_AVX
#ifdef USE_SSE
    constexpr int _et_res_cmp = 0xFFFF; // early termination mask
    transpose_decode_16x8_to_16x16_sse( \
        tag_encoded_vector, row_indices, 0, a);
    transpose_decode_16x8_to_16x16_sse( \
        tag_encoded_vector, column_indices, 0, b);
#endif
#endif
#ifdef USE_AVX
    constexpr int _et_res_cmp = 0xFFFFFFFF; // early termination mask
    {
    __m128i a1[_tag_max_length], a2[_tag_max_length];
    __m128i b1[_tag_max_length], b2[_tag_max_length];

    transpose_decode_16x8_to_16x16_sse( \
        tag_encoded_vector, row_indices, 0, a1);
    transpose_decode_16x8_to_16x16_sse( \
        tag_encoded_vector, row_indices, 16, a2);
    for(int i=0; i != _tag_max_length; ++i) {
        const __m256i ai = _mm256_castsi128_si256(a1[i]);
        a[i] = _mm256_insertf128_si256(ai, a2[i], 1);
    }
    transpose_decode_16x8_to_16x16_sse( \
        tag_encoded_vector, column_indices, 0, b1);
    transpose_decode_16x8_to_16x16_sse( \
        tag_encoded_vector, column_indices, 16, b2);
    for(int i=0; i != _tag_max_length; ++i) {
        const __m256i bi = _mm256_castsi128_si256(b1[i]);
        b[i] = _mm256_insertf128_si256(bi, b2[i], 1);
    }
    }
#endif

    const simd_t one  = simd_set1(1);
    const simd_t et_v = simd_set1(et + 1);

   	simd_t diag [_tag_max_length + 1];
   	simd_t diag2[_tag_max_length + 1];

	diag2[0] = simd_setzero();

  	int i, j, k;

  	for (k = 1 ;; ++k) {

	    int startRow = k > _tag_max_length ? k - _tag_max_length : 1;
		int endRow =   k > _tag_max_length ? _tag_max_length : k - 1;

		for (i = endRow; i >= startRow; --i) {

			j = k - i;

            simd_t above = simd_add(diag2[i-1], one);
            simd_t left  = simd_add(diag2[i], one);

            simd_t delta = simd_cmpeq(a[i-1], b[j-1]);
            delta = simd_add(delta, one);
            simd_t left_above = simd_add(diag[i-1], delta);

            diag[i] = simd_min(simd_min(above, left_above), left);
		}

		if (k <= _tag_max_length) {
            diag[0] = simd_set1(k);
            diag[k] = simd_set1(k);
		} else if (k == 2 * _tag_max_length) {
            simd_store(reinterpret_cast<simd_t*>(&result), diag[startRow]);
            return result;
        }

#ifdef LEVEN_EARLY_TERMINATION
        const simd_t et_cmp = simd_cmpgt(diag[k/2], et_v);
        const int et_res = simd_movemask(et_cmp);

        // TODO adopt value to SSE/AVX

        if (et_res == _et_res_cmp) {
            simd_store(reinterpret_cast<simd_t*>(&result), diag[k/2]);
/*
            std::cout << "ET SIMD LD = ";
            for(int x=0 ; x != SIMD_WIDTH; ++x)
                std::cout << ' ' << setfill('0') << std::setw(2) << std::hex << ((int)result[x]);
            std::cout << dec << '\n';
*/
            return result;

        }
#endif // LEVEN_EARLY_TERMINATION

		// swap diagonals
		std::swap(diag, diag2);
    }

#undef simd_set1
#undef simd_setzero
#undef simd_add
#undef simd_min
#undef simd_cmpeq
#undef simd_cmpgt
#undef simd_movemask
#undef simd_store
}


// OpenMP SIMD ----------------------------------------------------------

template<typename T = unsigned char>
inline constexpr auto min_omp_simd(const T a, const T b) {
    return a < b ? a : b;
}

template<typename T = unsigned char>
inline constexpr auto cmp_omp_simd(const T a, const T b) {
    return a == b ? 0 : 1;
}

// 47 cycles on Westmere using sse, 16x16
// 79 cycles on Broadwell using avx, 32x16
template<typename T, size_t t_dim, int t_simd_width>
inline void transpose_decode_omp_simd( \
    const tag_encoded_vector_t& tag_encoded_vector,
    const indices_array_t& indices,
    T out[t_dim][t_simd_width]) {

    // collect data
    alignas(t_simd_width) tag_encoded_t in[t_simd_width];
    for(auto i=0; i!=t_simd_width; ++i) {
        in[i] = tag_encoded_vector[indices[i]];
    }

    const T* c_array = reinterpret_cast<const T*>(&in[0]);

    #pragma omp simd aligned(c_array, out : t_simd_width)
    for(int s = 0; s < t_simd_width; ++s) {
        for(size_t j = 0; j < t_dim/2; ++j) {
            out[2*j  ][s] = c_array[8 * s + j] & 0xF;
            out[2*j+1][s] = c_array[8 * s + j] >> 4;
        }
    }
}

// Westmere 1074.77 cycles (SSE 16)
// C2D 1214.77 cycles (SSE 16)
// Broadwell 987.472 (SSE 16)
// Broadwell 997.047 (AVX2 32)
template<typename T, size_t t_dim, int t_simd_width>
inline ld_array_t leven_omp_simd( \
    const tag_encoded_vector_t& tag_encoded_vector, \
    const indices_array_t& row_indices, \
    const indices_array_t& column_indices) {

    alignas(t_simd_width) T M[t_dim+1][t_dim+1][t_simd_width];
    alignas(t_simd_width) T a[t_dim][t_simd_width];
    alignas(t_simd_width) T b[t_dim][t_simd_width];

    transpose_decode_omp_simd<T, t_dim, t_simd_width>(tag_encoded_vector, row_indices, a);
    transpose_decode_omp_simd<T, t_dim, t_simd_width>(tag_encoded_vector, column_indices, b);

    for(size_t index = 0; index < t_dim+1; ++index) {

        #pragma omp simd aligned(M : t_simd_width)
        for(int s = 0; s < t_simd_width; ++s) {
            M[0][index][s] = index;
            M[index][0][s] = index;
        }
    }

    for(size_t row = 1; row < t_dim+1; ++row) {

        #pragma omp simd aligned(a, b, M : t_simd_width)
        for(size_t s = 0; s < t_simd_width; ++s) {

            for(size_t col = 1; col < t_dim+1; ++col) {
                const T delta = cmp_omp_simd<T>(a[row - 1][s], b[col - 1][s]);
                M[row][col][s] = \
                    min_omp_simd<T>( \
                        min_omp_simd<T>(M[row-1][col][s] + 1, M[row][col-1][s] + 1), \
                    M[row-1][col-1][s] + delta);
            }
        }
    }

    ld_array_t result;

    #pragma omp simd aligned(M : t_simd_width)
    for(size_t s = 0; s < t_simd_width; ++s)
        result[s] = M[t_dim][t_dim][s];

    return result;
}




#endif // USE_SSE


#endif // _ECLEVEN_H_
