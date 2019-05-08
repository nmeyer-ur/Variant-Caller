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

#ifndef _ECCLUSTER_H_
#define _ECCLUSTER_H_

#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <array>
#include <algorithm>
#include <map>
#include "ecdefs.hpp"


class Levenshtein {

public:
    Levenshtein(const Parameters& p, const tag_encoded_vector_t& u, distance_matrix_t& d):
        params(p), tag_encoded_vector(u), M(d) {

        dim = (int)tag_encoded_vector.size();
        init_tag_base_freq_vector();
        eval_M();
    }

    void update_statistics(Statistics& ext_statistics) {

        ext_statistics.leven_complete_evaluations += statistics.leven_complete_evaluations;
        ext_statistics.leven_total_evaluations += statistics.leven_total_evaluations;
    }

private:
    Statistics statistics;
    const Parameters& params;
    const tag_encoded_vector_t& tag_encoded_vector;
    counter_vector_t tag_base_freq_vector;
    distance_matrix_t& M;
    int dim;

    void init_tag_base_freq_vector(void) {

        tag_base_freq_vector.reserve(dim);
        for(const auto tag : tag_encoded_vector) {
            tag_base_freq_vector.push_back(tag_base_freq(tag));
        }
    }


// put this block here because it is not inlined correctly

    /*  TAG bases frequencies encoding

        Microbenchmarks (g++ 5.4/6.x -O3 -march=native)

        Intel Core2 Duo P8700 @ 2.53GHz             36 cycles
        Intel Xeon E5645 Westmere @ 2.40GHz         30 cycles
        Intel Xeon E5-2620 v4 Broadwell @ 2.10GHz   17 cycles

        faster than any SSE intrinsics implementation
    */
    counter_t tag_base_freq(tag_encoded_t u) {

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

        faster than any SSE intrinsics implementation
    */
    leven_dist_t leven_lower_bound_by_base_freq(counter_t u, counter_t v) {

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
    leven_dist_t tag_hamming(tag_encoded_t u, tag_encoded_t v) {

        counter_t result = 0;
        counter_t mask = 0xF;
    	const counter_t x = u ^ v;

        for(int i=0; i != _tag_max_length; ++i) {
    		result += ((x & mask) != 0);
    		mask = mask << 4;
        }
        return (leven_dist_t)result;
    }

// ------------------------------------- inlining issue --------------------

    void eval_M(void) {

        // SIMD support
        indices_array_t row_array;
        indices_array_t column_array;
        ld_array_t ld_array;
        int array_index = 0;

        statistics.leven_total_evaluations += (dim*dim - dim) / 2;
        M(0,0) = 0;

        for(auto i=1; i != dim; ++i) {

            M(i,i) = 0;

            const auto tag1_base_freq = tag_base_freq_vector[i];
            for(auto j=0; j != i; ++j) {

                // Levenshtein distance lower bound check
                leven_dist_t ld = 0;

#ifdef LEVEN_LOWER_BOUND
                const auto tag2_base_freq = tag_base_freq_vector[j];
                ld = leven_lower_bound_by_base_freq(tag1_base_freq, tag2_base_freq);

                if (ld > params.ld) {
                    M(i,j) = M(j,i) = ld;
                    continue;
                }
#endif

                // Levenshtein distance evaluation
                row_array[array_index]    = i;
                column_array[array_index] = j;
                ++array_index;

                if (array_index == SIMD_WIDTH) {

#if SIMD_WIDTH > 1
/*
                    ld_array = leven_sse(tag_encoded_vector, \
                        row_array, column_array, (int)params.ld);
*/
                    ld_array = leven_omp_simd<unsigned char, _tag_max_length, SIMD_WIDTH> \
                        (tag_encoded_vector, row_array, column_array);

// test
/*
                    ld_array_t ld_array0;
                    ld_array0 = leven_sse(tag_encoded_vector, \
                        row_array, column_array, 100);

                    ld_array_t ld_array2;
                    for(auto k=0; k != SIMD_WIDTH; ++k) {
                        ld_array2[k] = (char)leven_gen(
                            tag_encoded_vector[row_array[k]],
                            tag_encoded_vector[column_array[k]],
                            (int)params.ld);
                    }

                    ld_array_t ld_array3;
                    for(auto k=0; k != SIMD_WIDTH; ++k) {
                        ld_array3[k] = (char)leven_gen(
                            tag_encoded_vector[row_array[k]],
                            tag_encoded_vector[column_array[k]],
                            100);
                    }

            using namespace std;
            cout << "SSE LD = ";
            for(int x=0 ; x != SIMD_WIDTH; ++x)
                cout << ' ' << setfill('0') << setw(2) << hex << ((int)ld_array0[x]);
            cout << dec << '\n';

            for(int x=0 ; x != SIMD_WIDTH; ++x)
                if (ld_array[x] > params.ld) ld_array[x] = params.ld + 1;

            cout << "SSE LD = ";
            for(int x=0 ; x != SIMD_WIDTH; ++x)
                cout << ' ' << setfill('0') << setw(2) << hex << ((int)ld_array[x]);
            cout << dec << '\n';

            cout << "GEN LD = ";
            for(int x=0 ; x != SIMD_WIDTH; ++x)
                cout << ' ' << setfill('0') << setw(2) << hex << ((int)ld_array3[x]);
            cout << dec << '\n';

            cout << "GEN LD = ";
            for(int x=0 ; x != SIMD_WIDTH; ++x)
                cout << ' ' << setfill('0') << setw(2) << hex << ((int)ld_array2[x]);
            cout << dec << '\n';

                    assert(ld_array == ld_array2);
                    assert(ld_array0 == ld_array3);
*/

// ----


#else // generic version
                    for(auto k=0; k != array_index; ++k) {
                        ld_array[k] = (char)leven_gen(
                            tag_encoded_vector[row_array[k]],
                            tag_encoded_vector[column_array[k]],
                            (int)params.ld);
                    }
#endif

                    for(auto k=0; k != SIMD_WIDTH; ++k) {
                        const auto row = row_array[k];
                        const auto col = column_array[k];
                        M(row, col) = M(col, row) = (leven_dist_t)ld_array[k];
                    }

                    array_index = 0;
                    statistics.leven_complete_evaluations += SIMD_WIDTH;
                }
            }
        }

        // remainders of Levenshtein distance matrix
        if (array_index > 0) {
#if SIMD_WIDTH > 1
            // complete array
            for(auto k=array_index; k < SIMD_WIDTH; ++k) {
                row_array[k]    = row_array[array_index - 1];
                column_array[k] = column_array[array_index - 1];
            }
            ld_array = leven_sse(tag_encoded_vector, \
                row_array, column_array, (int)params.ld);
#else // generic
            for(auto k=0; k != array_index; ++k) {
                ld_array[k] = (char)leven_gen(
                    tag_encoded_vector[row_array[k]],
                    tag_encoded_vector[column_array[k]],
                    (int)params.ld);
            }
#endif
            for(auto k=0; k != array_index; ++k) {
                const auto row = row_array[k];
                const auto col = column_array[k];
                M(row, col) = M(col, row) = (leven_dist_t)ld_array[k];
            }
            statistics.leven_complete_evaluations += array_index;
        }
    }
};


class Cluster {

public:
    cluster_map_t cluster_map_final;

    Cluster(const Parameters& p, const read_vector_t& rv):
        params(p), read_vector(rv), clustering_completed(false) {

        build_initial_cluster_map();
    }

    void cluster(void) {

        if (params.cluster_mode == "none")
            cluster_map_tmp = cluster_map_initial;
        else if (params.cluster_mode == "dbscan")
//            fun();
            fun2();
        else
            assert(!"Invalid cluster mode. Exiting.");

        post_process_cluster();
        clustering_completed = true;
    }

    bool is_valid() { return (cluster_map_final.size() > 0); }

    void display_stats(void) {

        if (clustering_completed == false) {
            _display_stats(cluster_map_initial);
            if (params.verbosity >= 4)
                _display_cluster(cluster_map_initial);
        } else {
            _display_stats(cluster_map_final);
            if (params.verbosity >= 4)
                _display_cluster(cluster_map_final);
        }
    }

    const cluster_map_t& cluster_map(void) { return cluster_map_final; }

    void update_statistics(Statistics& ext_statistics) {

        ext_statistics.tag_before_clustering += cluster_map_initial.size();
        ext_statistics.tag_before_clustering_vector.push_back(cluster_map_initial.size());
        ext_statistics.tag_after_clustering += cluster_map_final.size();
        if (cluster_map_final.size() > 0)
            ext_statistics.tag_after_clustering_vector.push_back(cluster_map_final.size());

        // update reads statistics
        for(const auto& mapping : cluster_map_initial)
            ext_statistics.reads_before_clustering += (counter_t)mapping.second.size();
        for(const auto& mapping : cluster_map_final)
            ext_statistics.reads_after_clustering += (counter_t)mapping.second.size();

        // update TAG statistics
        for(const auto& mapping : cluster_map_final) {
            const auto tag = mapping.second[0]->tag_encoded;
            ext_statistics.tag_counter_map[tag] += 1;
        }

        ext_statistics.leven_complete_evaluations += statistics.leven_complete_evaluations;
        ext_statistics.leven_total_evaluations += statistics.leven_total_evaluations;
        ext_statistics.temp_variable += statistics.temp_variable;
    }

private:
    Statistics statistics;
    const Parameters& params;
    const read_vector_t& read_vector;
    bool clustering_completed;
    cluster_map_t cluster_map_initial;
    cluster_map_t cluster_map_tmp;

    void _display_cluster(const cluster_map_t& cluster_map) {

        for(const auto& mapping : cluster_map) {
            const auto tag = tag_encoded_to_string(mapping.second[0]->tag_encoded);
            std::cout << tag << " <-- ";
            for(const auto& read : mapping.second) {
                std::cout << tag_encoded_to_string(read->tag_encoded) << ' ';
            }
            std::cout << '\n';
        }
    }

    void _display_stats(const cluster_map_t& cluster_map) {

        const size_t cluster_size = cluster_map.size();
        size_t min_size = 0xfffffff;
        size_t max_size = 0;

        std::map<int, int> freq;
        for(const auto& mapping : cluster_map) {
            const size_t s = mapping.second.size();
            min_size = std::min(s, min_size);
            max_size = std::max(s, max_size);
            freq[s] += 1;
        }

        // general statistics
        if (clustering_completed == false)
            std::cout << "TAGs before clustering: " << cluster_size << '\n';
        else
            std::cout << "TAGs after clustering and read count filtering: " << cluster_size << '\n';

        if (cluster_size == 0)
            return;

        std::cout << "Unique TAGs " << cluster_size << ", read count ";
        if (min_size != max_size)
            std::cout << min_size << "-" << max_size << '\n';
        else
            std::cout << max_size << '\n';


        // frequency count, highest to lowest
        std::cout << "TAG count\tread count" << '\n';
        for(auto it = freq.crbegin(); it != freq.crend(); ++it)
            std::cout << it->second << "\t\t" << it->first << '\n';
    }

    // build mapping between TAGs and reads
    // TAG_1 : [read_1, read_2, ...]
    // TAG_2 : [read_5, read_8, ...]
    void build_initial_cluster_map(void) {

        for(const auto& read : read_vector) {
            const auto tag = read->tag_encoded;
            cluster_map_initial[tag].push_back(read);
        }
    }

    // post processing: create the final mapping
    // discard read vectors with insufficient read count
    void post_process_cluster() {

        for(const auto& mapping : cluster_map_tmp) {
            const auto& read_vector = mapping.second;

            const bool check_cov = check_coverage(read_vector, params.min_reads);
            if (check_cov == true) {
                const auto tag = read_vector[0]->tag_encoded;
                cluster_map_final[tag] = read_vector;
            }
        }
    }

    void fun() {

        auto size = cluster_map_initial.size();
        if (size < 2) {
            cluster_map_tmp = cluster_map_initial;
            return;
        }

        // get the keys of the map
        std::vector<tag_encoded_t> tag_encoded_vector;
        tag_encoded_vector.reserve(size);

        std::vector<tag_base_freq_array_t> tag_base_freq_array;
        tag_base_freq_array.reserve(size);

        for(const auto& mapping : cluster_map_initial) {
            auto& tag_encoded = mapping.first;
            tag_encoded_vector.push_back(tag_encoded);
            tag_base_freq_array.push_back(tag_encoded_to_base_freq_array(tag_encoded));
        }

        // evaluate Levensthein distance
        distance_matrix_t d;
        d.resize(size, size);

        d(0,0) = 0;
        for(decltype(size) i=1; i<size; ++i) {
            d(i,i) = 0;
            auto& tag1 = tag_encoded_vector[i];
            auto& tag1_base_freq = tag_base_freq_array[i];

            for(decltype(size) j=0; j<i; ++j) {
                auto& tag2 = tag_encoded_vector[j];
                auto& tag2_base_freq = tag_base_freq_array[j];

                leven_dist_t ld = 0;
#ifdef LEVEN_LOWER_BOUND
                ld = (leven_dist_t)leven_lower_bound_array(tag1_base_freq, tag2_base_freq);
#endif
                if (ld <= params.ld) {
#ifndef LEVEN_EARLY_TERMINATION
// diagonal, full
                    ld = (leven_dist_t)leven_1a<leven_dist_t>(tag1, tag2);
#else
// row, early termination
                    ld = (leven_dist_t)leven_2a<leven_dist_t>(tag1, tag2);
#endif
                    statistics.leven_complete_evaluations += 1;
                }

                d(i,j) = d(j,i) = ld;
                statistics.leven_total_evaluations += 1;
/*
                std::cout << tag_encoded_to_string(tag1) << ' '
                    << tag_encoded_to_string(tag2) << ' '
                    << ld << ' ' << lower << '\n';
*/
            }
        }

        DBSCAN<eigen_vector_t, distance_matrix_t, leven_dist_t> cl_dbscan(params.ld);
        cl_dbscan.fit_precomputed(d);

        std::vector<int> labels = cl_dbscan.get_labels();

        // extend
        for(decltype(size) i=0; i<size; ++i) {
            auto tag = tag_encoded_vector[i];
            auto bin = labels[i];
            if (bin != -1)  {
                cluster_map_tmp[bin].insert( \
                    cluster_map_tmp[bin].end(), \
                    cluster_map_initial[tag].begin(), \
                    cluster_map_initial[tag].end());
            }
        }
    }

    void fun2() {

        int size = (int)cluster_map_initial.size();
        if (size < 2) {
            cluster_map_tmp = cluster_map_initial;
            return;
        }

        // get the TAGs
        tag_encoded_vector_t tag_encoded_vector;
        tag_encoded_vector.reserve(size);

        for(const auto& mapping : cluster_map_initial) {
            const auto tag_encoded = mapping.first;
            tag_encoded_vector.push_back(tag_encoded);
        }

        // evaluate Levensthein distance
        distance_matrix_t d;
        d.resize(size, size);

        Levenshtein lev(params, tag_encoded_vector, d);
        lev.update_statistics(statistics);

        DBSCAN<eigen_vector_t, distance_matrix_t, leven_dist_t> cl_dbscan(params.ld);
        cl_dbscan.fit_precomputed(d);

        std::vector<int> labels = cl_dbscan.get_labels();

        // extend
        for(decltype(size) i=0; i != size; ++i) {
            const auto tag = tag_encoded_vector[i];
            const auto bin = labels[i];
            if (bin != -1)  {
                cluster_map_tmp[bin].insert( \
                    cluster_map_tmp[bin].end(), \
                    cluster_map_initial[tag].begin(), \
                    cluster_map_initial[tag].end());
            }
        }
    }


    // adopted to 16 characters TAG from
    // https://en.wikipedia.org/wiki/Levenshtein_distance
    // tested with char and int
    // int is slightly faster than char
    template<typename T>
    T leven_1a(tag_encoded_t a, tag_encoded_t b) {

        std::array<T, _tag_max_length+1> v0, v1;
        std::array<T, _tag_max_length> s, t;

        for(T i=0; i<_tag_max_length; i++) {

            v0[i] = i;
            s[i] = (T)(a & 0xf);
            a = a >> 4;
            t[i] = (T)(b & 0xf);
            b = b >> 4;
        }
        v0[_tag_max_length] = _tag_max_length;

        for(T i=0; i <_tag_max_length; i++) {

            v1[0] = i + 1;

            for(T j=0; j<_tag_max_length; j++) {
                T cost = (s[i] == t[j]) ? 0 : 1;
                v1[j + 1] = std::min(std::min(v1[j] + 1, v0[j + 1] + 1), v0[j] + cost);
            }

            // copy row
            for(T j=0; j<=_tag_max_length; j++)
                v0[j] = v1[j];
        }
        return v1[_tag_max_length];
    }

    /**
     * row-based variant, slower than diagonal-based variant
     */
    template<typename T>
    T leven_2a(tag_encoded_t a, tag_encoded_t b) {

        std::array<T, _tag_max_length> v, s, t; // _attribute_ ((aligned (16)));
        T ret = 0;

        T length_a = 0, length_b = 0;

        for(T i=0; i<_tag_max_length; ++i) {
            v[i] = i+1;
            s[i] = (T)(a & 0xf);
            if (s[i] != _bin_N_base) ++length_a;
            a = a >> 4;
            t[i] = (T)(b & 0xf);
            b = b >> 4;
            if (t[i] != _bin_N_base) ++length_b;
        }

        T length = (T)std::max((int)length_a, (int)length_b);

//        for(int i=0; i<_tag_max_length; ++i) {
        for(T i=0; i<length; ++i) {

            T tmp = i;
            ret = i + 1;

            for (T j=0; j<length; ++j) {
//            for (int j=0; j<_tag_max_length; ++j) {

                T tmp2 = (t[i] == s[j]) ? tmp : tmp + 1;
                tmp = v[j];
                // amazingly ugly and confusing
                ret = v[j] = tmp > ret ? tmp2 > ret ? ret + 1 : tmp2 : tmp2 > tmp ? tmp + 1 : tmp2;
            }

            // early termination
            if (v[i] > params.ld) return v[i];
        }
        return ret;
    }

    static inline int leven_lower_bound_array(const tag_base_freq_array_t& freq_a, const tag_base_freq_array_t& freq_b) {

        // sum over absolute values of the frequency difference
        int sum = 0;
        for(auto i=0; i<_num_bases; ++i)
            sum += std::abs(freq_a[i] - freq_b[i]);

        // correct for double counting
        return (sum >> 1);
    }

    static inline tag_base_freq_array_t tag_encoded_to_base_freq_array(tag_encoded_t a) {

        tag_base_freq_array_t freq = {0};

        for(auto i=0; i<_tag_max_length; ++i) {
            char base_a = (char)(a & 0xf);
            freq[base_a] += 1;
            a = a >> 4;
        }

        return freq;
    }

};

#endif // _ECCLUSTER_H_
