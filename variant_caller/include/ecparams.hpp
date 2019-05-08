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

#ifndef _ECPARAMS_H_
#define _ECPARAMS_H_

#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <bitset>
#include <assert.h>
#include "ecdefs.hpp"


class Parameters {

public:
    std::string alignment_file;
    std::string ref_fasta_file;
    std::string cluster_mode;
    mapq_t min_mapq;
    int min_reads;
    leven_dist_t ld;
    int verbosity;
    bool stats;
    bool tmpvc;
    bool forward;
    bool reverse;
    tag_mask_t tag_mask;
    orientation_t orientation;
    flag_t flag_mask;
    flag_t pair_mask;
    float base_acceptance_threshold;
    int phred_score_threshold;
    int hts_threads;
    int ec_threads;
    bool region_check;
    int simd_w;
    std::vector<std::string> cmdline;
    std::vector<std::string> region_vector;

    Parameters(int argc, char* argv[]) {

        // command line string
        cmdline.assign(argv, argv + argc);

        // initialize default values
        alignment_file              = "none";
        ref_fasta_file              = "none";
        tag_mask                    = 0;
        cluster_mode                = "none";
        region_check                = false;
        base_acceptance_threshold   = 0.7;
        hts_threads                 = 0;
        ec_threads                  = omp_get_max_threads();
        phred_score_threshold       = 0;
        simd_w                      = SIMD_WIDTH;
        min_mapq                    = 0;
        min_reads                   = 5;
        ld                          = 2;
        stats                       = false;
        tmpvc                       = false;
        forward                     = true;
        reverse                     = true;
        orientation                 = 0;
        verbosity                   = 0;

        // Filter reads by flag
        flag_mask = _samflag_mask_read_unmapped | \
                    _samflag_mask_mate_unmapped | \
                    _samflag_mask_not_primary_alignment | \
                    _samflag_mask_read_fails_quality_check | \
                    _samflag_mask_read_is_pcr_or_optical_duplicate | \
                    _samflag_mask_supplementary_alignment;

        pair_mask = _samflag_mask_read_paired | \
                    _samflag_mask_read_mapped_in_proper_pair;

        const struct option longopts[] = {
            {"version",     no_argument,        0, 'V'},
            {"help",        no_argument,        0, 'h'},
            {"stats",       no_argument,        0, 's'},
            {"tmpvc",       no_argument,        0, 'T'},
            {"forward",     no_argument,        0, 'F'},
            {"reverse",     no_argument,        0, 'R'},
            {"mapq",        required_argument,  0, 'm'},
            {"phredscore",  required_argument,  0, 'p'},
            {"coverage",    required_argument,  0, 'c'},
            {"distance",    required_argument,  0, 'd'},
            {"input",       required_argument,  0, 'i'},
            {"reference",   required_argument,  0, 'f'},
            {"clustermode", required_argument,  0, 'C'},
            {"tag",         required_argument,  0, 'u'},
            {"region",      required_argument,  0, 'r'},
            {"orientation", required_argument,  0, 'O'},
            {"threads",     required_argument,  0, 't'},
            {"htsthreads",  required_argument,  0, 'H'},
            {0,0,0,0},
        };

        //opterr = 1;
        int opt = 0;
        int index = 0;

        while(opt != -1) {

            std::string tmp;

            opt = getopt_long(argc, argv, \
                "hvVTFRsi:f:r:C:m:c:d:u:r:t:H:p:b:O:", longopts, &index);
            if (opt == -1) break;

            switch(opt) {
                case 'V':
                    display_version();
                    exit(EXIT_SUCCESS);
                    break;
                case 'h':
                    display_help();
                    exit(EXIT_SUCCESS);
                    break;
                case 's': stats = true; break;
                case 'T': tmpvc = true; break;
                case 'm': min_mapq = (mapq_t)std::stoi(optarg); break;
                case 'c': min_reads = std::stoi(optarg); break;
                case 'd': ld = (leven_dist_t)std::stoi(optarg); break;
                case 'i': alignment_file = std::string(optarg); break;
                case 'f': ref_fasta_file = std::string(optarg); break;
                case 'C': cluster_mode = std::string(optarg); break;
                case 'v': ++verbosity; break;
                case 'F': reverse = false; break;
                case 'R': forward = false; break;
                case 't': ec_threads = std::stoi(optarg); break;
                case 'H': hts_threads = std::stoi(optarg); break;
                case 'p': phred_score_threshold = std::stoi(optarg); break;
                case 'b': base_acceptance_threshold = std::stof(optarg); break;
                case 'r': region_vector.push_back(std::string(optarg)); break;
                case 'u':
                    tmp = std::string(optarg);
                    if (tmp == "F1") tag_mask.set(_tag_f1_index);
                    else if (tmp == "R1") tag_mask.set(_tag_r1_index);
                    else if (tmp == "F2") tag_mask.set(_tag_f2_index);
                    else if (tmp == "R2") tag_mask.set(_tag_r2_index);
                    else {
                        std::cout << "Caught invalid TAG. Exiting.\n";
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'O':
                    tmp = std::string(optarg);
                    if (tmp == "all") {
                        orientation.set(_orientation_index_FR);
                        orientation.set(_orientation_index_RF);
                        orientation.set(_orientation_index_FF);
                        orientation.set(_orientation_index_RR);
                    }
                    else if (tmp == "FR")
                        orientation.set(_orientation_index_FR);
                    else if (tmp == "RF")
                        orientation.set(_orientation_index_RF);
                    else if (tmp == "FF")
                        orientation.set(_orientation_index_FF);
                    else if (tmp == "RR")
                        orientation.set(_orientation_index_RR);
                    else {
                        std::cout << "Caught invalid orientation. Exiting.\n";
                        exit(EXIT_FAILURE);
                    }
                    break;
                default:
                    std::cout << "Caught invalid parameter. Exiting.\n";
                    exit(EXIT_FAILURE);
                    break;
            }
        }

        // File checks
        if (access(alignment_file.c_str(), F_OK) == -1) {
            std::cout << "Input BAM file not found. Exiting.\n";
            exit(EXIT_FAILURE);
        }

        if (access(ref_fasta_file.c_str(), F_OK) == -1) {
            std::cout << "Reference file not found. Exiting.\n";
            exit(EXIT_FAILURE);
        }

        // Boundary checks
        //if ((min_mapq < 0) || (min_mapq > 60)) {
        if (min_mapq > 60) {
            std::cout << "MAPQ score invalid. Exiting.\n";
            exit(EXIT_FAILURE);
        }

        if (phred_score_threshold < 0) {
            std::cout << "Phred score threshold invalid. Exiting.\n";
            exit(EXIT_FAILURE);
        }

        if ((base_acceptance_threshold <= 0.5) || (base_acceptance_threshold > 1.0)) {
            std::cout << "Base acceptance threshold invalid. Exiting.\n";
            exit(EXIT_FAILURE);
        }

        if (min_reads < 1) {
            std::cout << "Minimum reads invalid. Exiting.\n";
            exit(EXIT_FAILURE);
        }

        if (ld < 0) {
            std::cout << "Levensthein distance invalid. Exiting.\n";
            exit(EXIT_FAILURE);
        }

        // forward / reverse reads
        if (!(forward || reverse)) {
            std::cout << "No reads. Exiting.\n";
            exit(EXIT_FAILURE);
        }

        // TAG
        if (tag_mask.none())
            tag_mask.set(_tag_f1_index);

        // orientation
        if (orientation.none()) {
            orientation.set(_orientation_index_FR);
            orientation.set(_orientation_index_RF);
            orientation.set(_orientation_index_FF);
            orientation.set(_orientation_index_RR);
        }

		if (orientation.test(_orientation_index_FF) && (orientation.count() == 1) && !forward) {
            std::cout << "Illegal orientation and read direction. Exiting.\n";
            exit(EXIT_FAILURE);

		}

		if (orientation.test(_orientation_index_RR) && (orientation.count() == 1) && !reverse) {
            std::cout << "Illegal orientation and read direction. Exiting.\n";
            exit(EXIT_FAILURE);
		}

        // region selection
        auto& rv = region_vector;
        if (rv.size() == 0) {
            rv.push_back("all");
        } else {
            // check for "all"
            if (std::find(rv.begin(), rv.end(), "all") == rv.end()) {
                std::sort(rv.begin(), rv.end());
                rv.erase(std::unique(rv.begin(), rv.end()), rv.end());
                region_check = true;
            } else {
                rv.clear();
                rv.push_back("all");
            }
        }

        if (verbosity > 2) {
            ec_threads = 1;
            std::cout << "Running single-threaded.\n";
        }

        if ((ec_threads < 1) || (hts_threads < 0)) {
            std::cout << "Number of threads invalid. Exiting.\n";
            exit(EXIT_FAILURE);
        }
    }

    void display_params(void) {

        std::string orientation_str;
        if (orientation.test(_orientation_index_FR))
            orientation_str += "FR ";
        if (orientation.test(_orientation_index_RF))
            orientation_str += "RF ";
        if (orientation.test(_orientation_index_FF))
            orientation_str += "FF ";
        if (orientation.test(_orientation_index_RR))
            orientation_str += "RR ";

        std::string tag_str;
        if (tag_mask.test(_tag_f1_index))
            tag_str += "F1 ";
        if (tag_mask.test(_tag_f2_index))
            tag_str += "F2 ";
        if (tag_mask.test(_tag_r1_index))
            tag_str += "R1 ";
        if (tag_mask.test(_tag_r2_index))
            tag_str += "R2 ";

        std::cout << "Runtime parameters\n\n";
        std::cout << "Command                    ";
        for(const auto& c : cmdline)
            std::cout << c << " ";
        std::cout << '\n';
        std::cout << "SIMD width                 " << simd_w;
        switch(simd_w) {
            case  1: std::cout << " (generic scalar)\n"; break;
            case 16: std::cout << " (SSE)\n"; break;
            case 32: std::cout << " (AVX)\n"; break;
            case 64: std::cout << " (AVX512)\n"; break;
            default: std::cout << " (unknown)\n"; break;
        }
        std::cout << "Error correction threads   " << ec_threads << '\n';
        std::cout << "HTSlib I/O threads         " << hts_threads << '\n';
        std::cout << "BAM file                   " << alignment_file << '\n';
        std::cout << "Reference fasta file       " << ref_fasta_file << '\n';
        std::cout << "Region                     ";
        for(const auto& r : region_vector)
            std::cout << r << ' ';
        std::cout << '\n';
        std::cout << "Cluster mode               " << cluster_mode << '\n';
        std::cout << "Max. Levensthein distance  " << (int)ld << '\n';
        std::cout << "Min. MAPQ score            " << (int)min_mapq << '\n';
        std::cout << "Phred score threshold      " << (int)phred_score_threshold << '\n';
        std::cout << "Error correction threshold " << base_acceptance_threshold << '\n';
        std::cout << "Min. reads per allele      " << min_reads << '\n';
        std::cout << "Read orientations          " << orientation_str << '\n';
        std::cout << "Read selection             " \
            << (forward ? "forward " : "") << (reverse ? "reverse" : "") << '\n';
        std::cout << "TAG                        " << tag_str << '\n';
        std::cout << "TAG mask                   " << tag_mask << '\n';
        std::cout << "Statistics output          " << (stats ? "yes" : "no") << '\n';
        std::cout << "display intermediate VC    " << (tmpvc ? "yes" : "no") << '\n';
        std::cout << "Verbosity level            " << verbosity << '\n';
    }

    void display_version(void) {

        std::cout << "ec " << _version << '\n';
        std::cout << "Allele-specific variant caller\n";
        std::cout << "Nils Meyer, University of Regensburg, Department of Physics\n";
    }

    void display_help(void) {

        const char* usage =
        "usage: ec <options>\n"
        "\n"
        "required arguments:\n"
        "  -i INPUT, --input INPUT\n"
        "                        input alignment file, BAM format\n"
        "  -f REFERENCE, --reference REFERENCE\n"
        "                        reference FASTA file\n"
        "\n"
        "optional arguments:\n"
        "  -h, --help            show this help message\n"
        "  -s, --stats           display statistics (default: false)\n"
        "  -t THREADS, --threads THREADS\n"
        "                        number of threads used for error correction (default: max)\n"
        "  -H HTSTHREADS, --htsthreads HTSTHREADS\n"
        "                        thread number passed to htslib (default: 0)\n"
        "                        *for performance testing only*\n"
        "  -m MAPQ, --mapq MAPQ  minimum MAPQ score (0-60) (default: 0)\n"
        "  -p PHREDSCORE, --phredscore PHREDSCORE\n"
        "                        convert bases below phred score threshold to 'N' (default: 0)\n"
        "  -b THRESHOLD, --base THRESHOLD\n"
        "                        error correction threshold (>0.5-1.0) (default: 0.7)\n"
        "  -c COVERAGE, --coverage COVERAGE\n"
        "                        minimum coverage per block of overlapping reads (>0) (default: 5)\n"
        "  -d DISTANCE, --distance DISTANCE\n"
        "                        maximum Levenshtein distance for TAG clustering (default: 2)\n"
        "  -C {none, dbscan}, --clustermode {none, dbscan} (default: none)\n"
        "  -r REGION, --region REGION\n"
        "                        region/contig without position information, additive (default: all)\n"
        "  -u {F1, F2, R1, R2}, --tag {F1, F2, R1, R2}\n"
        "                        TAG, additive (default: F1)\n"
        "  -O {FR, RF, FF, RR, all}, --orientation {FR, RF, FF, RR, all}\n"
        "                        read orientation, additive (default: all)\n"
        "  -F, --forward         evaluate forward reads only (default: false)\n"
        "  -R, --reverse         evaluate reverse reads only (default: false)\n"
        "  -v, --verbosity       increase output verbosity (default: none)\n"
        "  -T, --tmpvc           display variant calling in intermediate representation (default: false)\n"
        "  --version             version information"
        ;

        std::cout << usage << '\n';
    }
};

#endif // _ECPARAMS_H_
