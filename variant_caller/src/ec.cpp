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

#include <iostream>
#include <omp.h>
#include "Eigen/Core"
#include "ecdefs.hpp"
#include "ecparams.hpp"
#include "ecstats.hpp"
#include "ecsam.hpp"
#include "ecreads.hpp"
#include "ecblock.hpp"
#include "ecleven.hpp"
#include "ecdbscan.hpp"
#include "eccluster.hpp"
#include "ecallele.hpp"
#include "ecvcall.hpp"
#include "ecref.hpp"

/*
ec

Allele-specific variant caller
Nils Meyer, University of Regensburg, Department of Physics


Version history
===============

v0.12   crc32 support

v0.11   removed namespace std
        fixed block output
        introduced cli parameter for tmpvc output
        corrected tmpvc output order
        corrected naming convention mutation -> SNP
        introduced read pair orientation
        corrected TAG ordering
        region/contig selection via cli

v0.10   revised threading scheme: one thread per contig, no serial processing

TODO
====

    * vcf file output
    * correct insertions statistics in tmpvcf

IMPROVEMENTS
============

    * cluster_map_tmp implementation yields correct results,
      but the implementation mixes TAGs and labels
*/


void process_block(const Parameters& params,
    const std::string region,
    Block& block,
    Reference_file& reference_file,
    std::string& output_str,
    Statistics& statistics) {

    const auto verbosity = params.verbosity;

    const bool block_valid = block.is_valid();
    if (block_valid == false)
        return;

    if (verbosity >= 2)
        block.display_stats();

    block.update_statistics(statistics);

    const position_t start = block.start;
    const position_t end = block.end;

    read_vector_t& rv = block.read_vector;
    Cluster cluster(params, rv);

    if (verbosity >= 3)
        cluster.display_stats();

    cluster.cluster();
    cluster.update_statistics(statistics);

    if (verbosity >= 3)
        cluster.display_stats();

    if (!cluster.is_valid())
        return;

    const cluster_map_t& cluster_map = cluster.cluster_map();

    base_vector_t reference;
    reference_file.fetch(reference, region, start, end);

    Variant_caller var_call(params, reference, \
        cluster_map, region, start, end);

    var_call.variant_call();
    var_call.tmpvcf();

    if (verbosity >= 4) {
        var_call.display_tmpvcf();
    }

    var_call.gather_consensus_stats();
    var_call.update_statistics(statistics);
    output_str += var_call.master_str;
}



int main(int argc, char* argv[]) {

    Parameters params(argc, argv);
    params.display_params();

    const auto verbosity = params.verbosity;

    Alignment_file global_alignment_file(params);
    const metadata_t global_metadata = global_alignment_file.get_metadata();

    const int global_num_regions = global_metadata.size();

    if (global_num_regions == 0) {
        std::cout << "No region to process. Exiting.\n";
        exit(EXIT_FAILURE);
    }

    Statistics global_statistics;
    std::vector<std::string> output_string_vector(global_num_regions);

    if (verbosity >= 1) {
        constexpr int w = 12;
        std::cout << "\n\nNumber of selected regions: " << global_num_regions << '\n';
        std::cout << "\n  Region      read count      region index\n";
        for(auto i=0; i != global_num_regions; ++i) {
            std::cout << std::setw(w) << global_metadata[i].name << "\t" \
            << std::setw(w) << global_metadata[i].num_reads << "\t" \
            << std::setw(w) << global_metadata[i].index << '\n';
        }
        std::cout << "\n\nProgress\n\n";
    }

    // parallel processing
    const auto threads = params.ec_threads;
    auto global_region_index = 0; // region id

    omp_set_num_threads(threads);

    #pragma omp parallel \
        shared(global_statistics, output_string_vector, global_region_index) \
        firstprivate(params, global_metadata, global_num_regions)
    {
        Reference_file reference_file(params);

        std::vector<bam1_t*> aln_vector;
        aln_vector.reserve(_init_prefetch_vector_size);
        auto statistics = Statistics();
        std::string output_str;
        output_str.reserve(_init_output_string_size);

        Alignment_file alignment_file(params);
        bam_hdr_t* aln_header = alignment_file.aln_header();
        Read_preprocessing preprocessing = Read_preprocessing(params, aln_header);

        int region_index = -1;
        int output_region_index = -1;
        std::string region_name;
        bool all_regions_eof = false;

        // loop until all reads are processed
        while(true) {

            // select region
            #pragma omp critical (init_region)
            {
                if (global_region_index < global_num_regions) {
                    if (params.verbosity >= 1 && region_index >= 0) {
                        std::cout << region_name << " ... done\n";
                    }
                    region_index = global_region_index;
                    region_name = global_metadata[region_index].name;
                    output_region_index = global_metadata[region_index].index;
                    if (params.verbosity >= 1)
                        std::cout << region_name << " ... processing\n";
                    alignment_file.select_region(region_name);
                    ++global_region_index;
                } else
                    all_regions_eof = true;
            } // end omp critical

            // break if all regions done
            if (all_regions_eof == true)
                break;

            //output_str.clear();
//            aln_vector.clear();
//            block.clear();
            output_str.clear();

            Block block = Block(params, region_name);
            position_t section_start = 0xfffffff;
            position_t section_end = 0;
            bool eof = false;

            while(!eof) {

                // fetch reads
                {
                    while(!eof) {
                        bool new_section = false;

                        const bool aln_pending = alignment_file.fetch();
                        if (aln_pending) {

                            const position_t read_start = alignment_file.read_start();
                            const position_t read_end = alignment_file.read_end();

                            // check section boundaries
                            // new section
                            if (read_start < section_start) {
                                new_section = true;
                                // process reads if prefetch limit is reached
                                if (aln_vector.size() >= _num_prefetch_reads)
                                    break;
                            }

                            // if a new section starts reset section
                            if (new_section == true) {
                                section_start = read_start;
                                section_end = read_end;
                            }

                            // no new section, check boundaries
                            if (read_start >= (section_end + _section_gap)) {
                                // sufficient reads for processing
                                if (aln_vector.size() >= _num_prefetch_reads)
                                    break;
                            }

                            section_end = std::max(section_end, read_end);
                            auto& aln = alignment_file.aln();
                            alignment_file.accept_aln();
                            aln_vector.push_back(aln);
                        }
                        eof = alignment_file.eof();
                    }
                }

                // read processing
                for(auto aln : aln_vector) {
                    // dealloc if aln is not valid
                    if (preprocessing.is_valid(aln) == false)
                        continue;

                    // aln is ready for processing
                    const auto tag_encoded = preprocessing.tag_encoded;
					const auto orientation = preprocessing.orientation;
					const auto is_forward  = preprocessing.is_forward;

                    read_ptr_t read = std::make_shared<Read>\
                        (params, aln, aln_header, tag_encoded, orientation, is_forward);

//                    read->display_stats();
//                    read->display_read();

                    bool block_completed = !block.append_read(read);
                    if (block_completed == true) {
                        process_block(params, region_name, block, reference_file, output_str, statistics);
                        block.clear();
                    }
                }
                // finalize leftover
                process_block(params, region_name, block, reference_file, output_str, statistics);
                aln_vector.clear();
                block.clear();
            }

            #pragma omp critical (finalize_output)
            {
                output_string_vector[output_region_index] = output_str;
            }
        } // all regions done

        // update global statistics
        #pragma omp critical (finalize_statistics)
        {
            reference_file.update_statistics(global_statistics);
            alignment_file.update_statistics(global_statistics);
            global_statistics += statistics;
        }

    } // omp parallel

    // display output
    if (params.tmpvc || params.verbosity >= 2)
        for(const auto& str : output_string_vector)
            fwrite(str.data(), 1, str.size(), stdout);

    if (params.stats == true) {
        global_statistics.display_stats(global_metadata);
        global_statistics.display_crc32();
    }

    return EXIT_SUCCESS;
}
