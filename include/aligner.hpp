#pragma once

#include <arg_parse.hpp>
#include <file_to_map.hpp>

#include <thread>
#include <vector>
#include <iostream>

#include <seqan3/range/views/kmer_hash.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <seqan3/range/views/complement.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/split.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/sliding.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/distance.hpp>

#include <range/v3/algorithm/count.hpp>
#include <range/v3/algorithm/max.hpp>
#include <range/v3/view/set_algorithm.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>
#include <range/v3/algorithm/find_first_of.hpp>

#include <range/v3/algorithm/copy_backward.hpp>

void get_unique_kmers(  const cmd_arguments args, 
                        const std::vector<std::string> group_names,
                        const std::vector<std::size_t> group_scaffolds,
                        std::vector<std::vector<std::size_t>> & unique_kmers,
                        std::vector<std::size_t> & total_kmers);
int speq_run_1(cmd_arguments args);
//int speq_run_2(cmd_arguments args);
