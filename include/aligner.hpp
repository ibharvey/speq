#pragma once

#include <arg_parse.hpp>

#include <thread>
#include <vector>
#include <iostream>

#include <seqan3/search/configuration/all.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <seqan3/range/views/complement.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/sliding.hpp>

#include <range/v3/algorithm/copy_backward.hpp>


int speq_run_1(cmd_arguments args);
//int speq_run_2(cmd_arguments args);
