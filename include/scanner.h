#pragma once

#include <arg_parse.h>
#include <file_to_map.h>

#include <thread>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

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
#include <range/v3/view/concat.hpp>
#include <range/v3/distance.hpp>
#include <range/v3/view/set_algorithm.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>


namespace speq
{
    void load_index( const speq::args::cmd_arguments args);
    int scan_1( const speq::args::cmd_arguments args);

    int scan_2(const speq::args::cmd_arguments args);
}

