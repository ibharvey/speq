#pragma once

#include <arg_parse.h>
#include <file_to_map.h>

#include <seqan3/core/debug_stream.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <vector>
#include <string>
#include <fstream>

namespace speq
{
    namespace kmer_hash
    {
        void get_unique_kmers(  const speq::args::cmd_arguments args, 
                            const std::vector<std::string> group_names,
                            const std::vector<int> group_scaffolds,
                            std::vector<std::vector<std::size_t>> & unique_kmers,
                            std::vector<std::size_t> & total_kmers);

        int index(  const speq::args::cmd_arguments args);
    }
    
}