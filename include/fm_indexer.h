#pragma once

#include <arg_parse.h>
#include <file_to_map.h>

#include <ThreadPool.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/parallel.hpp>
#include <seqan3/search/configuration/mode.hpp>
#include <seqan3/search/configuration/output.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/complement.hpp>

#include <seqan3/range/views/async_input_buffer.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/sliding.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/repeat_n.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/distance.hpp>


#include <vector>
#include <unordered_set>
#include <string>
#include <fstream>
#include <algorithm>
#include <future>
#include <tuple>



namespace speq
{
    namespace fm
    {
        int index(  speq::args::cmd_arguments & args);

        seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> generate_fm_index(speq::args::cmd_arguments & args);
        

        void async_count_unique_kmers_per_group(
            speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & index
        );

        void count_unique_kmers_per_group(  
            speq::args::cmd_arguments & args, 
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & index
        );

        void fast_count_unique_kmers_per_group(
            speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds
        );

    }
    
}
