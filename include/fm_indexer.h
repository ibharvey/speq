#pragma once

#include <arg_parse.h>
#include <file_to_map.h>

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/chrono.hpp>
#include <cereal/types/string.hpp>

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
#include <chrono>



namespace speq
{
    namespace fm
    {
        int index(  speq::args::cmd_arguments & args);

        
        void generate_fm_index(
            speq::args::cmd_arguments & args,
            const std::filesystem::file_time_type & scaffold_file_time,
            const std::filesystem::file_time_type & group_file_time,
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds
        );
        
    }
    
}
