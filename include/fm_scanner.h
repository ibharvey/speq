#pragma once

#include <arg_parse.h>

#include <ThreadPool.h>

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/parallel.hpp>
#include <seqan3/search/configuration/mode.hpp>
#include <seqan3/search/configuration/output.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/range/views/zip.hpp>

#include <seqan3/range/views/async_input_buffer.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/sliding.hpp>
#include <range/v3/distance.hpp>
#include <range/v3/begin_end.hpp>
#include <range/v3/view/concat.hpp>

#include <vector>
#include <fstream>
#include <string>
#include <atomic>
#include <future>



namespace speq
{
    namespace scan
    {
        int one(        speq::args::cmd_arguments & args);
        int async_one(  speq::args::cmd_arguments & args);
        int two(        speq::args::cmd_arguments & args);
        int async_two(  speq::args::cmd_arguments & args);
    }
}