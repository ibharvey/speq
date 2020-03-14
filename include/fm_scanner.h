#pragma once

#include <arg_parse.h>

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

#include <seqan3/io/sequence_file/input.hpp>


namespace speq
{
    namespace scan
    {
        int one(    const speq::args::cmd_arguments args);
        int two(    const speq::args::cmd_arguments args);
    }
}