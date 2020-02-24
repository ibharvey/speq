#pragma once

#include <arg_parse.hpp>

#include <thread>
#include <vector>

#include <seqan3/search/configuration/all.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <seqan3/range/views/complement.hpp>
#include <range/v3/view/repeat.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/repeat_n.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/view/join.hpp>

#include <range/v3/algorithm/copy_backward.hpp>


int speq_run_1(cmd_arguments args);
int speq_run_2(cmd_arguments args);
