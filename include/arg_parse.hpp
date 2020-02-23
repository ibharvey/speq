#pragma once

#include <seqan3/std/filesystem>
#include <seqan3/argument_parser/all.hpp>

struct cmd_arguments
{
    std::filesystem::path in_file_path{};
    std::filesystem::path out_file_path{"output.fa"};
    bool reverse{};
    bool overwrite{};
    bool force{};
    bool use_gzip{};
    bool use_bzip{};
    u_int threads{std::thread::hardware_concurrency()};
};

cmd_arguments initialize_argument_parser(const std::string name, int argc, char ** argv);
void check_arguments(cmd_arguments & args);