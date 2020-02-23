#pragma once

#include <seqan3/std/filesystem>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/std/filesystem>

#include <iostream>
#include <string>

struct cmd_arguments
{
    std::filesystem::path in_file_path{};
    std::filesystem::path out_file_path{"output.fa"};
    bool reverse{};
    bool use_gzip{};
    bool use_bzip{};
    u_int threads{std::thread::hardware_concurrency()};
};

cmd_arguments initialize_argument_parser(const std::string name, int argc, char ** argv);
int check_arguments(cmd_arguments & args, int & err_handle);