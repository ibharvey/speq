#pragma once

#include <seqan3/std/filesystem>
#include <seqan3/argument_parser/all.hpp>

struct cmd_arguments
{
    std::filesystem::path in_file_reads_path_1{};
    std::filesystem::path in_file_reads_path_2{};
    std::filesystem::path in_file_references{};
    std::filesystem::path in_file_references_groups{};
    std::filesystem::path out_file_path{"output.fa"};
    bool force{};
    u_int threads{std::thread::hardware_concurrency()};
    u_int chunk{10000};
    u_int8_t kmer{17};
};

cmd_arguments initialize_argument_parser(const std::string name, int argc, char ** argv);
void check_arguments(cmd_arguments & args);
void check_in_file(std::filesystem::path & a_path);