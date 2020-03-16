#pragma once

#include <seqan3/std/filesystem>
#include <seqan3/argument_parser/all.hpp>

namespace speq
{
    namespace args
    {
        struct cmd_arguments
        {
            std::filesystem::path in_file_reads_path_1{};
            std::filesystem::path in_file_reads_path_2{};
            std::filesystem::path in_file_references{};
            std::filesystem::path in_file_references_groups{};
            std::filesystem::path io_file_index{"output.idx"};
            std::filesystem::path out_file_path{"output.txt"};
            bool is_force{};
            unsigned int threads{std::thread::hardware_concurrency()};
            unsigned int chunk{10000};
            unsigned int kmer{30};
            bool is_indexer{false};
            bool is_scanner{false};
            bool is_parsed{false};
        };

        cmd_arguments initialize_argument_parser(int argc, char ** argv);
        cmd_arguments get_scan_arguments(seqan3::argument_parser & parser);
        cmd_arguments get_index_arguments(seqan3::argument_parser & parser);
        cmd_arguments get_all_arguments(seqan3::argument_parser & parser);
        
        void check_out_file(std::filesystem::path & a_path, const bool is_force);
        void check_in_file(std::filesystem::path & a_path);
    }
}

