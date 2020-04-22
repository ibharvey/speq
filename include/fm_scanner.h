#pragma once

#include <arg_parse.h>

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/chrono.hpp>
#include <cereal/types/string.hpp>

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
#include <seqan3/range/views/complement.hpp>

#include <seqan3/range/views/async_input_buffer.hpp>

#include <range/v3/view/transform.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/sliding.hpp>
#include <range/v3/distance.hpp>
#include <range/v3/begin_end.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/algorithm/min.hpp>


#include <vector>
#include <fstream>
#include <string>
#include <atomic>
#include <future>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <iostream>



namespace speq
{
    namespace scan
    {

        int one(        speq::args::cmd_arguments & args);
        int async_one(  speq::args::cmd_arguments & args);
        int _async_one_with_global_error_rate(  speq::args::cmd_arguments & args, double percent_perfect);
        int _async_one_with_local_error_rate(   speq::args::cmd_arguments & args);

        int two(        speq::args::cmd_arguments & args);
        int async_two(  speq::args::cmd_arguments & args);
        int _async_two_with_global_error_rate(  speq::args::cmd_arguments & args, double percent_perfect);
        int _async_two_with_local_error_rate(   speq::args::cmd_arguments & args);

        void _async_count_unique_kmers_per_group(
            speq::args::cmd_arguments & args,
            const std::filesystem::file_time_type & index_file_time,
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & index,
            std::vector<std::size_t> & unique_kmers,
            std::vector<std::size_t> & total_kmers        
        );
        
        std::vector<double> unique_to_percent(
            std::vector<double> unique_in_reads,
            std::size_t total_in_reads,
            std::vector<std::size_t> unique_in_refs,
            std::vector<std::size_t> total_in_refs
        );

        std::vector<double> unique_to_percent(
            std::vector<double> unique_in_reads,
            std::size_t total_in_reads,
            std::vector<double> unique_in_refs,
            std::vector<double> total_in_refs
        );

        std::vector<double> _async_one_global_estimate_kmer_per_group(
            const std::vector<double> & percent_per_group,
            const speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & double_group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & fm_index,
            const double percent_perfect
        );

        std::vector<double> _async_one_local_estimate_kmer_per_group(
            const std::vector<double> & percent_per_group,
            const speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & double_group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & fm_index
        );

        std::vector<double> _async_two_global_estimate_kmer_per_group(
            const std::vector<double> & percent_per_group,
            const speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & double_group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & fm_index,
            const double percent_perfect
        );

        std::vector<double> _async_two_local_estimate_kmer_per_group(
            const std::vector<double> & percent_per_group,
            const speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & double_group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & fm_index
        );

    }
}
