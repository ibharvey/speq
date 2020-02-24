#include "aligner.hpp"

int speq_run_1(cmd_arguments args)
{
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};

    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                    seqan3::search_cfg::mode{seqan3::search_cfg::all_best} |
                    seqan3::search_cfg::parallel{args.threads};

    // Duplicate the sequences, making the second one a reverse complement sequence.
    auto both_ways = fin | ranges::views::for_each([](auto c)
    {
        auto rev = std::move(c) | seqan3::views::complement | std::views::reverse;
        return ranges::view::zip(c,rev) | ranges::view::join;
    }) | ranges::to<std::vector>();

    for(auto & x : both_ways)
    {
        seqan3::debug_stream << seqan3::get<seqan3::field::seq>(x) << "\n";
    }
    
    return 0;
}

int speq_run_2(cmd_arguments args)
{
    seqan3::sequence_file_input fin1{args.in_file_reads_path_1};
    seqan3::sequence_file_input fin2{args.in_file_reads_path_2};

    auto fin = ranges::views::zip(fin1, fin2) | ranges::view::join;

    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                    seqan3::search_cfg::mode{seqan3::search_cfg::all_best} |
                    seqan3::search_cfg::parallel{args.threads};

    // Duplicate the sequences, making the second one a reverse complement sequence.
    auto both_ways = fin | ranges::views::for_each([](auto c)
    {
        auto rev = std::move(c) | seqan3::views::complement | std::views::reverse;
        auto fnr = ranges::view::zip(c,rev) | ranges::view::join;
        
    }) | ranges::to<std::vector>();

    for(auto & x : both_ways)
    {
        seqan3::debug_stream << seqan3::get<seqan3::field::seq>(x) << "\n";
    }
    
    return 0;
}