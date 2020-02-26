#include "aligner.hpp"



int speq_run_1(cmd_arguments args)
{
    // Build the reference index


    // Organize the reference sequences by strain/variant type
    std::map<size_t,std::string> map_back = file_to_map(args.in_file_references_groups);

    // Input the reads
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};
    auto chunk_fin = fin | ranges::views::chunk(args.chunk);

    // FM-Index search read kmers against the 
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
                    seqan3::search_cfg::mode{seqan3::search_cfg::all} |
                    seqan3::search_cfg::parallel{args.threads};

    // For each chunk of the input read file
    for(auto it = std::begin(chunk_fin); it != std::end(chunk_fin); it++)
    {
        std::vector<std::vector<seqan3::dna5>> v_chunk;
        for(auto jt = std::begin(*it); jt != std::end(*it); jt++)
        {
            // Get all kmers for one sequence and it's reverse complement
            auto fwd_seq = seqan3::get<seqan3::field::seq>(*jt);
            auto rev_seq = fwd_seq | seqan3::views::complement | std::views::reverse;
            auto fwd_kmers = fwd_seq | ranges::views::sliding(args.kmer);
            for(auto k = ranges::begin(fwd_kmers); k != ranges::end(fwd_kmers); k++)
            {
                v_chunk.insert(v_chunk.end(), *k | ranges::to<std::vector>());
            }
            
            // I flip these kmers so when I am comparing search results I can iterate concurrently.
            auto rev_kmers = rev_seq | ranges::views::sliding(args.kmer) | std::views::reverse;
            for(auto k = ranges::begin(rev_kmers); k != ranges::end(rev_kmers); k++)
            {
                v_chunk.push_back(*k | ranges::to<std::vector>());
            }
        }
        seqan3::debug_stream << "Num Kmers in chunk: " << v_chunk.size() << "\n";
    }
    return 0;
}
/*
int speq_run_2(cmd_arguments args)
{
    seqan3::sequence_file_input fin1{args.in_file_reads_path_1};
    seqan3::sequence_file_input fin2{args.in_file_reads_path_2};

    auto fin = ranges::views::concat(fin1, fin2);

    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                    seqan3::search_cfg::mode{seqan3::search_cfg::all_best} |
                    seqan3::search_cfg::parallel{args.threads};

    // Duplicate the sequences, making the second one a reverse complement sequence.
    auto chunk_fin = fin | ranges::views::chunk(10000);
    /*
    auto both_ways = chunk_fin | ranges::views::for_each([](auto chunk)
    {
        auto post_one_chunk = chunk | ranges::views::for_each([](auto c)
        {
            auto rev = std::move(c) | seqan3::views::complement | std::views::reverse;
            auto fnr = ranges::view::zip(c,rev) | ranges::view::join;
        }

        
    }) | ranges::to<std::vector>();

    for(auto & x : both_ways)
    {
        seqan3::debug_stream << seqan3::get<seqan3::field::seq>(x) << "\n";
    }
    
    return 0;
}
*/