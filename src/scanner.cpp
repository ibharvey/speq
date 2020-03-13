#include <scanner.h>


using seqan3::operator""_dna5;

void speq::load_index( const speq::args::cmd_arguments args)
{

}


int speq::scan_1(   const speq::args::cmd_arguments args)
{
    // Need total_kmers_count at the end of this for percent quantification
    // Need group_names to know how many groups I am parsing from the index file
    // Need unique_kmers, which **is** the kmr index file
    // Input the reads
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};
    auto chunk_fin = fin | ranges::views::chunk(args.chunk);

    // For each chunk of the input read file
    size_t ambiguous_reads = 0;
    std::vector<std::size_t> hit_unique_kmers_per_group;
    hit_unique_kmers_per_group.resize(group_names.size(),0);
    size_t count_out = 0;
    for(auto it = std::begin(chunk_fin); it != std::end(chunk_fin); it++)
    {
        for(auto jt = std::begin(*it); jt != std::end(*it); jt++)
        {
            // Sort the kmers within this read
            auto read_kmers = seqan3::get<seqan3::field::seq>(*jt)  | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{args.kmer}})
                                                                    | ranges::to<std::vector>();
            read_kmers = ranges::actions::sort(read_kmers);
            /* Compare these to the reference unique kmers
                There should only be one type with hits!
                If there are more than one group hit, this may be:
                    1) Recombination between two nucleotide sequences
                    2) Another species not accounted for in the reference file
            */
            bool already_hit = false; bool is_ambiguous = false;
            std::vector<std::size_t> hits;
            std::size_t group_index_counter = 0; std::size_t group_hit_index = 0;
            for(auto group_it = std::begin(unique_kmers); group_it != std::end(unique_kmers); ++group_it)
            {
                auto temp_intersection = ranges::views::set_intersection(read_kmers, *group_it);
                if(ranges::distance(temp_intersection) > 0)
                {
                    if(already_hit)
                    {
                        // The read is hitting at least two different reference species
                        // so we are calling it as ambiguous.
                        ambiguous_reads++;
                        is_ambiguous = true;
                        break;
                    }
                    else
                    {
                        hits = temp_intersection | ranges::to<std::vector>();
                        group_hit_index = group_index_counter;
                        already_hit = true;
                    }
                }
                group_index_counter++;
            }
            // If you hit one and only one group
            if(!is_ambiguous && already_hit)
            {
                // Add the number of unique kmers hit to that group's total
                // TODO: If you are only using 'hits' for the size, just keep a variable for the size
                hit_unique_kmers_per_group[group_hit_index] += ranges::distance(hits);
            }
        }
        seqan3::debug_stream << "\r" << ++count_out * args.chunk << " reads finished";
    }
    seqan3::debug_stream << "\n";
    seqan3::debug_stream << "Hits per group: " << hit_unique_kmers_per_group << "\n";
    auto unique_kmers_count = ranges::views::transform(unique_kmers, [] (auto c){return ranges::distance(c);})
                                | ranges::to<std::vector>();
    seqan3::debug_stream << "Unique Kmers per group: " << unique_kmers_count << "\n";
    seqan3::debug_stream << "Total Kmers per group: " << total_kmers_count << "\n";
    seqan3::debug_stream << "Abiguous reads: " << ambiguous_reads << "\n";
    return 0;
}
/*
int speq_run_2(speq::args::cmd_arguments args)
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

// std::vector<seqan3::dna5> text{"ACGTAGC"_dna5};
// auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}}) | ranges::to<std::vector>();
// hashes = ranges::actions::sort(hashes);
// seqan3::debug_stream << hashes << "\n";

// std::vector<seqan3::dna5> text2{"GCTACGTGCT"_dna5};
// auto hashes2 = text2 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}}) | ranges::to<std::vector>();
// hashes2 = ranges::actions::sort(hashes2);
// seqan3::debug_stream << hashes2 << "\n";

// auto hashes3 = ranges::actions::unique(hashes2);
// seqan3::debug_stream << hashes3 << "\n";

// auto hashes4 = ranges::views::set_difference(hashes3, hashes);
// seqan3::debug_stream << hashes4 << "\n";

// std::vector<seqan3::dna5> text25{"GCTACGTNGCT"_dna5};
// auto hashes25 = text25 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}}) | ranges::to<std::vector>();
// seqan3::debug_stream << hashes25 << "\n";
// hashes25 = ranges::actions::sort(hashes25);
// return 0;
