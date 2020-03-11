#include "aligner.hpp"


using seqan3::operator""_dna5;


void get_unique_kmers(  const cmd_arguments args,
                        const std::vector<std::string> group_names,
                        const std::vector<int> group_scaffolds,
                        std::vector<std::vector<std::size_t>> & unique_kmers,
                        std::vector<std::size_t> & total_kmers_count)
{
    // Input the template scaffolds (genomes)
    seqan3::sequence_file_input fgenomes{args.in_file_references};
    // Kmer hash and unique sort each sequence
    std::vector<std::vector<std::vector<std::size_t>>> unique_hash_reference;
    unique_hash_reference.resize(group_names.size());
    std::vector<std::vector<std::vector<std::size_t>>> all_hash_reference;
    all_hash_reference.resize(group_names.size());
    std::size_t isolate_index = 0;
    for(auto & rec : fgenomes)
    {
        // If this scaffold was indicated for use in the groupings file...
        if(group_scaffolds[isolate_index] != -1)
        {
            // All N-containing sequences are arbitrary, and should be excluded
            // auto splitN = seqan3::get<seqan3::field::seq>(rec)  | ranges::views::split("N"_dna5) 
            //                                                     | ranges::to<std::vector>();
            // auto a_hash= ranges::views::join(ranges::views::transform(splitN, [args](auto c)
            // {
            //     return c    | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{args.kmer}})
            //                 | ranges::to<std::vector>();
            // })) | ranges::to<std::vector>();
            // auto b_hash= ranges::views::join(ranges::views::transform(splitN, [args](auto c)
            // {
            //     return seqan3::views::complement(std::views::reverse(c))
            //                 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{args.kmer}})
            //                 | ranges::to<std::vector>();
            // })) | ranges::to<std::vector>();

            auto a_hash = seqan3::get<seqan3::field::seq>(rec)  | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{args.kmer}})
                                                                | ranges::to<std::vector>();
            auto b_hash = seqan3::get<seqan3::field::seq>(rec)  | seqan3::views::complement
                                                                | std::views::reverse
                                                                | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{args.kmer}})
                                                                | ranges::to<std::vector>();
            a_hash = ranges::actions::sort(a_hash) | ranges::to<std::vector>();
            b_hash = ranges::actions::sort(b_hash) | ranges::to<std::vector>();
            auto gi_index = group_scaffolds[isolate_index];
            all_hash_reference[gi_index].push_back(a_hash);
            all_hash_reference[gi_index].push_back(b_hash);
            auto a_unique_hash = ranges::actions::unique(a_hash) | ranges::to<std::vector>();
            auto b_unique_hash = ranges::actions::unique(b_hash) | ranges::to<std::vector>();
            unique_hash_reference[gi_index].push_back(a_unique_hash);
            unique_hash_reference[gi_index].push_back(b_unique_hash);
        }
        ++isolate_index;
    }

    // For each strain-group
    std::vector<std::vector<std::size_t>> total_union_kmers;
    for(auto it = std::begin(all_hash_reference); it != std::end(all_hash_reference); ++it)
    {
        // Make a temp kmer vector
        std::vector<std::size_t> temp_kmer_vec;
        for(auto jt = std::begin(*it); jt != std::end(*it); ++jt)
        {
            // Concatenate the forward and reverse kmers
            auto forward = *jt++;
            auto reverse = *jt;
            auto concat_temp = ranges::views::concat(forward,reverse) | ranges::to<std::vector>();
            concat_temp = ranges::actions::sort(concat_temp) | ranges::to<std::vector>();
            // set_union with the current temp kmer vec
            temp_kmer_vec = ranges::views::set_union(temp_kmer_vec, concat_temp) | ranges::to<std::vector>();
        }
        // Add the distance of the temp kmer vector 
        total_union_kmers.push_back(temp_kmer_vec);
        total_kmers_count.push_back(ranges::distance(temp_kmer_vec));
    }

    seqan3::debug_stream << total_kmers_count << "\n";
    std::vector<std::vector<std::size_t>> single_unique_kmers;


    for(auto it = std::begin(unique_hash_reference); it != std::end(unique_hash_reference); ++it)
    {
        std::vector<std::size_t> temp_unique_kmer_vec;
        auto changer(*it);
        for(auto jt = std::begin(unique_hash_reference); jt != std::end(unique_hash_reference); ++jt)
        {
            if(it != jt)
            {
                // set_difference each hash_vector in each group against all other hash_vectors in **other** groups
                
                for(auto it2 = changer.begin(); it2 != changer.end(); ++it2)
                {
                    for(auto jt2 = std::begin(*jt); jt2 != std::end(*jt); ++jt2)
                    {
                        *it2 = ranges::views::set_difference(*it2, *jt2) | ranges::to<std::vector>();
                    }
                }
            }
        }
        // Set the union of the unique kmers for a given group into a single vector
        for(auto it2 = changer.begin(); it2 != changer.end(); ++it2)
        {
            temp_unique_kmer_vec = ranges::views::set_union(temp_unique_kmer_vec, *it2) | ranges::to<std::vector>();
        }
        // Push into the major vector
        //      Outer vector represents groups
        //          Inner vector represents sorted unique kmers
        single_unique_kmers.push_back(temp_unique_kmer_vec);
    }


    // Increase the number of each group-specific kmer to 
    // the max number actually found in any one scaffold within this group

    // for(std::size_t i = 0; i < single_unique_kmers.size(); ++i)
    // {
    //         seqan3::debug_stream << "guk: Preexpansion unique kmers for a group: " << single_unique_kmers[i].size() << "\n";
    //     single_unique_kmers[i] = ranges::views::for_each(single_unique_kmers[i], [all_hash_reference,i](auto c)
    //         {
    //             auto temp = all_hash_reference[i]   | ranges::views::transform([c](auto d){return ranges::count(d,c);});
    //             std::size_t num_repeats = ranges::max(temp);
    //             return ranges::yield_from(ranges::views::repeat_n(c, num_repeats)); 
    //         }
    //     ) | ranges::to<std::vector>();
    //     seqan3::debug_stream    << "guk: Increased to max kmers to " << single_unique_kmers[i].size() << "\n";
    // }

    unique_kmers.resize(single_unique_kmers.size());
    auto all_hash_it = std::begin(all_hash_reference);
    for(std::size_t i = 0; i < single_unique_kmers.size(); ++i)
    {
        seqan3::debug_stream << "guk: Preexpansion unique kmers for a group: " << single_unique_kmers[i].size() << "\n";
        unique_kmers[i].reserve(single_unique_kmers[i].size());
        for(auto jt = std::begin(single_unique_kmers[i]); jt != std::end(single_unique_kmers[i]); ++jt)
        {
            // Check each scaffold within each group to find the one with the most of this kmer
            std::size_t max_repeats = 0;
            for(auto kt = all_hash_it->begin(); kt != all_hash_it->end(); ++kt)
            {
                auto temp = std::lower_bound(kt->begin(), kt->end(), *jt);
                std::size_t this_max = std::count(temp, kt->end(), *jt);
                if(this_max > max_repeats) max_repeats = this_max;
            }
    
            if(max_repeats == 0)
            {
                seqan3::debug_stream << "WTF\t";
            }
            
            for(std::size_t k = 0; k < max_repeats; ++k)
            {
                unique_kmers[i].push_back(*jt);
            }
        }
        seqan3::debug_stream    << "guk: Increased to max kmers to " << unique_kmers[i].size() << "\n";
        ++all_hash_it;
    }

    return;
}

int speq_run_1(cmd_arguments args)
{
    seqan3::debug_stream << "Starting Run:\n";
    // Organize the reference sequences by strain/variant type
    std::vector<std::string> group_names;
    std::vector<int> group_scaffolds;
    file_to_map(args.in_file_references_groups, group_names, group_scaffolds);
    seqan3::debug_stream << "Mapped file.:" << group_names << " | " << group_scaffolds <<  "\n";

    // Find unique kmers in each scaffold/genome-set
    std::vector<std::vector<std::size_t>> unique_kmers;
    std::vector<std::size_t> total_kmers_count;
    get_unique_kmers(args, group_names, group_scaffolds, unique_kmers, total_kmers_count);
    seqan3::debug_stream << "Got Unique Kmers.\n";

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
    seqan3::debug_stream << hit_unique_kmers_per_group << "\n";
    seqan3::debug_stream << total_kmers_count << "\n";
    seqan3::debug_stream << ambiguous_reads << "\n";
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