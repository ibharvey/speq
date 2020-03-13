#include <fm_indexer.h>

void speq::fm::count_unique_kmers_per_group(
            speq::args::cmd_arguments args,
            const std::vector<std::string> group_names,
            const std::vector<int> group_scaffolds)
{
    //################ Input the references for indexing ################
    // Pull all sequences from the reference lazily
    seqan3::sequence_file_input ref_recs{args.in_file_references};
    auto ref_seqs = ref_recs | std::views::transform([] (auto s) 
    { 
        return seqan3::get<seqan3::field::seq>(s);
    });
    // Get both the forward and reverse complements in the index file
    auto fnr_ref_seqs = ranges::view::for_each(ref_seqs,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();

    for(auto it = std::begin(fnr_ref_seqs); it < std::end(fnr_ref_seqs); ++it)
    {
        *(++it) = seqan3::views::complement(std::views::reverse(*it)) | ranges::to<std::vector>();
    }
    // Do the same doublign for group_scaffolds
    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();
    //################ Indexing ##################
    // Generate the fm_index against the sequences
    seqan3::fm_index index{fnr_ref_seqs};

    // Output the index for scanning later
    {
        args.io_file_index.replace_extension(".idx");
        std::ofstream os{args.io_file_index, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }
    //############# Counting Unique Kmers ##############
    // Check each sequence
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
                seqan3::search_cfg::output{seqan3::search_cfg::text_position} |
                seqan3::search_cfg::mode{seqan3::search_cfg::all} |
                seqan3::search_cfg::parallel{1};
    std::size_t scaffold_counter = 0;
    std::vector<std::future<std::vector<std::vector<std::size_t>>>> pool_results;
    ThreadPool a_pool(args.threads);
    for(auto it = std::begin(fnr_ref_seqs); it != std::end(fnr_ref_seqs); ++it)
    {
        auto a_seq = *it;
        pool_results.emplace_back(
            a_pool.enqueue
            ([a_seq, double_group_scaffolds, group_names, scaffold_counter, args, index, config]
                {
                    std::vector<std::size_t> a_total_kmers_per_group(group_names.size(),0);
                    std::vector<std::size_t> a_unique_kmers_per_group(group_names.size(),0);
                    auto this_group = double_group_scaffolds[scaffold_counter];
                    auto ref_to_kmers = a_seq   | ranges::views::sliding(args.kmer)
                                                | ranges::to<std::vector>();
                                                
                    auto results = search(ref_to_kmers, index, config);
                    for(size_t ri = 0; ri < results.size(); ++ri)
                    {
                        std::vector<std::size_t> hits_per_scaffold(double_group_scaffolds.size(),0);
                    
                        for(auto & [idx, pos] : results[ri])
                            hits_per_scaffold[idx] += 1;
                        
                        // Check if only one group of scaffolds got hits
                        bool is_unique = true;
                        for(int hi = 0; hi < hits_per_scaffold.size(); ++hi)
                        {
                            if(hits_per_scaffold[hi] != 0)
                            {
                                if(this_group != double_group_scaffolds[hi])
                                {
                                    is_unique = false;
                                }
                            }
                        
                        }
                        a_total_kmers_per_group[this_group]++;
                        if(is_unique)
                        {
                            a_unique_kmers_per_group[this_group]++;
                        }
                    }
                    std::vector<std::vector<std::size_t>> output;
                    output.push_back(a_total_kmers_per_group);
                    output.push_back(a_unique_kmers_per_group);
                    return output;
                }
            )
        );
        scaffold_counter++;
    }
    std::vector<std::size_t> total_kmers_per_group(group_names.size(),0);
    std::vector<std::size_t> unique_kmers_per_group(group_names.size(),0);
    for(auto && result : pool_results)
    {
        auto temp_result = result.get();
        auto temp_total = temp_result[0];
        auto temp_unique = temp_result[1];
        for(std::size_t i = 0; i < group_names.size(); ++i)
        {
            total_kmers_per_group[i] += temp_total[i];
            unique_kmers_per_group[i] += temp_unique[i];
        }
    }



    // std::vector<std::unordered_set<std::size_t>> previously_hit_per_scaffold(double_group_scaffolds.size());
    // for(auto it = std::begin(fnr_ref_seqs); it != std::end(fnr_ref_seqs); ++it)
    // {
    //     int this_group = double_group_scaffolds[scaffold_counter];
    //     auto already_seen_hits = previously_hit_per_scaffold[scaffold_counter];

    //     auto ref_to_kmers = *it | ranges::views::sliding(args.kmer)
    //                             | ranges::to<std::vector>();

    //     // Remove all kmers from this search that were already found previously
    //     seqan3::debug_stream << "now you see me " << ref_to_kmers.size() << "\n";
    //     seqan3::debug_stream << "removing " << previously_hit_per_scaffold[scaffold_counter].size() << "\n";
    //     for(auto jt = previously_hit_per_scaffold[scaffold_counter].begin();
    //         jt != previously_hit_per_scaffold[scaffold_counter].end();
    //         ++jt)
    //     {
    //         ref_to_kmers.erase(ref_to_kmers.begin() + *jt);
    //     }
    //     seqan3::debug_stream << "now you don't " << ref_to_kmers.size() << "\n";
    //     previously_hit_per_scaffold[scaffold_counter].clear();
    //     // FM-index search
    //     auto config = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
    //                                                 seqan3::search_cfg::output{seqan3::search_cfg::text_position} |
    //                                                 seqan3::search_cfg::mode{seqan3::search_cfg::all} |
    //                                                 seqan3::search_cfg::parallel{args.threads};
    //     seqan3::debug_stream << "click " << double_group_scaffolds[scaffold_counter] << "\n";
    //     auto results = search(ref_to_kmers, index, config);
    //     // For each hit: determine if unique to one group
    //     //          If not, add to the previously_hit_per_scaffold vector for future reference seqs
    //     for(size_t ri = 0; ri < results.size(); ++ri)
    //     {
    //         // Check if this kmer was already hit in this result set
    //         auto check_it = std::find(  previously_hit_per_scaffold[scaffold_counter].begin(), 
    //                                     previously_hit_per_scaffold[scaffold_counter].end(),
    //                                     ri);
    //         if(check_it == previously_hit_per_scaffold[scaffold_counter].end())
    //         {
    //             std::vector<std::size_t> hits_per_scaffold(double_group_scaffolds.size(),0);
    //             for(auto & [idx, pos] : results[ri])
    //             {
    //                 hits_per_scaffold[idx] += 1;
    //                 // I'm going to count all instances of this hit from this first search hit
    //                 // Speed up the algorithm by not re-searching against known kmers
    //                 if(idx > scaffold_counter)
    //                 {
    //                     // If this kmer is found in a future scaffold
                        
    //                     previously_hit_per_scaffold[idx].insert(pos);
    //                     if(pos == 331)
    //                         seqan3::debug_stream    << "future: " 
    //                                                 << scaffold_counter << ", "
    //                                                 << idx << ", " 
    //                                                 << pos << "\n";
    //                 }
    //                 else if(idx == scaffold_counter)
    //                 {
    //                     // If the position found through search is different than the current position
    //                     if(pos > ri)
    //                     {
    //                         previously_hit_per_scaffold[idx].insert(pos);
    //                         if(pos == 331)
    //                         {
    //                             seqan3::debug_stream << results[ri] << "\n";
    //                             seqan3::debug_stream    << "now: " 
    //                                                     << scaffold_counter << ", "
    //                                                     << idx << ", " 
    //                                                     << pos << ", "
    //                                                     << ri  << "\n";
    //                         }
                                
    //                         //seqan3::debug_stream << "Add doublet" << ri << ", " << idx << ", " << pos << "\t";
    //                     }
    //                 }
                    
    //             }


    //             // Check if only one group of scaffolds got hits
    //             int first_group_index = -1;
    //             int count_group_hits = 0;
    //             for(int hi = 0; hi < hits_per_scaffold.size(); ++hi)
    //             {
    //                 total_kmers_per_group[double_group_scaffolds[hi]] += hits_per_scaffold[hi];
    //                 if(first_group_index != -2)
    //                 {
    //                     if(hits_per_scaffold[hi] != 0)
    //                     {
    //                         if(first_group_index == -1)
    //                         {
    //                             first_group_index = double_group_scaffolds[hi];
    //                             count_group_hits += hits_per_scaffold[hi];
    //                         }
    //                         else if(first_group_index == double_group_scaffolds[hi])
    //                         {
    //                             count_group_hits += hits_per_scaffold[hi];
    //                         }
    //                         else
    //                         {
    //                             // Hit 2 or more groups
    //                             // Not a unique kmer
    //                             first_group_index = -2;
    //                         }
    //                     }
    //                 }
                    
    //             }
    //             if(first_group_index >= 0)
    //             {
    //                 unique_kmers_per_group[first_group_index] += count_group_hits;
    //             }
    //         }
    //         else
    //         {
    //             // Remove this value from the previously hit scaffolds list
    //             previously_hit_per_scaffold[scaffold_counter].erase(check_it);
    //             //seqan3::debug_stream << "Removed: " << *check_it << "\n";
    //         }
    //     }
    //     // if(previously_hit_per_scaffold[scaffold_counter].size() != 0)
    //     // {
    //     //     seqan3::debug_stream << previously_hit_per_scaffold[scaffold_counter] << "\n";
    //     // }
    //     assert(previously_hit_per_scaffold[scaffold_counter].size() == 0);
    //     seqan3::debug_stream << "boom\n";
    //     scaffold_counter++;
    // }
    //################## Cereal output other info ##################
    {
        args.io_file_index.replace_extension(".vec");
        std::ofstream os{args.io_file_index, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(   group_names, 
                    double_group_scaffolds, 
                    unique_kmers_per_group,
                    total_kmers_per_group);
    }
    seqan3::debug_stream << unique_kmers_per_group << "\n";
    seqan3::debug_stream << total_kmers_per_group << "\n";
    // Unnecessary given that args is not passed by reference
    // ....but comforting
    args.io_file_index.replace_extension(".idx");
    return;
}


int speq::fm::index(const speq::args::cmd_arguments args)
{
    // Organize the reference sequences by strain/variant type
    std::vector<int> group_scaffolds;
    std::vector<std::string> group_names;
    speq::file_to_map(args.in_file_references_groups, group_names, group_scaffolds);
    
    // Count unique kmers in each scaffold/genome-set
    speq::fm::count_unique_kmers_per_group(   
        args, 
        group_names,
        group_scaffolds);
    return 0;
}