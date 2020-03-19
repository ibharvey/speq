#include <fm_indexer.h>

// This function needs to be re-written to avoid 
// containers. Would dramatically help both
// CPU and memory overhead.
seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection>
speq::fm::generate_fm_index(speq::args::cmd_arguments & args)
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
    return index;
}

void speq::fm::async_count_unique_kmers_per_group(
            speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & index)
{
    // Double up the group_scaffolds to match the fm_index values
    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();
    
    //############# Counting Unique Kmers ##############
    // Get the reference sequences
    seqan3::sequence_file_input fin{args.in_file_references};
    auto fin_group = ranges::views::zip(group_scaffolds, fin);
    auto v = fin_group | seqan3::views::async_input_buffer(args.threads);
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
                seqan3::search_cfg::output{seqan3::search_cfg::text_position}|
                seqan3::search_cfg::mode{seqan3::search_cfg::all};
    auto worker = [=,&v]()
    {
        std::vector<std::size_t> a_total_kmers_per_group(group_names.size(),0);
        std::vector<std::size_t> a_unique_kmers_per_group(group_names.size(),0);
        auto search_and_count = [=, &a_total_kmers_per_group, &a_unique_kmers_per_group]
            (auto & kmers, auto & this_group)
            {
                auto results = search(kmers, index, config);
                for(auto rit = ranges::begin(results); rit != ranges::end(results); ++rit)
                {
                    auto result = *rit;
                    std::vector<std::size_t> hits_per_scaffold(double_group_scaffolds.size(),0);
                    for(auto & r : result)
                        hits_per_scaffold[r.first] += 1;
                    bool is_unique = true;
                    for(std::size_t hi = 0; hi < hits_per_scaffold.size(); ++hi)
                    {
                        if(hits_per_scaffold[hi] != 0)
                        {
                            if(this_group != double_group_scaffolds[hi])
                            {
                                is_unique = false;
                            }
                        }
                    }
                    ++a_total_kmers_per_group[this_group];
                    if(is_unique) ++a_unique_kmers_per_group[this_group];
                }
                return;
            };
        for(auto && [this_group, record] : v)
        {
            auto for_kmers = seqan3::get<seqan3::field::seq>(record)  | ranges::views::sliding(args.kmer); 
            search_and_count(for_kmers, this_group);
            auto rev_kmers = seqan3::get<seqan3::field::seq>(record)    | ranges::views::reverse
                                                                        | seqan3::views::deep{seqan3::views::complement}
                                                                        | ranges::views::sliding(args.kmer); 
            search_and_count(rev_kmers, this_group);
        }
        auto output = std::make_tuple(a_unique_kmers_per_group, a_total_kmers_per_group);
        return output;
    };

    std::vector<std::future<std::tuple<std::vector<std::size_t>, std::vector<std::size_t>>>> futures;
    for(std::size_t thread = 0; thread < args.threads; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

    std::vector<std::size_t> unique_kmers(group_names.size(),0);
    std::vector<std::size_t> total_kmers(group_names.size(),0);
    for(auto &f : futures)
    {
        auto output = f.get();
        auto a_unique = std::get<0>(output);
        auto a_total = std::get<1>(output);
        for(std::size_t i = 0; i < group_names.size(); ++i)
        {
            unique_kmers[i] += a_unique[i];
            total_kmers[i] += a_total[i];
        }
    }
    {
        auto holder = args.io_file_index.filename();
        std::string fn = args.io_file_index.stem();
        args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
        std::ofstream os{args.io_file_index, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(args.kmer);
        oarchive(group_names);
        oarchive(double_group_scaffolds);
        oarchive(unique_kmers);
        oarchive(total_kmers);
        args.io_file_index.replace_filename(holder);
    }
    seqan3::debug_stream << unique_kmers << "\n";
    seqan3::debug_stream << total_kmers << "\n";
    return;

}


void speq::fm::count_unique_kmers_per_group(
            speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & index)
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
    // Double up the group_scaffolds to match the fm_index values
    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();

    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
                seqan3::search_cfg::output{seqan3::search_cfg::text_position} |
                seqan3::search_cfg::mode{seqan3::search_cfg::all};
    
    std::size_t scaffold_counter = 0;
    std::vector<std::future<std::vector<std::vector<std::size_t>>>> pool_results;
    ThreadPool a_pool(args.threads);
    for(auto it = std::begin(fnr_ref_seqs); it != std::end(fnr_ref_seqs); ++it)
    {
        auto a_seq = *it;
        pool_results.emplace_back(
            a_pool.enqueue
            // I would think passing by reference would keep the memory costs of the program lower
            //([&]
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
                    
                        for(auto & res : results[ri])
                            hits_per_scaffold[res.first] += 1;
                        
                        // Check if only one group of scaffolds got hits
                        bool is_unique = true;
                        for(std::size_t hi = 0; hi < hits_per_scaffold.size(); ++hi)
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
    {
        auto holder = args.io_file_index.filename();
        std::string fn = args.io_file_index.stem();
        args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
        std::ofstream os{args.io_file_index, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(args.kmer);
        oarchive(group_names);
        oarchive(double_group_scaffolds);
        oarchive(unique_kmers_per_group);
        oarchive(total_kmers_per_group);
        args.io_file_index.replace_filename(holder);
    }
    seqan3::debug_stream << unique_kmers_per_group << "\n";
    seqan3::debug_stream << total_kmers_per_group << "\n";
    return;
}

/* 

TODO: The premise of the faster version is that the fm_search returns
    all hits for all scaffolds for each kmer. So we shouldn't have to
    search against a previous hit, just annotate it correctly the 
    first time you see it, and ignore all future instances.

Current bug: Failed assertion that after completing all kmers in a 
    scaffold, the ignore_vector should be empty.

TODO2: Filter the output index files to only include unique kmers,
    to save space and make the scanner faster.
*/
void speq::fm::fast_count_unique_kmers_per_group(
            speq::args::cmd_arguments & args,
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds)
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
    std::vector<std::unordered_set<std::size_t>> previously_hit_per_scaffold(double_group_scaffolds.size());
    std::vector<std::size_t> total_kmers_per_group(group_names.size(),0);
    std::vector<std::size_t> unique_kmers_per_group(group_names.size(),0);
    for(auto it = std::begin(fnr_ref_seqs); it != std::end(fnr_ref_seqs); ++it)
    {
        auto already_seen_hits = previously_hit_per_scaffold[scaffold_counter];

        auto ref_to_kmers = *it | ranges::views::sliding(args.kmer)
                                | ranges::to<std::vector>();

        // Remove all kmers from this search that were already found previously
        seqan3::debug_stream << "now you see me " << ref_to_kmers.size() << "\n";
        seqan3::debug_stream << "removing " << previously_hit_per_scaffold[scaffold_counter].size() << "\n";
        for(auto jt = previously_hit_per_scaffold[scaffold_counter].begin();
            jt != previously_hit_per_scaffold[scaffold_counter].end();
            ++jt)
        {
            ref_to_kmers.erase(ref_to_kmers.begin() + *jt);
        }
        seqan3::debug_stream << "now you don't " << ref_to_kmers.size() << "\n";
        previously_hit_per_scaffold[scaffold_counter].clear();
        // FM-index search
        auto config = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
                                                    seqan3::search_cfg::output{seqan3::search_cfg::text_position} |
                                                    seqan3::search_cfg::mode{seqan3::search_cfg::all} |
                                                    seqan3::search_cfg::parallel{args.threads};
        seqan3::debug_stream << "click " << double_group_scaffolds[scaffold_counter] << "\n";
        auto results = search(ref_to_kmers, index, config);
        // For each hit: determine if unique to one group
        //          If not, add to the previously_hit_per_scaffold vector for future reference seqs
        for(size_t ri = 0; ri < results.size(); ++ri)
        {
            // Check if this kmer was already hit in this result set
            auto check_it = std::find(  previously_hit_per_scaffold[scaffold_counter].begin(), 
                                        previously_hit_per_scaffold[scaffold_counter].end(),
                                        ri);
            if(check_it == previously_hit_per_scaffold[scaffold_counter].end())
            {
                std::vector<std::size_t> hits_per_scaffold(double_group_scaffolds.size(),0);
                for(auto & [idx, pos] : results[ri])
                {
                    hits_per_scaffold[idx] += 1;
                    // I'm going to count all instances of this hit from this first search hit
                    // Speed up the algorithm by not re-searching against known kmers
                    if(idx > scaffold_counter)
                    {
                        // If this kmer is found in a future scaffold
                        
                        previously_hit_per_scaffold[idx].insert(pos);
                        if(pos == 331)
                            seqan3::debug_stream    << "future: " 
                                                    << scaffold_counter << ", "
                                                    << idx << ", " 
                                                    << pos << "\n";
                    }
                    else if(idx == scaffold_counter)
                    {
                        // If the position found through search is different than the current position
                        if(pos > ri)
                        {
                            previously_hit_per_scaffold[idx].insert(pos);
                            if(pos == 331)
                            {
                                seqan3::debug_stream << results[ri] << "\n";
                                seqan3::debug_stream    << "now: " 
                                                        << scaffold_counter << ", "
                                                        << idx << ", " 
                                                        << pos << ", "
                                                        << ri  << "\n";
                            }
                                
                            //seqan3::debug_stream << "Add doublet" << ri << ", " << idx << ", " << pos << "\t";
                        }
                    }
                    
                }


                // Check if only one group of scaffolds got hits
                int first_group_index = -1;
                int count_group_hits = 0;
                for(std::size_t hi = 0; hi < hits_per_scaffold.size(); ++hi)
                {
                    total_kmers_per_group[double_group_scaffolds[hi]] += hits_per_scaffold[hi];
                    if(first_group_index != -2)
                    {
                        if(hits_per_scaffold[hi] != 0)
                        {
                            if(first_group_index == -1)
                            {
                                first_group_index = double_group_scaffolds[hi];
                                count_group_hits += hits_per_scaffold[hi];
                            }
                            else if(first_group_index == double_group_scaffolds[hi])
                            {
                                count_group_hits += hits_per_scaffold[hi];
                            }
                            else
                            {
                                // Hit 2 or more groups
                                // Not a unique kmer
                                first_group_index = -2;
                            }
                        }
                    }
                    
                }
                if(first_group_index >= 0)
                {
                    unique_kmers_per_group[first_group_index] += count_group_hits;
                }
            }
            else
            {
                // Remove this value from the previously hit scaffolds list
                previously_hit_per_scaffold[scaffold_counter].erase(check_it);
                //seqan3::debug_stream << "Removed: " << *check_it << "\n";
            }
        }
        // if(previously_hit_per_scaffold[scaffold_counter].size() != 0)
        // {
        //     seqan3::debug_stream << previously_hit_per_scaffold[scaffold_counter] << "\n";
        // }
        assert(previously_hit_per_scaffold[scaffold_counter].size() == 0);
        seqan3::debug_stream << "boom\n";
        scaffold_counter++;
    }
    //################## Cereal output other info ##################
    {
        auto holder = args.io_file_index.filename();
        std::string fn = args.io_file_index.stem();
        args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
        std::ofstream os{args.io_file_index, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(args.kmer);
        oarchive(group_names);
        oarchive(double_group_scaffolds);
        oarchive(unique_kmers_per_group);
        oarchive(total_kmers_per_group);
        args.io_file_index.replace_filename(holder);
    }
    seqan3::debug_stream << unique_kmers_per_group << "\n";
    seqan3::debug_stream << total_kmers_per_group << "\n";
    // Unnecessary given that args is not passed by reference
    // ....but comforting
    args.io_file_index.replace_extension(".idx");
    return;
}


int speq::fm::index( speq::args::cmd_arguments & args)
{
    // Organize the reference sequences by strain/variant type
    std::vector<int> group_scaffolds;
    std::vector<std::string> group_names;
    speq::file_to_map(args.in_file_references_groups, group_names, group_scaffolds);
    
    // If the index_file already exists and the user isn't forcing a re-write
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> index;
    if(std::filesystem::exists(args.io_file_index) && !args.is_force)
    {
        // Just input the index from the existing file
        {
            args.io_file_index.replace_extension(".idx");
            std::ifstream is{args.io_file_index, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index);
        }
    }else
    {
        // Generate the fm index from the sequences
        index = speq::fm::generate_fm_index(args);
    }
    auto holder = args.io_file_index.filename();
    std::string fn = args.io_file_index.stem();
    args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
    if(!std::filesystem::exists(args.io_file_index) || args.is_force)
    {
        args.io_file_index.replace_filename(holder);
        // Count unique kmers in each scaffold/genome-set
        speq::fm::async_count_unique_kmers_per_group(   
            args, 
            group_names,
            group_scaffolds,
            index);
    }
    else
    {
        // Don't need to remake kmer.dat files that already exist
        args.io_file_index.replace_filename(holder);
    }
    
    
    return 0;
}