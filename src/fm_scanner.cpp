#include <fm_scanner.h>

int speq::scan::async_one(    speq::args::cmd_arguments & args, double percent_perfect)
{
    // Choose whether to use the global or read-specific error rates
    if(percent_perfect == 0.0)
    {
        return speq::scan::_async_one_with_local_error_rate(args);
    }
    else
    {
        return speq::scan::_async_one_with_global_error_rate(args, percent_perfect);
    }
}

int speq::scan::async_two(  speq::args::cmd_arguments & args, double percent_perfect)
{
    if(percent_perfect == 0.0)
    {
        return speq::scan::_async_two_with_local_error_rate(args);
    }
    else
    {
        return speq::scan::_async_two_with_global_error_rate(args, percent_perfect);
    }
}


int speq::scan::_async_one_with_global_error_rate(  speq::args::cmd_arguments & args, double percent_perfect)
{
    // Import the index sequence file (.idx)
    std::vector<std::string> group_names{};
    std::vector<int> group_scaffolds{};
    std::filesystem::path stored_scaffold_file_name;
    std::filesystem::file_time_type scaffold_file_time;
    std::filesystem::file_time_type group_file_time;
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> fm_index;
    {
        args.io_file_index.replace_extension(".idx");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        std::string stored_scaffold_file_name_str;
        iarchive(stored_scaffold_file_name_str);
        stored_scaffold_file_name = std::filesystem::path(stored_scaffold_file_name_str);
        iarchive(scaffold_file_time);
        iarchive(group_file_time);
        iarchive(group_names);
        iarchive(group_scaffolds);
        iarchive(fm_index);
    }

    if(args.in_file_references != stored_scaffold_file_name)
    {
        // Ignore the user input because a different file was used for constructing the index file
        if(std::filesystem::exists(stored_scaffold_file_name))
        {
            args.in_file_references = stored_scaffold_file_name;
        }
        else
        {
            throw std::logic_error{"SPeQ could not find the sequence file used to make the specified index file."};
        }
        
    }

    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();

    std::filesystem::file_time_type old_index_file_time = std::filesystem::last_write_time(args.io_file_index);
    auto holder = args.io_file_index.filename();
    std::string fn = args.io_file_index.stem();
    args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
    auto second_holder = args.io_file_index.filename();
    // Import or build the index data file (.dat)
    std::vector<std::size_t> total_kmers_per_group(group_names.size(),0);
    std::vector<std::size_t> unique_kmers_per_group(group_names.size(),0);
    std::filesystem::file_time_type index_file_time;

    if(std::filesystem::exists(args.io_file_index))
    {
        {
            std::ifstream is{args.io_file_index, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index_file_time);
        }
        args.io_file_index.replace_filename(holder);
        std::time_t dftime = decltype(index_file_time)::clock::to_time_t(index_file_time);
        if(old_index_file_time != index_file_time)
        {
            speq::scan::_async_count_unique_kmers_per_group(
                args,
                old_index_file_time,
                group_names,
                group_scaffolds,
                fm_index,
                unique_kmers_per_group,
                total_kmers_per_group
            );
        }
        else
        {
            {
                args.io_file_index.replace_filename(second_holder);
                std::ifstream is{args.io_file_index, std::ios::binary};
                cereal::BinaryInputArchive iarchive{is};
                iarchive(index_file_time);
                iarchive(unique_kmers_per_group);
                iarchive(total_kmers_per_group);
                args.io_file_index.replace_filename(holder);
            }
        }
        
    }
    else
    {
        args.io_file_index.replace_filename(holder);
        speq::scan::_async_count_unique_kmers_per_group(
            args,
            old_index_file_time,
            group_names,
            group_scaffolds,
            fm_index,
            unique_kmers_per_group,
            total_kmers_per_group
        );
    }

    // Import the reads
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
            seqan3::search_cfg::output{seqan3::search_cfg::text_position};
    auto v = fin | seqan3::views::async_input_buffer(args.threads * 2);

    std::atomic<std::size_t> ambiguous_reads = 0;
    std::atomic<std::size_t> total_kmers = 0;
    auto worker = [&v, &ambiguous_reads, &total_kmers, fm_index, 
                    group_names, double_group_scaffolds, config, args] ()
    {
        std::vector<std::size_t> unique_kmers(group_names.size(),0);
        for(auto & record : v)
        {
            bool is_ambiguous = false;
            int which_group = -1;
            
            auto do_a_count = [&](auto & results, auto & qmers)
            {
                auto q_it = ranges::begin(qmers);
                for(auto r_it = ranges::begin(results); r_it != ranges::end(results); ++r_it)
                {
                    auto result = *r_it;
                    auto qesult = *q_it;
                    if(ranges::min(qesult) > args.phred_cutoff)
                    {
                        ++total_kmers;
                        int which_hit  = -1;
                        for(auto & res : result)
                        {
                            if(which_hit == -1)
                            {
                                which_hit = double_group_scaffolds[res.first];
                            }
                            else if(which_hit != double_group_scaffolds[res.first])
                            {
                                which_hit = -2;
                                break;
                            }
                        }
                        if(which_hit >= 0) 
                        {
                            ++unique_kmers[which_hit];
                            // For now, still count unique kmers on an ambiguous read
                            //      because I don't have a justification to ignore it.
                            if(which_group >= 0 && which_group != which_hit)
                            {
                                is_ambiguous = true;
                            }
                            else
                            {
                                which_group = which_hit;
                            }
                        }  
                    }
                    ++q_it; 
                }
            };
            // Pull the sequences and quality scores
            auto f_seq = seqan3::get<seqan3::field::seq>(record);
            auto f_qual = seqan3::get<seqan3::field::qual>(record);
            // auto rc_seq = seqan3::get<seqan3::field::seq>(record)   | ranges::views::reverse | seqan3::views::deep{seqan3::views::complement};
            // auto rc_qual = seqan3::get<seqan3::field::qual>(record)   | ranges::views::reverse;
            // Produce kmers of the sequences and quality scores
            auto f_kmers = f_seq | ranges::views::sliding(args.kmer);
            auto f_qmers = f_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            // auto rc_kmers = rc_seq | ranges::views::sliding(args.kmer);
            // auto rc_qmers = rc_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            // Search against the kmers
            auto f_results = search(f_kmers, fm_index, config);
            // auto rc_results = search(rc_kmers, fm_index, config);
            // Count kmers based on kmer-specific quality scores
            do_a_count(f_results, f_qmers);
            // do_a_count(rc_results, rc_qmers);
            if(is_ambiguous) ++ambiguous_reads;
        }
        return unique_kmers;
    };

    std::vector<std::future<std::vector<std::size_t>>> futures;
    for(std::size_t thread = 0; thread < args.threads - 1; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

    std::vector<std::size_t> unique_totals_int(group_names.size(),0);
    for(auto &e : futures)
    {
        auto a_unique = e.get();
        for(std::size_t i = 0; i < a_unique.size(); ++i)
        {
            unique_totals_int[i] += a_unique[i];
            
        }
    }
    std::vector<double> unique_totals(group_names.size(),0.0);
    for(std::size_t i = 0; i < group_names.size(); ++i)
        unique_totals[i] += static_cast<double>(unique_totals_int[i]) / percent_perfect;
    seqan3::debug_stream << unique_totals << "\n";
    seqan3::debug_stream << total_kmers << "\t" << ambiguous_reads << "\n";
    return 1;
}

int speq::scan::_async_one_with_local_error_rate(   speq::args::cmd_arguments & args)
{
    // Import the index sequence file (.idx)
    std::vector<std::string> group_names{};
    std::vector<int> group_scaffolds{};
    std::filesystem::path stored_scaffold_file_name;
    std::filesystem::file_time_type scaffold_file_time;
    std::filesystem::file_time_type group_file_time;
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> fm_index;
    {
        args.io_file_index.replace_extension(".idx");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        std::string stored_scaffold_file_name_str;
        iarchive(stored_scaffold_file_name_str);
        stored_scaffold_file_name = std::filesystem::path(stored_scaffold_file_name_str);
        iarchive(scaffold_file_time);
        iarchive(group_file_time);
        iarchive(group_names);
        iarchive(group_scaffolds);
        iarchive(fm_index);
    }
    
    if(args.in_file_references != stored_scaffold_file_name)
    {
        // Ignore the user input because a different file was used for constructing the index file
        if(std::filesystem::exists(stored_scaffold_file_name))
        {
            args.in_file_references = stored_scaffold_file_name;
        }
        else
        {
            throw std::logic_error{"SPeQ could not find the sequence file used to make the specified index file."};
        }
        
    }

    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();

    std::filesystem::file_time_type old_index_file_time = std::filesystem::last_write_time(args.io_file_index);
    auto holder = args.io_file_index.filename();
    std::string fn = args.io_file_index.stem();
    args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
    auto second_holder = args.io_file_index.filename();
    // Import or build the index data file (.dat)
    std::vector<std::size_t> total_kmers_per_group(group_names.size(),0);
    std::vector<std::size_t> unique_kmers_per_group(group_names.size(),0);
    std::filesystem::file_time_type index_file_time;

    if(std::filesystem::exists(args.io_file_index))
    {
        {
            std::ifstream is{args.io_file_index, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index_file_time);
        }
        args.io_file_index.replace_filename(holder);
        std::time_t dftime = decltype(index_file_time)::clock::to_time_t(index_file_time);
        if(old_index_file_time != index_file_time)
        {
            speq::scan::_async_count_unique_kmers_per_group(
                args,
                old_index_file_time,
                group_names,
                group_scaffolds,
                fm_index,
                unique_kmers_per_group,
                total_kmers_per_group
            );
        }
        else
        {
            {
                args.io_file_index.replace_filename(second_holder);
                std::ifstream is{args.io_file_index, std::ios::binary};
                cereal::BinaryInputArchive iarchive{is};
                iarchive(index_file_time);
                iarchive(unique_kmers_per_group);
                iarchive(total_kmers_per_group);
                args.io_file_index.replace_filename(holder);
            }
        }
        
    }
    else
    {
        args.io_file_index.replace_filename(holder);
        speq::scan::_async_count_unique_kmers_per_group(
            args,
            old_index_file_time,
            group_names,
            group_scaffolds,
            fm_index,
            unique_kmers_per_group,
            total_kmers_per_group
        );
    }

    // Import the reads
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
            seqan3::search_cfg::output{seqan3::search_cfg::text_position};
    auto v = fin | seqan3::views::async_input_buffer(args.threads * 2);

    std::atomic<std::size_t> ambiguous_reads = 0;
    std::atomic<std::size_t> total_kmers = 0;
    auto worker = [&v, &ambiguous_reads, &total_kmers, fm_index, 
                    group_names, double_group_scaffolds, config, args] ()
    {
        std::vector<double> unique_kmers(group_names.size(),0);
        for(auto & record : v)
        {
            bool is_ambiguous = false;
            int which_group = -1;
            
            auto do_a_count = [&](auto & results, auto & qmers)
            {
                auto q_it = ranges::begin(qmers);
                for(auto r_it = ranges::begin(results); r_it != ranges::end(results); ++r_it)
                {
                    auto result = *r_it;
                    auto qesult = *q_it;
                    if(ranges::min(qesult) > args.phred_cutoff)
                    {
                        ++total_kmers;
                        int which_hit  = -1;
                        for(auto & res : result)
                        {
                            if(which_hit == -1)
                            {
                                which_hit = double_group_scaffolds[res.first];
                            }
                            else if(which_hit != double_group_scaffolds[res.first])
                            {
                                which_hit = -2;
                                break;
                            }
                        }
                        if(which_hit >= 0) 
                        {
                            // Average the Phred scores to get a probability that this kmer is perfect
                            auto qavg = ranges::accumulate(qesult, 1.0, [](double a, double b){return a / (1.0 - 1.0/pow(10.0,b/10.0));});
                            unique_kmers[which_hit] += qavg;
                            // For now, still count unique kmers on an ambiguous read
                            //      because I don't have a justification to ignore it.
                            if(which_group >= 0 && which_group != which_hit)
                            {
                                is_ambiguous = true;
                            }
                            else
                            {
                                which_group = which_hit;
                            }
                        }  
                    }
                    ++q_it; 
                }
            };
            // Pull the sequences and quality scores
            auto f_seq = seqan3::get<seqan3::field::seq>(record);
            auto f_qual = seqan3::get<seqan3::field::qual>(record);
            // auto rc_seq = seqan3::get<seqan3::field::seq>(record)   | ranges::views::reverse | seqan3::views::deep{seqan3::views::complement};
            // auto rc_qual = seqan3::get<seqan3::field::qual>(record)   | ranges::views::reverse;
            // Produce kmers of the sequences and quality scores
            auto f_kmers = f_seq | ranges::views::sliding(args.kmer);
            auto f_qmers = f_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);

            // auto rc_kmers = rc_seq | ranges::views::sliding(args.kmer);
            // auto rc_qmers = rc_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            // Search against the kmers
            auto f_results = search(f_kmers, fm_index, config);
            // auto rc_results = search(rc_kmers, fm_index, config);
            // Count kmers based on kmer-specific quality scores
            do_a_count(f_results, f_qmers);
            // do_a_count(rc_results, rc_qmers);
            if(is_ambiguous) ++ambiguous_reads;
        }
        return unique_kmers;
    };

    std::vector<std::future<std::vector<double>>> futures;
    for(std::size_t thread = 0; thread < args.threads - 1; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

    std::vector<double> unique_totals(group_names.size(),0);
    for(auto &e : futures)
    {
        auto a_unique = e.get();
        for(std::size_t i = 0; i < a_unique.size(); ++i)
        {
            unique_totals[i] += a_unique[i];
        }
    }
    seqan3::debug_stream << unique_totals << "\n";
    seqan3::debug_stream << total_kmers << "\t" << ambiguous_reads << "\n";
    return 1;
}

int speq::scan::_async_two_with_global_error_rate(    speq::args::cmd_arguments & args, double percent_perfect)
{
    // Import the index sequence file (.idx)
    std::vector<std::string> group_names{};
    std::vector<int> group_scaffolds{};
    std::filesystem::path stored_scaffold_file_name;
    std::filesystem::file_time_type scaffold_file_time;
    std::filesystem::file_time_type group_file_time;
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> fm_index;
    {
        args.io_file_index.replace_extension(".idx");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        std::string stored_scaffold_file_name_str;
        iarchive(stored_scaffold_file_name_str);
        stored_scaffold_file_name = std::filesystem::path(stored_scaffold_file_name_str);
        iarchive(scaffold_file_time);
        iarchive(group_file_time);
        iarchive(group_names);
        iarchive(group_scaffolds);
        iarchive(fm_index);
    }

    if(args.in_file_references != stored_scaffold_file_name)
    {
        // Ignore the user input because a different file was used for constructing the index file
        if(std::filesystem::exists(stored_scaffold_file_name))
        {
            args.in_file_references = stored_scaffold_file_name;
        }
        else
        {
            throw std::logic_error{"SPeQ could not find the sequence file used to make the specified index file."};
        }
        
    }

    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();

    std::filesystem::file_time_type old_index_file_time = std::filesystem::last_write_time(args.io_file_index);
    auto holder = args.io_file_index.filename();
    std::string fn = args.io_file_index.stem();
    args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
    // Import or build the index data file (.dat)
    std::vector<std::size_t> total_kmers_per_group(group_names.size(),0);
    std::vector<std::size_t> unique_kmers_per_group(group_names.size(),0);
    std::filesystem::file_time_type index_file_time;
    if(std::filesystem::exists(args.io_file_index))
    {
        {
            std::ifstream is{args.io_file_index, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index_file_time);
            iarchive(unique_kmers_per_group);
            iarchive(total_kmers_per_group);
        }
        if(old_index_file_time != index_file_time)
        {
            speq::scan::_async_count_unique_kmers_per_group(
                args,
                index_file_time,
                group_names,
                group_scaffolds,
                fm_index,
                unique_kmers_per_group,
                total_kmers_per_group
            );
        }
    }
    else
    {
        speq::scan::_async_count_unique_kmers_per_group(
            args,
            index_file_time,
            group_names,
            group_scaffolds,
            fm_index,
            unique_kmers_per_group,
            total_kmers_per_group
        );
    }
    args.io_file_index.replace_filename(holder);

    // Import the reads
    seqan3::sequence_file_input fin1{args.in_file_reads_path_1};
    seqan3::sequence_file_input fin2{args.in_file_reads_path_2};

    auto combined = seqan3::views::zip(fin1, fin2);
    auto v = combined | seqan3::views::async_input_buffer(args.threads * 2);

    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
            seqan3::search_cfg::output{seqan3::search_cfg::text_position};

    std::atomic<std::size_t> ambiguous_reads = 0;
    std::atomic<std::size_t> total_kmers = 0;
    auto worker = [&v, &ambiguous_reads, &total_kmers, fm_index, 
                    group_names, double_group_scaffolds, config, args] ()
    {
        std::vector<std::size_t> unique_kmers(group_names.size(),0);
        for(auto && [rec1, rec2] : v)
        {
            bool is_ambiguous = false;
            int which_group = -1;
            auto seq1 = seqan3::get<seqan3::field::seq>(rec1);
            auto kmers1 = seq1 | ranges::views::sliding(args.kmer);
            auto seq2 = seqan3::get<seqan3::field::seq>(rec2);
            auto kmers2 = seq2 | ranges::views::sliding(args.kmer);
            auto all_kmers = ranges::views::concat(kmers1, kmers2);
            auto results = search(all_kmers, fm_index, config);
            for(auto r_it = ranges::begin(results); r_it != ranges::end(results); ++r_it)
            {
                ++total_kmers;
                auto result = *r_it;
                int which_hit  = -1;
                for(auto & res : result)
                {
                    if(which_hit == -1)
                    {
                        which_hit = double_group_scaffolds[res.first];
                    }
                    else if(which_hit != double_group_scaffolds[res.first])
                    {
                        which_hit = -2;
                        break;
                    }
                }
                if(which_hit >= 0) 
                {
                    ++unique_kmers[which_hit];
                    // For now, still count unique kmers on an ambiguous read
                    //      because I don't have a justification to ignore it.
                    if(which_group >= 0 && which_group != which_hit)
                    {
                        is_ambiguous = true;
                    }
                    else
                    {
                        which_group = which_hit;
                    }
                }   
            }
            if(is_ambiguous) ++ambiguous_reads;
        }
        return unique_kmers;
    };

    std::vector<std::future<std::vector<std::size_t>>> futures;
    for(std::size_t thread = 0; thread < args.threads; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

    std::vector<std::size_t> unique_totals(group_names.size(),0);
    for(auto &e : futures)
    {
        auto a_unique = e.get();
        for(std::size_t i = 0; i < a_unique.size(); ++i)
        {
            unique_totals[i] += a_unique[i];
        }
    }

    seqan3::debug_stream << unique_totals << "\n";
    seqan3::debug_stream << total_kmers << "\t" << ambiguous_reads << "\n";
    return 1;
}

int speq::scan::_async_two_with_local_error_rate(    speq::args::cmd_arguments & args)
{
    // Import the index sequence file (.idx)
    std::vector<std::string> group_names{};
    std::vector<int> group_scaffolds{};
    std::filesystem::path stored_scaffold_file_name;
    std::filesystem::file_time_type scaffold_file_time;
    std::filesystem::file_time_type group_file_time;
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> fm_index;
    {
        args.io_file_index.replace_extension(".idx");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        std::string stored_scaffold_file_name_str;
        iarchive(stored_scaffold_file_name_str);
        stored_scaffold_file_name = std::filesystem::path(stored_scaffold_file_name_str);
        iarchive(scaffold_file_time);
        iarchive(group_file_time);
        iarchive(group_names);
        iarchive(group_scaffolds);
        iarchive(fm_index);
    }

    if(args.in_file_references != stored_scaffold_file_name)
    {
        // Ignore the user input because a different file was used for constructing the index file
        if(std::filesystem::exists(stored_scaffold_file_name))
        {
            args.in_file_references = stored_scaffold_file_name;
        }
        else
        {
            throw std::logic_error{"SPeQ could not find the sequence file used to make the specified index file."};
        }
        
    }

    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();

    std::filesystem::file_time_type old_index_file_time = std::filesystem::last_write_time(args.io_file_index);
    auto holder = args.io_file_index.filename();
    std::string fn = args.io_file_index.stem();
    args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
    auto second_holder = args.io_file_index.filename();
    // Import or build the index data file (.dat)
    std::vector<std::size_t> total_kmers_per_group(group_names.size(),0);
    std::vector<std::size_t> unique_kmers_per_group(group_names.size(),0);
    std::filesystem::file_time_type index_file_time;

    if(std::filesystem::exists(args.io_file_index))
    {
        {
            std::ifstream is{args.io_file_index, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index_file_time);
        }
        args.io_file_index.replace_filename(holder);
        std::time_t dftime = decltype(index_file_time)::clock::to_time_t(index_file_time);
        if(old_index_file_time != index_file_time)
        {
            speq::scan::_async_count_unique_kmers_per_group(
                args,
                old_index_file_time,
                group_names,
                group_scaffolds,
                fm_index,
                unique_kmers_per_group,
                total_kmers_per_group
            );
        }
        else
        {
            {
                args.io_file_index.replace_filename(second_holder);
                std::ifstream is{args.io_file_index, std::ios::binary};
                cereal::BinaryInputArchive iarchive{is};
                iarchive(index_file_time);
                iarchive(unique_kmers_per_group);
                iarchive(total_kmers_per_group);
                args.io_file_index.replace_filename(holder);
            }
        }
        
    }
    else
    {
        args.io_file_index.replace_filename(holder);
        speq::scan::_async_count_unique_kmers_per_group(
            args,
            old_index_file_time,
            group_names,
            group_scaffolds,
            fm_index,
            unique_kmers_per_group,
            total_kmers_per_group
        );
    }

    // Import the reads
    seqan3::sequence_file_input fin1{args.in_file_reads_path_1};
    seqan3::sequence_file_input fin2{args.in_file_reads_path_2};

    auto combined = seqan3::views::zip(fin1, fin2);
    auto v = combined | seqan3::views::async_input_buffer(args.threads * 2);

    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
            seqan3::search_cfg::output{seqan3::search_cfg::text_position};

    std::atomic<std::size_t> ambiguous_reads = 0;
    std::atomic<std::size_t> total_kmers = 0;
    auto worker = [&v, &ambiguous_reads, &total_kmers, fm_index, 
                    group_names, double_group_scaffolds, config, args] ()
    {
        std::vector<std::size_t> unique_kmers(group_names.size(),0);
        for(auto && [rec1, rec2] : v)
        {
            bool is_ambiguous = false;
            int which_group = -1;
            auto seq1 = seqan3::get<seqan3::field::seq>(rec1);
            auto kmers1 = seq1 | ranges::views::sliding(args.kmer);
            auto seq2 = seqan3::get<seqan3::field::seq>(rec2);
            auto kmers2 = seq2 | ranges::views::sliding(args.kmer);
            auto all_kmers = ranges::views::concat(kmers1, kmers2);
            auto results = search(all_kmers, fm_index, config);
            for(auto r_it = ranges::begin(results); r_it != ranges::end(results); ++r_it)
            {
                ++total_kmers;
                auto result = *r_it;
                int which_hit  = -1;
                for(auto & res : result)
                {
                    if(which_hit == -1)
                    {
                        which_hit = double_group_scaffolds[res.first];
                    }
                    else if(which_hit != double_group_scaffolds[res.first])
                    {
                        which_hit = -2;
                        break;
                    }
                }
                if(which_hit >= 0) 
                {
                    ++unique_kmers[which_hit];
                    // For now, still count unique kmers on an ambiguous read
                    //      because I don't have a justification to ignore it.
                    if(which_group >= 0 && which_group != which_hit)
                    {
                        is_ambiguous = true;
                    }
                    else
                    {
                        which_group = which_hit;
                    }
                }   
            }
            if(is_ambiguous) ++ambiguous_reads;
        }
        return unique_kmers;
    };

    std::vector<std::future<std::vector<std::size_t>>> futures;
    for(std::size_t thread = 0; thread < args.threads; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

    std::vector<std::size_t> unique_totals(group_names.size(),0);
    for(auto &e : futures)
    {
        auto a_unique = e.get();
        for(std::size_t i = 0; i < a_unique.size(); ++i)
        {
            unique_totals[i] += a_unique[i];
        }
    }

    seqan3::debug_stream << unique_totals << "\n";
    seqan3::debug_stream << total_kmers << "\t" << ambiguous_reads << "\n";
    return 1;
}



void speq::scan::_async_count_unique_kmers_per_group(
            speq::args::cmd_arguments & args,
            const std::filesystem::file_time_type & index_file_time,
            const std::vector<std::string> & group_names,
            const std::vector<int> & group_scaffolds,
            const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & index,
            std::vector<std::size_t> & unique_kmers,
            std::vector<std::size_t> & total_kmers)
{
    // Double up the group_scaffolds to match the fm_index values
    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();
    
    //############# Counting Unique Kmers ##############
    // Get the reference sequences
    seqan3::debug_stream << "firsta" << args.in_file_references << "\n";
    seqan3::sequence_file_input fin{args.in_file_references};
    seqan3::debug_stream << "firstb\n";
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
    for(std::size_t thread = 0; thread < args.threads - 1; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

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
        oarchive(index_file_time);
        oarchive(unique_kmers);
        oarchive(total_kmers);
        args.io_file_index.replace_filename(holder);
    }
    seqan3::debug_stream << unique_kmers << "\n";
    seqan3::debug_stream << total_kmers << "\n";
    return;

}

// There is a bug I have yet to find in this
int speq::scan::one(    speq::args::cmd_arguments & args)
{
    // Import the index sequence file (.idx)
    std::vector<std::string> group_names{};
    std::vector<int> group_scaffolds{};
    std::filesystem::path stored_scaffold_file_name;
    std::filesystem::file_time_type scaffold_file_time;
    std::filesystem::file_time_type group_file_time;
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> fm_index;
    {
        args.io_file_index.replace_extension(".idx");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        std::string stored_scaffold_file_name_str;
        iarchive(stored_scaffold_file_name_str);
        stored_scaffold_file_name = std::filesystem::path(stored_scaffold_file_name_str);
        iarchive(scaffold_file_time);
        iarchive(group_file_time);
        iarchive(group_names);
        iarchive(group_scaffolds);
        iarchive(fm_index);
    }

    if(args.in_file_references != stored_scaffold_file_name)
    {
        // Ignore the user input because a different file was used for constructing the index file
        if(std::filesystem::exists(stored_scaffold_file_name))
        {
            args.in_file_references = stored_scaffold_file_name;
        }
        else
        {
            throw std::logic_error{"SPeQ could not find the sequence file used to make the specified index file."};
        }
        
    }

    auto double_group_scaffolds = ranges::view::for_each(group_scaffolds,[](auto c) 
    {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();

    std::filesystem::file_time_type old_index_file_time = std::filesystem::last_write_time(args.io_file_index);
    auto holder = args.io_file_index.filename();
    std::string fn = args.io_file_index.stem();
    args.io_file_index.replace_filename(fn + "_" + std::to_string(args.kmer) + "mer.dat");
    // Import or build the index data file (.dat)
    std::vector<std::size_t> total_kmers_per_group(group_names.size(),0);
    std::vector<std::size_t> unique_kmers_per_group(group_names.size(),0);
    std::filesystem::file_time_type index_file_time;
    if(std::filesystem::exists(args.io_file_index))
    {
        {
            std::ifstream is{args.io_file_index, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index_file_time);
            iarchive(unique_kmers_per_group);
            iarchive(total_kmers_per_group);
        }
        if(old_index_file_time != index_file_time)
        {
            speq::scan::_async_count_unique_kmers_per_group(
                args,
                index_file_time,
                group_names,
                group_scaffolds,
                fm_index,
                unique_kmers_per_group,
                total_kmers_per_group
            );
        }
    }
    else
    {
        speq::scan::_async_count_unique_kmers_per_group(
            args,
            index_file_time,
            group_names,
            group_scaffolds,
            fm_index,
            unique_kmers_per_group,
            total_kmers_per_group
        );
    }
    args.io_file_index.replace_filename(holder);

    // Import the reads
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};
    // seqan3::debug_stream << *fin.begin() << '\n';
    auto chunk_fin = fin | seqan3::views::chunk(args.chunk);
    // auto tcf = *chunk_fin.begin();
    // auto tcf2 = *tcf.begin();
    // seqan3::debug_stream << tcf2 << '\n';
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
                seqan3::search_cfg::output{seqan3::search_cfg::text_position} |
                seqan3::search_cfg::mode{seqan3::search_cfg::all} |
                seqan3::search_cfg::parallel{1};

    std::size_t ambiguous_reads = 0;
    std::size_t total_kmers = 0;
    std::vector<std::size_t> hit_unique_kmers_per_group;
    hit_unique_kmers_per_group.resize(group_names.size(),0);
    std::vector<std::future<std::vector<std::vector<std::size_t>>>> pool_results;
    ThreadPool a_pool(args.threads);
    for(auto it = ranges::begin(chunk_fin); it != ranges::end(chunk_fin); ++it)
    {
        auto a_chunk = *it;
        pool_results.emplace_back(
            a_pool.enqueue([a_chunk, group_names, args, fm_index, double_group_scaffolds]
            {
                seqan3::debug_stream << "go in \t";
                std::size_t a_total_kmers = 0;
                std::size_t a_ambiguous_reads = 0;
                std::vector<std::size_t> a_unique_kmers_per_group(group_names.size(),0);
                
                // Get the sequences from the reads in this chunk
                auto it_seqs = a_chunk | std::views::transform([] (auto s)
                {
                    return seqan3::get<seqan3::field::seq>(s);
                });
                // }) | ranges::to<std::vector>();
                // For each read in this chunk
                int which_group = -1;
                seqan3::debug_stream << "for each seq\n";
                for(auto jt = ranges::begin(it_seqs); jt != ranges::end(it_seqs); jt++)
                {
                    bool is_ambiguous = false;
                    auto a_jt_seq = *jt;
                    // Split the read into kmers
                    auto seq_to_kmer = a_jt_seq  | ranges::views::sliding(args.kmer)
                                            // | ranges::to<std::vector>();
                                            ;
                    seqan3::debug_stream << seq_to_kmer << "\n\n";
                    a_total_kmers += ranges::distance(seq_to_kmer);
                    // Search against the indexed genomes
                    seqan3::debug_stream << "start search\n";
                    auto results = search(seq_to_kmer, fm_index);
                    seqan3::debug_stream << "end search\n" << results << "\n";
                    // For each kmer search
                    for(std::size_t ri = 0; ri < results.size(); ++ri)
                    {
                        int which_hit = -1;
                        for(auto & res : results[ri])
                        {
                            if(which_hit == -1)
                            {
                                which_hit = double_group_scaffolds[res.first];
                            }
                            else if(which_hit != double_group_scaffolds[res.first])
                            {
                                which_hit = -2;
                                break;
                            }
                        }
                        if(which_hit >= 0) 
                        {
                            seqan3::debug_stream << "hit\t";
                            a_unique_kmers_per_group[which_hit]++;
                            // For now, still count unique kmers on an ambiguous read
                            //      because I don't have a justification to ignore it.
                            if(which_group >= 0 && which_group != which_hit)
                            {
                                is_ambiguous = true;
                            }
                            else
                            {
                                which_group = which_hit;
                            }
                        }
                            
                    }
                    if(is_ambiguous) ++a_ambiguous_reads;
                }
                std::vector<std::size_t> total_and_amb;
                total_and_amb.push_back(a_total_kmers);
                total_and_amb.push_back(a_ambiguous_reads);
                std::vector<std::vector<std::size_t>> output;
                output.push_back(a_unique_kmers_per_group);
                output.push_back(total_and_amb);
                return output;
            })
        );
    }
    // For each chunk processed
    for(auto && result : pool_results)
    {
        auto temp_result = result.get();
        auto temp_unique = temp_result[0];
        total_kmers += temp_result[1][0];
        ambiguous_reads += temp_result[1][1];
        for(std::size_t i = 0; i < group_names.size(); ++i)
            hit_unique_kmers_per_group[i] += temp_unique[i];
    }
    seqan3::debug_stream << "Read hits: " << hit_unique_kmers_per_group << "\n";
    seqan3::debug_stream << "Total kmers checked: " << total_kmers << "\n";
    seqan3::debug_stream << "Ambiguous reads: " << ambiguous_reads << "\n";
    return 0;
}

int speq::scan::two(    speq::args::cmd_arguments & args)
{
    seqan3::sequence_file_input fin_1{args.in_file_reads_path_1};
    seqan3::sequence_file_input fin_2{args.in_file_reads_path_2};

    return 0;
}