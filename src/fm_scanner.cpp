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
                    if(ranges::min(qesult) > static_cast<int>(args.phred_cutoff))
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
    
    auto percent_each_group = speq::scan::unique_to_percent(
        unique_totals,
        total_kmers,
        unique_kmers_per_group,
        total_kmers_per_group
    );

    seqan3::debug_stream << percent_each_group << "\n";
    seqan3::debug_stream << unique_totals << "\n";
    seqan3::debug_stream << total_kmers << "\t" << ambiguous_reads << "\n";
    std::vector<double> diff_perc(group_names.size(),1.0);
    while(*std::max_element(diff_perc.begin(), diff_perc.end()) > args.precision)
    {
        // Calculate the experimental "total_kmers_per_group"
        std::vector<double> next_tkpg = speq::scan::_async_one_estimate_kmer_per_group
        (
            percent_each_group,
            args,
            group_names,
            double_group_scaffolds,
            fm_index
        );
        // Use the experimental total_kmers_per_group to calculate 
        // the experimental percent_each_group
        auto next_percent_each_group = speq::scan::unique_to_percent(
            unique_totals,
            total_kmers,
            unique_totals,
            next_tkpg
        );

        
        for(std::size_t i = 0; i < group_names.size(); ++i)
        {
            diff_perc[i] = abs(next_percent_each_group[i] - percent_each_group[i]);
        }
        percent_each_group = next_percent_each_group;
        seqan3::debug_stream << percent_each_group << "\n";
        seqan3::debug_stream << next_tkpg << "\n\n";
    }
    return 1;
}

std::vector<double> speq::scan::unique_to_percent(
    std::vector<double> unique_in_reads,
    std::size_t total_in_reads,
    std::vector<std::size_t> unique_in_refs,
    std::vector<std::size_t> total_in_refs
)
{
    std::vector<double> d_unique_in_refs(unique_in_refs.size(),0.0);
    std::copy(unique_in_refs.begin(), 
        unique_in_refs.end(), 
        d_unique_in_refs.begin()
    );

    std::vector<double> d_total_in_refs(total_in_refs.size(),0.0);
    std::copy(total_in_refs.begin(), 
        total_in_refs.end(), 
        d_total_in_refs.begin()
    );
    return speq::scan::unique_to_percent(
        unique_in_reads,
        total_in_reads,
        d_unique_in_refs,
        d_total_in_refs
    );
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
                    if(ranges::min(qesult) > static_cast<int>(args.phred_cutoff))
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
    auto percent_each_group = speq::scan::unique_to_percent(
        unique_totals,
        total_kmers,
        unique_kmers_per_group,
        total_kmers_per_group
    );

    seqan3::debug_stream << percent_each_group << "\n";
    seqan3::debug_stream << unique_totals << "\n";
    seqan3::debug_stream << total_kmers << "\t" << ambiguous_reads << "\n";
    std::vector<double> diff_perc(group_names.size(),1.0);
    while(*std::max_element(diff_perc.begin(), diff_perc.end()) > args.precision)
    {
        // Calculate the experimental "total_kmers_per_group"
        std::vector<double> next_tkpg = speq::scan::_async_one_estimate_kmer_per_group
        (
            percent_each_group,
            args,
            group_names,
            double_group_scaffolds,
            fm_index
        );
        // Use the experimental total_kmers_per_group to calculate 
        // the experimental percent_each_group
        auto next_percent_each_group = speq::scan::unique_to_percent(
            unique_totals,
            total_kmers,
            unique_totals,
            next_tkpg
        );
        
        for(std::size_t i = 0; i < group_names.size(); ++i)
        {
            diff_perc[i] = abs(next_percent_each_group[i] - percent_each_group[i]);
        }
        percent_each_group = next_percent_each_group;
        seqan3::debug_stream << "Percent of each group: " << percent_each_group << "\n";
        seqan3::debug_stream << "Total Kmers per group: " << next_tkpg << "\n\n";
    }
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
    auto second_holder = args.io_file_index.stem();
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
    auto worker = [&, fm_index, group_names, double_group_scaffolds, config, args] ()
    {
        std::vector<std::size_t> unique_kmers(group_names.size(),0);
        auto do_a_count = [&](auto & results, auto & qmers, int & which_group, bool & is_ambiguous)
        {
            auto q_it = ranges::begin(qmers);
            for(auto r_it = ranges::begin(results); r_it != ranges::end(results); ++r_it)
            {
                auto result = *r_it;
                auto qesult = *q_it;
                if(ranges::min(qesult) > static_cast<int>(args.phred_cutoff))
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
        for(auto && [rec1, rec2] : v)
        {
            bool is_ambiguous = false;
            int which_group = -1;
            
            auto f_seq = seqan3::get<seqan3::field::seq>(rec1);
            auto f_qual = seqan3::get<seqan3::field::qual>(rec1);
            auto f_kmers = f_seq | ranges::views::sliding(args.kmer);
            auto f_qmers = f_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            auto f_results = search(f_kmers, fm_index, config);
            do_a_count(f_results, f_qmers, which_group, is_ambiguous);

            auto r_seq = seqan3::get<seqan3::field::seq>(rec2);
            auto r_qual = seqan3::get<seqan3::field::qual>(rec2);
            auto r_kmers = r_seq | ranges::views::sliding(args.kmer);
            auto r_qmers = r_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            auto r_results = search(r_kmers, fm_index, config);
            do_a_count(r_results, r_qmers, which_group, is_ambiguous);

            if(is_ambiguous) ++ambiguous_reads;
        }

        return unique_kmers;
    };

    std::vector<std::future<std::vector<std::size_t>>> futures;
    for(std::size_t thread = 0; thread < args.threads - 1; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

    std::vector<double> unique_totals(group_names.size(),0.0);
    for(auto &e : futures)
    {
        auto a_unique = e.get();
        for(std::size_t i = 0; i < a_unique.size(); ++i)
        {
            unique_totals[i] += static_cast<double>(a_unique[i]) / percent_perfect;
        }
    }
    
    auto percent_each_group = speq::scan::unique_to_percent(
        unique_totals,
        total_kmers,
        unique_kmers_per_group,
        total_kmers_per_group
    );

    seqan3::debug_stream << percent_each_group << "\n";
    seqan3::debug_stream << unique_totals << "\n";
    seqan3::debug_stream << total_kmers << "\t" << ambiguous_reads << "\n";

    std::vector<double> diff_perc(group_names.size(),1.0);
    while(*std::max_element(diff_perc.begin(), diff_perc.end()) > args.precision)
    {
        // Calculate the experimental "total_kmers_per_group"
        std::vector<double> next_tkpg = speq::scan::_async_two_estimate_kmer_per_group
        (
            percent_each_group,
            args,
            group_names,
            double_group_scaffolds,
            fm_index
        );
        // Use the experimental total_kmers_per_group to calculate 
        // the experimental percent_each_group
        auto next_percent_each_group = speq::scan::unique_to_percent(
            unique_totals,
            total_kmers,
            unique_totals,
            next_tkpg
        );

        
        for(std::size_t i = 0; i < group_names.size(); ++i)
        {
            diff_perc[i] = abs(next_percent_each_group[i] - percent_each_group[i]);
        }
        percent_each_group = next_percent_each_group;
        seqan3::debug_stream << percent_each_group << "\n";
        seqan3::debug_stream << next_tkpg << "\n\n";
    }

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
        std::vector<double> unique_kmers(group_names.size(),0);
        
        auto do_a_count = [&](auto & results, auto & qmers, int & which_group, bool & is_ambiguous)
        {
            auto q_it = ranges::begin(qmers);
            for(auto r_it = ranges::begin(results); r_it != ranges::end(results); ++r_it)
            {
                auto result = *r_it;
                auto qesult = *q_it;
                if(ranges::min(qesult) > static_cast<int>(args.phred_cutoff))
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
        for(auto && [rec1, rec2] : v)
        {
            bool is_ambiguous = false;
            int which_group = -1;
            
            auto f_seq = seqan3::get<seqan3::field::seq>(rec1);
            auto f_qual = seqan3::get<seqan3::field::qual>(rec1);
            auto f_kmers = f_seq | ranges::views::sliding(args.kmer);
            auto f_qmers = f_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            auto f_results = search(f_kmers, fm_index, config);
            do_a_count(f_results, f_qmers, which_group, is_ambiguous);

            auto r_seq = seqan3::get<seqan3::field::seq>(rec2);
            auto r_qual = seqan3::get<seqan3::field::qual>(rec2);
            auto r_kmers = r_seq | ranges::views::sliding(args.kmer);
            auto r_qmers = r_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            auto r_results = search(r_kmers, fm_index, config);
            do_a_count(r_results, r_qmers, which_group, is_ambiguous);

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
    auto percent_each_group = speq::scan::unique_to_percent(
        unique_totals,
        total_kmers,
        unique_kmers_per_group,
        total_kmers_per_group
    );

    seqan3::debug_stream << percent_each_group << "\n";
    seqan3::debug_stream << unique_totals << "\n";
    seqan3::debug_stream << total_kmers << "\t" << ambiguous_reads << "\n";

    std::vector<double> diff_perc(group_names.size(),1.0);
    while(*std::max_element(diff_perc.begin(), diff_perc.end()) > args.precision)
    {
        // Calculate the experimental "total_kmers_per_group"
        std::vector<double> next_tkpg = speq::scan::_async_two_estimate_kmer_per_group
        (
            percent_each_group,
            args,
            group_names,
            double_group_scaffolds,
            fm_index
        );
        // Use the experimental total_kmers_per_group to calculate 
        // the experimental percent_each_group
        auto next_percent_each_group = speq::scan::unique_to_percent(
            unique_totals,
            total_kmers,
            unique_totals,
            next_tkpg
        );

        
        for(std::size_t i = 0; i < group_names.size(); ++i)
        {
            diff_perc[i] = abs(next_percent_each_group[i] - percent_each_group[i]);
        }
        percent_each_group = next_percent_each_group;
        seqan3::debug_stream << percent_each_group << "\n";
        seqan3::debug_stream << next_tkpg << "\n\n";
    }
    return 1;
}

std::vector<double> speq::scan::_async_one_estimate_kmer_per_group(
    const std::vector<double> & percent_per_group,
    const speq::args::cmd_arguments & args,
    const std::vector<std::string> & group_names,
    const std::vector<int> & double_group_scaffolds,
    const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & fm_index
)
{
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};

    auto v = fin | seqan3::views::async_input_buffer(args.threads * 2);
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
            seqan3::search_cfg::output{seqan3::search_cfg::text_position};
    auto worker = [&, fm_index, group_names, double_group_scaffolds, config, args]()
    {
        std::vector<double> hits_per_group(group_names.size(),0.0);
        auto do_a_count = [&](auto & results, auto & qmers)
        {
            auto q_it = ranges::begin(qmers);
            for(auto r_it = ranges::begin(results); r_it != ranges::end(results); ++r_it)
            {
                auto result = *r_it;
                auto qesult = *q_it;
                if(ranges::min(qesult) > static_cast<int>(args.phred_cutoff))
                {
                    std::vector<double> hits_per_group_int(group_names.size(), 0.0);
                    for(auto & res : result)
                    {
                        hits_per_group_int[double_group_scaffolds[res.first]] += 1.0;
                    }
                    // Based on the percent_per_group and number of hits for this kmer per group
                    //  we can estimate the representation of this kmer between each group
                    std::vector<double> a_hit_per_group(group_names.size(),0.0);
                    for(std::size_t i = 0; i < group_names.size(); ++i)
                    {
                        a_hit_per_group[i] = hits_per_group_int[i] * percent_per_group[i];
                    }
                    // The total sum of the normalized groups
                    double normalized_sum = std::accumulate(a_hit_per_group.begin(), a_hit_per_group.end(), 0.0);
                    if(normalized_sum > 0.0)
                    {
                        // Divide the norm_hits by norm_sum to get the fraction of this kmer assigned to each group
                        for(std::size_t i = 0; i < group_names.size(); ++i)
                        {
                            a_hit_per_group[i] /= normalized_sum;
                            hits_per_group[i] += a_hit_per_group[i];
                        }
                    }
                }
                ++q_it;
            }
        };
        
        for(auto & rec : v)
        {
            
            auto f_seq = seqan3::get<seqan3::field::seq>(rec);
            auto f_qual = seqan3::get<seqan3::field::qual>(rec);
            auto f_kmers = f_seq | ranges::views::sliding(args.kmer);
            auto f_qmers = f_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            auto f_results = search(f_kmers, fm_index, config);
            do_a_count(f_results, f_qmers);
        }
        return hits_per_group;
    };
    
    std::vector<std::future<std::vector<double>>> futures;
    for(std::size_t thread = 0; thread < args.threads - 1; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

    std::vector<double> hits_per_group_totals(group_names.size(), 0.0);
    for(auto &e : futures)
    {
        auto a_hits_per_group = e.get();
        for(std::size_t i = 0; i < a_hits_per_group.size(); ++i)
        {
            hits_per_group_totals[i] += a_hits_per_group[i];
        }
    }
    return hits_per_group_totals;
}

std::vector<double> speq::scan::_async_two_estimate_kmer_per_group(
    const std::vector<double> & percent_per_group,
    const speq::args::cmd_arguments & args,
    const std::vector<std::string> & group_names,
    const std::vector<int> & double_group_scaffolds,
    const seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & fm_index
)
{
    seqan3::sequence_file_input fin1{args.in_file_reads_path_1};
    seqan3::sequence_file_input fin2{args.in_file_reads_path_2};

    auto combined = seqan3::views::zip(fin1, fin2);
    auto v = combined | seqan3::views::async_input_buffer(args.threads * 2);
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
            seqan3::search_cfg::output{seqan3::search_cfg::text_position};
    auto worker = [&, fm_index, group_names, double_group_scaffolds, config, args]()
    {
        std::vector<double> hits_per_group(group_names.size(),0.0);
        auto do_a_count = [&](auto & results, auto & qmers)
        {
            auto q_it = ranges::begin(qmers);
            for(auto r_it = ranges::begin(results); r_it != ranges::end(results); ++r_it)
            {
                auto result = *r_it;
                auto qesult = *q_it;
                if(ranges::min(qesult) > static_cast<int>(args.phred_cutoff))
                {
                    std::vector<double> hits_per_group_int(group_names.size(), 0.0);
                    for(auto & res : result)
                    {
                        ++hits_per_group_int[double_group_scaffolds[res.first]];
                    }
                    // Based on the percent_per_group and number of hits for this kmer per group
                    //  we can estimate the representation of this kmer between each group
                    std::vector<double> a_hit_per_group(group_names.size(),0.0);
                    for(std::size_t i = 0; i < group_names.size(); ++i)
                    {
                        a_hit_per_group[i] = hits_per_group_int[i] * percent_per_group[i];
                    }
                    // The total sum of the normalized groups
                    double normalized_sum = std::accumulate(a_hit_per_group.begin(), a_hit_per_group.end(), 0.0);
                    // Divide the norm_hits by norm_sum to get the fraction of this kmer assigned to each group
                    for(std::size_t i = 0; i < group_names.size(); ++i)
                    {
                        a_hit_per_group[i] /= normalized_sum;
                        hits_per_group[i] += a_hit_per_group[i];
                    }
                }
                ++q_it;
            }
        };
        
        for(auto && [rec1, rec2] : v)
        {
            
            auto f_seq = seqan3::get<seqan3::field::seq>(rec1);
            auto f_qual = seqan3::get<seqan3::field::qual>(rec1);
            auto f_kmers = f_seq | ranges::views::sliding(args.kmer);
            auto f_qmers = f_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            auto f_results = search(f_kmers, fm_index, config);
            do_a_count(f_results, f_qmers);

            auto r_seq = seqan3::get<seqan3::field::seq>(rec2);
            auto r_qual = seqan3::get<seqan3::field::qual>(rec2);
            auto r_kmers = r_seq | ranges::views::sliding(args.kmer);
            auto r_qmers = r_qual | ranges::views::transform([](auto q) { return q.to_phred();}) | ranges::views::sliding(args.kmer);
            auto r_results = search(r_kmers, fm_index, config);
            do_a_count(r_results, r_qmers);
        }
        return hits_per_group;
    };
    
    std::vector<std::future<std::vector<double>>> futures;
    for(std::size_t thread = 0; thread < args.threads - 1; ++thread)
    {
        futures.push_back(std::async(std::launch::async, worker));
    }

    std::vector<double> hits_per_group_totals(group_names.size(), 0.0);
    for(auto &e : futures)
    {
        auto a_hits_per_group = e.get();
        for(std::size_t i = 0; i < a_hits_per_group.size(); ++i)
        {
            hits_per_group_totals[i] += a_hits_per_group[i];
        }
    }
    return hits_per_group_totals;
}

std::vector<double> speq::scan::unique_to_percent(
    std::vector<double> unique_in_reads,
    std::size_t total_in_reads,
    std::vector<double> unique_in_refs,
    std::vector<double> total_in_refs
)
{
    std::vector<double> output(unique_in_refs.size(),0.0);
    for(std::size_t i = 0; i < unique_in_refs.size(); ++i)
    {
        // Get the percentage of each reference genome that is unique
        double percent_uniques = unique_in_refs[i] / total_in_refs[i];
        output[i] = unique_in_reads[i] / static_cast<double>(total_in_reads) / percent_uniques;
    }
    return output;
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
