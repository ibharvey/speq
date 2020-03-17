#include <fm_scanner.h>

int speq::scan::async_one(    speq::args::cmd_arguments & args)
{
    // Import the index sequence file (.idx)
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> fm_index;
    {
        args.io_file_index.replace_extension(".idx");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(fm_index);
    }
    // Import the index data file (.vec)
    std::vector<std::size_t> total_kmers_per_group{};
    std::vector<std::size_t> unique_kmers_per_group{};
    std::vector<std::string> group_names{};
    std::vector<int> double_group_scaffolds{};
    {
        args.io_file_index.replace_extension(".vec");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(args.kmer);
        iarchive(group_names);
        iarchive(double_group_scaffolds);
        iarchive(unique_kmers_per_group);
        iarchive(total_kmers_per_group);

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
            auto seq = seqan3::get<seqan3::field::seq>(record);
            auto kmers = seq | ranges::views::sliding(args.kmer);
            auto results = search(kmers, fm_index, config);
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

// There is a bug I have yet to find in this
int speq::scan::one(    speq::args::cmd_arguments & args)
{
    // Import the index sequence file (.idx)
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> fm_index;
    {
        args.io_file_index.replace_extension(".idx");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(fm_index);
    }
    // Import the index data file (.vec)
    std::vector<std::size_t> total_kmers_per_group{};
    std::vector<std::size_t> unique_kmers_per_group{};
    std::vector<std::string> group_names{};
    std::vector<int> double_group_scaffolds{};
    {
        args.io_file_index.replace_extension(".vec");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(args.kmer);
        iarchive(group_names);
        iarchive(double_group_scaffolds);
        iarchive(unique_kmers_per_group);
        iarchive(total_kmers_per_group);

    }
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
                        for(auto & [idx, pos] : results[ri])
                        {
                            if(which_hit == -1)
                            {
                                which_hit = double_group_scaffolds[idx];
                            }
                            else if(which_hit != double_group_scaffolds[idx])
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