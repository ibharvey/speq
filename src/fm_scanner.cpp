#include <fm_scanner.h>




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
    std::vector<std::size_t> total_kmers_per_group;
    std::vector<std::size_t> unique_kmers_per_group;
    std::vector<std::string> group_names;
    std::vector<int> double_group_scaffolds;
    {
        args.io_file_index.replace_extension(".vec");
        std::ifstream is{args.io_file_index, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(   group_names,
                    double_group_scaffolds,
                    unique_kmers_per_group,
                    total_kmers_per_group);
    }
    
    // Import the reads
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};
    auto chunk_fin = fin | ranges::views::chunk(args.chunk);
    auto config =   seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}} |
                seqan3::search_cfg::output{seqan3::search_cfg::text_position} |
                seqan3::search_cfg::mode{seqan3::search_cfg::all} |
                seqan3::search_cfg::parallel{1};

    std::size_t ambiguous_reads = 0;
    std::size_t total_reads = 0;
    std::vector<std::size_t> hit_unique_kmers_per_group;
    hit_unique_kmers_per_group.resize(group_names.size(),0);
    size_t count_out = 0;

    std::vector<std::future<std::vector<std::size_t>>> pool_results;
    ThreadPool a_pool(args.threads);
    for(auto it = std::begin(chunk_fin); it != std::end(chunk_fin); it++)
    {
        pool_results.emplace_back(
            a_pool.enqueue([&]
            {
                std::vector<std::size_t> a_unique_kmers_per_group(group_names.size(),0);
                
                // Get the sequences from the reads in this chunk
                auto it_seqs = *it | std::views::transform([] (auto s)
                {
                    return seqan3::get<seqan3::field::seq>(s);
                }) | ranges::to<std::vector>();
                // For each read in this chunk
                for(auto jt = std::begin(it_seqs); jt != std::end(it_seqs); jt++)
                {
                    // Split the read into kmers
                    auto seq_to_kmer = *jt  | ranges::views::sliding(args.kmer)
                                            | ranges::to<std::vector>();
                    // Search against the indexed genomes
                    auto results = search(seq_to_kmer, fm_index, config);
                    // For each kmer search
                    for(std::size_t ri = 0; ri < results.size(); ++ri)
                    {
                        int which_hit = -1;
                        for(auto & [idx, pos] : results[ri])
                        {
                            if(which_hit == -1)
                            {
                                which_hit == double_group_scaffolds[idx];
                            }
                            else if(which_hit != double_group_scaffolds[idx])
                            {
                                which_hit = -2;
                                break;
                            }
                        }
                        if(which_hit >= 0) 
                            a_unique_kmers_per_group[which_hit]++;
                    }
                }

                return a_unique_kmers_per_group;
            })
        );
    }
    // For each chunk processed
    for(auto && result : pool_results)
    {
        auto temp_result = result.get();
        for(std::size_t i = 0; i < group_names.size(); ++i)
            hit_unique_kmers_per_group[i] += temp_result[i];
    }
    seqan3::debug_stream << "Read hits: " << hit_unique_kmers_per_group << "\n";
    return 0;
}

int speq::scan::two(    speq::args::cmd_arguments & args)
{
    seqan3::sequence_file_input fin_1{args.in_file_reads_path_1};
    seqan3::sequence_file_input fin_2{args.in_file_reads_path_2};

    return 0;
}