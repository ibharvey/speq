#include <khash_indexer.h>


void speq::kmer_hash::get_unique_kmers(  const cmd_arguments args,
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
    
    // Make a cereal-stream for each genome to hold all_kmers
    std::vector<cereal::BinaryOutputArchive> all_kmer_hash_ref;
    std::vector<cereal::BinaryOutputArchive> unique_kmer_hash_ref;
    all_kmer_hash_ref.reserve(group_names.size());
    unique_kmer_hash_ref.reserve(group_names.size());
    for(std::size_t i = 0; i < group_names.size(); ++i)
    {
        // All
        std::ofstream os("all_hashes_" + group_names[i] + ".kmr", std::ios::binary);
        cereal::BinaryOutputArchive oar(os);
        all_kmer_hash_ref.push_back(oar);
        // Unique
        std::ofstream os2("unique_hashes_" + group_names[i] + ".kmr", std::ios::binary);
        cereal::BinaryOutputArchive oar2(os2);
        unique_kmer_hash_ref.push_back(oar2);
    }

    seqan3::debug_stream << "pull genomes\n";

    for(auto & rec : fgenomes)
    {
        // If this scaffold was indicated for use in the groupings file...
        if(group_scaffolds[isolate_index] != -1)
        {
            seqan3::debug_stream << seqan3::get<seqan3::field::id>(rec) << "\n";
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
            // all_hash_reference[gi_index].push_back(a_hash);
            // all_hash_reference[gi_index].push_back(b_hash);
            all_kmer_hash_ref[gi_index](a_hash);
            all_kmer_hash_ref[gi_index](b_hash);
            auto a_unique_hash = ranges::actions::unique(a_hash) | ranges::to<std::vector>();
            auto b_unique_hash = ranges::actions::unique(b_hash) | ranges::to<std::vector>();
            // unique_hash_reference[gi_index].push_back(a_unique_hash);
            // unique_hash_reference[gi_index].push_back(b_unique_hash);
            unique_kmer_hash_ref[gi_index](a_unique_hash);
            unique_kmer_hash_ref[gi_index](b_unique_hash);
        }
        ++isolate_index;
    }
    
    seqan3::debug_stream << "pulled genomes\n";

    // For each strain-group
    std::vector<std::vector<std::size_t>> total_union_kmers;
    // for(auto it = std::begin(all_hash_reference); it != std::end(all_hash_reference); ++it)
    for(auto it = all_kmer_hash_ref.begin(); it != all_kmer_hash_ref.end(); ++it)
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
    
    seqan3::debug_stream << "union of kmers\n";

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

    unique_kmers.resize(single_unique_kmers.size());
    auto all_hash_it = std::begin(all_hash_reference);
    for(std::size_t i = 0; i < single_unique_kmers.size(); ++i)
    {
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
            
            for(std::size_t k = 0; k < max_repeats; ++k)
            {
                unique_kmers[i].push_back(*jt);
            }
        }
        ++all_hash_it;
    }

    return;
}

int speq::kmer_hash::index(const cmd_arguments args)
{
    seqan3::debug_stream << "Starting Run:\n";
    // Organize the reference sequences by strain/variant type
    std::vector<int> group_scaffolds;
    std::vector<std::string> group_names;
    speq::file_to_map(args.in_file_references_groups, group_names, group_scaffolds);
    seqan3::debug_stream << "Mapped file.:" << group_names << " | " << group_scaffolds <<  "\n";

    // Find unique kmers in each scaffold/genome-set
    std::vector<std::vector<std::size_t>> unique_kmers;
    std::vector<std::size_t> total_kmers_count;
    speq::kmer_hash::get_unique_kmers(args, group_names, group_scaffolds, unique_kmers, total_kmers_count);
    seqan3::debug_stream << "Got Unique Kmers.\n";
    return 0;
}
    
