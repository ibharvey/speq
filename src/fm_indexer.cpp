#include <fm_indexer.h>


// This function needs to be re-written to avoid 
// loading all genomes into containers at once. 
// Would dramatically help both CPU and memory overhead.

void speq::fm::generate_fm_index(
    speq::args::cmd_arguments & args,
    const std::filesystem::file_time_type & scaffold_file_time,
    const std::filesystem::file_time_type & group_file_time,
    const std::vector<std::string> & group_names,
    const std::vector<int> & group_scaffolds
)
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
        oarchive(args.in_file_references.string());
        oarchive(scaffold_file_time);
        oarchive(group_file_time);
        oarchive(group_names);
        oarchive(group_scaffolds);
        oarchive(index);
    }
    return;
}


int speq::fm::index( speq::args::cmd_arguments & args)
{
    // If the index_file already exists and the user isn't forcing a re-write
    std::filesystem::path old_scaffold_file_name;
    std::filesystem::file_time_type old_group_file_time;
    std::filesystem::file_time_type old_scaffold_file_time;
    std::filesystem::file_time_type scaffold_file_time = std::filesystem::last_write_time(args.in_file_references);
    std::filesystem::file_time_type group_file_time = std::filesystem::last_write_time(args.in_file_references_groups);
    // Organize the reference sequences by strain/variant type
    std::vector<int> group_scaffolds;
    std::vector<std::string> group_names;
    speq::file_to_map(args.in_file_references_groups, group_names, group_scaffolds);
    if(std::filesystem::exists(args.io_file_index) && !args.is_force)
    {
        // Check if the index, or files that construct the index, were changed since last indexing
        {
            args.io_file_index.replace_extension(".idx");
            std::ifstream is{args.io_file_index, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            std::string old_scaffold_file_name_str;
            iarchive(old_scaffold_file_name_str);
            old_scaffold_file_name = std::filesystem::path(old_scaffold_file_name_str);
            iarchive(old_scaffold_file_time);
            iarchive(old_group_file_time);
        }
        
        if(
            old_scaffold_file_name != args.in_file_references || 
            old_group_file_time != group_file_time || 
            old_scaffold_file_time != scaffold_file_time
        )
        {
            speq::fm::generate_fm_index(
                args, 
                scaffold_file_time,
                group_file_time, 
                group_names,
                group_scaffolds
            );
        }
    }
    else
    {
        // Generate the fm index from the sequences
        speq::fm::generate_fm_index(
            args, 
            scaffold_file_time,
            group_file_time, 
            group_names,
            group_scaffolds
        );
    }
    return 0;
}