#include <fm_scanner.h>




int speq::scan::one(    const speq::args::cmd_arguments args)
{
    seqan3::sequence_file_input fin{args.in_file_reads_path_1};
    auto chunk_fin = fin | ranges::views::chunk(args.chunk);

    std::size_t ambiguous_reads = 0;
    std::size_t total_reads = 0;
    std::vector<std::size_t> hit_unique_kmers_per_group;


}

int speq::scan::two(    const speq::args::cmd_arguments args)
{
    seqan3::sequence_file_input fin_1{args.in_file_reads_path_1};
    seqan3::sequence_file_input fin_2{args.in_file_reads_path_2};
}