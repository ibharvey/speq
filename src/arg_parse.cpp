#include "arg_parse.hpp"

cmd_arguments initialize_argument_parser(const std::string name, int argc, char ** argv)
{
    cmd_arguments args{};
    seqan3::argument_parser parser{name, argc, argv};
    
    parser.info.app_name = name;
    parser.info.author = "Ian B Harvey";
    parser.info.short_description = "Species PErcent Quantifier";
    parser.info.synopsis = {"./speq [options]"};
    parser.info.version = "0.1";

    parser.add_section("Input/Output Arguments");
    parser.add_option(args.in_file_reads_path_1,'1',"foward", "A fast(aq) file of the forward reads.");
    parser.add_option(args.in_file_reads_path_2,'2',"reverse", "A fast(aq) file of the reverse reads.");
    parser.add_option(args.in_file_references,'r',"reference", "The reference sequences in FASTA format.");
    parser.add_option(args.in_file_references_groups,'g',"groups", "Groupings of reference sequences.");
    parser.add_option(args.out_file_path, 'o', "output", "Output file of percentages.");
    parser.add_flag(args.force, 'f', "force", "Force the program to overwrite an existing file.");

    parser.add_section("Misc");
    parser.add_option(args.chunk, 'c', "chunk", "Chunk size when pulling the input file(s).");
    parser.add_option(args.kmer, 'k', "kmer", "Size of the kmer used in searching the references.");
    parser.add_option(args.threads, 't', "threads", "Number of threads to use.", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1,static_cast<double>(std::thread::hardware_concurrency())});

    parser.parse();
    check_arguments(args);
    
    return args;
}

void check_in_file(std::filesystem::path & a_path)
{
    // Check the input file #1
    if(a_path.is_relative())
    {
        if(std::filesystem::exists(std::filesystem::current_path() / a_path))
        {
            a_path = std::filesystem::current_path() / a_path;
        }
        else
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "The relative-path file ", a_path, " was not found."));
        }
    }else{
        if(!std::filesystem::exists(a_path))
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "The full-path file ", a_path, " was not found."));
        }
    }
}

void check_arguments(cmd_arguments & args)
{
    // Check the input file #1
    check_in_file(args.in_file_reads_path_1);
    check_in_file(args.in_file_reads_path_2);
    check_in_file(args.in_file_references);
    check_in_file(args.in_file_references_groups);

    if(args.out_file_path.is_relative())
    {
        args.out_file_path = std::filesystem::current_path() / args.out_file_path;
        if(std::filesystem::exists(args.out_file_path) && !args.force)
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "Cowardly refusing to use an existing output file. Use '-f' to overwrite."));
        }
    }

    return;
}