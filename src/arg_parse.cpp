#include "arg_parse.hpp"

cmd_arguments initialize_argument_parser(const std::string name, int argc, char ** argv)
{
    cmd_arguments args{};
    seqan3::argument_parser parser{name, argc, argv};
    
    parser.info.app_name = name;
    parser.info.author = "Ian B Harvey";
    parser.info.short_description = "Homogenize Orientation by Pairwise Alignment.";
    parser.info.synopsis = {"./speq [options]"};
    parser.info.version = "0.1";

    parser.add_section("Input/Output Arguments");
    parser.add_option(args.in_file_reads_path_1,'1',"foward", "A fast(aq) file of the forward reads.");
    parser.add_option(args.in_file_reads_path_2,'2',"reverse", "A fast(aq) file of the reverse reads.");
    parser.add_option(args.in_file_references,'r',"reference", "The reference sequences in FASTA format.");
    parser.add_option(args.in_file_references_groups,'g',"groups", "Groupings of reference sequences.");
    parser.add_option(args.out_file_path, 'o', "output", "Output file of percentages.");

    parser.add_section("Misc");
    parser.add_option(args.chunk, 'c', "chunk", "Chunk size when pulling the input file(s).");
    parser.add_option(args.kmer, 'k', "kmer", "Size of the kmer used in searching the references.");
    parser.add_flag(args.force, 'f', "force", "Force the program to overwrite an existing file.");
    // TODO: Implement gzip and bzip functionality
    /*
    parser.add_flag(args.use_bzip,'j',"bzip", "BZip the output file.");
    parser.add_flag(args.use_gzip,'z',"gzip", "GZip the output file.");
    */
    parser.add_option(args.threads, 't', "threads", "Number of threads to use.", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1,static_cast<double>(std::thread::hardware_concurrency())});

    parser.parse();
    check_arguments(args);
    
    return args;
}


void check_arguments(cmd_arguments & args)
{
    // Check the input file
    if(args.in_file_reads_path_1.is_relative())
    {
        if(std::filesystem::exists(std::filesystem::current_path() / args.in_file_reads_path_1))
        {
            args.in_file_reads_path_1 = std::filesystem::current_path() / args.in_file_reads_path_1;
        }
        else
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "The relative-path file ", args.in_file_reads_path_1, " was not found."));
        }
    }else{
        if(!std::filesystem::exists(args.in_file_reads_path_1))
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "The full-path file ", args.in_file_reads_path_1, " was not found."));
        }
    }

    // Cannot both be overwriting and setting an output file path
    if(args.overwrite && args.out_file_path != "output.fa" && args.out_file_path != args.in_file_reads_path_1)
    {
        throw seqan3::validation_error(seqan3::detail::to_string( "Either designate an output file OR overwrite the input file."));
    }

    // Set the output file path if not already a good value.
    if(args.overwrite)
    {
        args.out_file_path = args.in_file_reads_path_1;
    }
    else if(args.out_file_path.is_relative())
    {
        args.out_file_path = std::filesystem::current_path() / args.out_file_path;
        if(std::filesystem::exists(args.out_file_path) && !args.force)
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "Cowardly refusing to use an existing output file. Use 'f' to overwrite."));
        }
    }

    return;
}