#include "arg_parse.h"

speq::args::cmd_arguments speq::args::get_all_arguments(seqan3::argument_parser & parser)
{
    speq::args::cmd_arguments args{};
    args.is_indexer = true;
    args.is_scanner = true;
 
    parser.info.author = "Ian B Harvey";
    parser.info.short_description = "Index reference sequences and scan against reads in one command.";
    parser.info.synopsis = {"speq all [options]"};
    parser.info.version = "0.1.0";

    parser.add_section("Input/Output Arguments");
    parser.add_option(args.in_file_reads_path_1,'1',"forward", "A fast(aq) file of the forward reads.");
    parser.add_option(args.in_file_reads_path_2,'2',"reverse", "A fast(aq) file of the reverse reads.");
    parser.add_option(args.in_file_references,'r',"reference", "The reference sequences in FASTA format.");
    parser.add_option(args.in_file_references_groups,'g',"groups", "Groupings of reference sequences.");
    parser.add_option(args.io_file_index, 'x', "index", "Index file-name prefix.");
    parser.add_option(args.out_file_path, 'o', "output", "Output file of percentages.");
    parser.add_flag(args.is_force, 'f', "force", "Force the program to overwrite an existing file.");

    parser.add_section("Misc");
    parser.add_option(args.kmer, 'k', "kmer", "Size of the kmer used in searching the references.");
    parser.add_option(args.phred_cutoff, '\0', "phred-cutoff", "Single nucleotide phred score cutoff value for a kmer to be included in analysis.");
    parser.add_option(args.fixed_accuracy, '\0', "fixed-accuracy", "Force program to use a specific accuracy per basepair. Default uses Phred.", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{0.0,1.0});
    //parser.add_option(args.chunk, 'c', "chunk", "Chunk size when pulling the input file(s).");
    parser.add_option(args.threads, 't', "threads", "Number of threads to use.", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{2,static_cast<double>(std::thread::hardware_concurrency())});

    try
    {
        parser.parse();
        // Check the input file #1
        speq::args::check_in_file(args.in_file_reads_path_1);
        speq::args::check_in_file(args.in_file_reads_path_2);
        speq::args::check_out_file(args.out_file_path, args.is_force);
        args.is_parsed = true;
        return args;
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "speq all | argument parsing error:  " << ext.what() << "\n";
        return args;
    }
}

speq::args::cmd_arguments speq::args::get_scan_arguments(seqan3::argument_parser & parser)
{
    speq::args::cmd_arguments args{};
    args.is_scanner = true;
    
    parser.info.author = "Ian B Harvey";
    parser.info.short_description = "Scan reads against a pre-built reference index.";
    parser.info.synopsis = {"speq scan [options]"};
    parser.info.version = "0.1.0";

    parser.add_section("Input/Output Arguments");
    parser.add_option(args.in_file_reads_path_1,'1',"forward", "A fast(aq) file of the forward reads.");
    parser.add_option(args.in_file_reads_path_2,'2',"reverse", "A fast(aq) file of the reverse reads.");
    parser.add_option(args.io_file_index, 'x', "index", "Input index file-name prefix.");
    parser.add_option(args.out_file_path, 'o', "output", "Output file of percentages.");
    parser.add_flag(args.is_force, 'f', "force", "Force the program to overwrite an existing file.");

    parser.add_section("Misc");
    parser.add_option(args.kmer, 'k', "kmer", "Size of the kmer used in searching the references.");
    parser.add_option(args.phred_cutoff, '\0', "phred-cutoff", "Single nucleotide phred score cutoff value for a kmer to be included in analysis.");
    parser.add_option(args.fixed_accuracy, '\0', "fixed-accuracy", "Force program to use a specific accuracy per basepair", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{0.0,1.0});
    //parser.add_option(args.chunk, 'c', "chunk", "Chunk size when pulling the input file(s).");
    parser.add_option(args.threads, 't', "threads", "Number of threads to use.", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{2,static_cast<double>(std::thread::hardware_concurrency())});

    try
    {
        parser.parse();
        // Check the input file #1
        speq::args::check_in_file(args.in_file_reads_path_1);
        speq::args::check_in_file(args.in_file_reads_path_2);
        speq::args::check_in_file(args.io_file_index);
        speq::args::check_out_file(args.out_file_path, args.is_force);
        args.is_parsed = true;
        return args;
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "speq scan | argument parsing error:  " << ext.what() << "\n";
        return args;
    }
}

speq::args::cmd_arguments speq::args::get_index_arguments(seqan3::argument_parser & parser)
{
    speq::args::cmd_arguments args{};
    args.is_indexer = true;

    parser.info.author = "Ian B Harvey";
    parser.info.short_description = "Index reference genomes/sequences.";
    parser.info.synopsis = {"speq index [options]"};
    parser.info.version = "0.1.0";

    parser.add_section("Input/Output Arguments");
    parser.add_option(args.in_file_references,'r',"reference", "The reference sequences in FASTA format.");
    parser.add_option(args.in_file_references_groups,'g',"groups", "Groupings of reference sequences.");
    parser.add_option(args.io_file_index, 'x', "index", "Output index file-name prefix.");
    parser.add_flag(args.is_force, 'f', "force", "Force the program to overwrite an existing file.");

    parser.add_section("Misc");
    //parser.add_option(args.chunk, 'c', "chunk", "Chunk size when pulling the input file(s).");
    parser.add_option(args.threads, 't', "threads", "Number of threads to use.", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{2,static_cast<double>(std::thread::hardware_concurrency())});

    try
    {
        parser.parse();
        // Check the input files
        speq::args::check_in_file(args.in_file_references);
        speq::args::check_in_file(args.in_file_references_groups);
        speq::args::check_out_file(args.out_file_path, args.is_force);
        args.is_parsed = true;
        return args;
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "speq index | argument parsing error:  " << ext.what() << "\n";
        return args;
    }
}

speq::args::cmd_arguments speq::args::initialize_argument_parser( int argc, char ** argv)
{
    seqan3::argument_parser top_level_parser{"speq", argc, argv, true, {"index", "scan", "all"}};
    
    top_level_parser.info.app_name = "speq";
    top_level_parser.info.author = "Ian B Harvey";
    top_level_parser.info.short_description = "Species PErcent Quantifier";
    top_level_parser.info.synopsis = {"speq [index|scan|all] [options]"};
    top_level_parser.info.version = "0.1.0";
    try
    {
        top_level_parser.parse();
    }
    catch(seqan3::argument_parser_error const & err)
    {
        seqan3::debug_stream << "SPeQ Error: " << err.what() << "\n";
        speq::args::cmd_arguments args;
        return args;
    }
    
    seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser();
    if(sub_parser.info.app_name == std::string_view{"speq-all"})
    {
        return speq::args::get_all_arguments(sub_parser);
    }
    else if(sub_parser.info.app_name == std::string_view{"speq-index"})
    {
        return speq::args::get_index_arguments(sub_parser);
    }
    else if(sub_parser.info.app_name == std::string_view{"speq-scan"})
    {
        return speq::args::get_scan_arguments(sub_parser);
    }
    else
    {
        throw std::logic_error{"SPeQ does not contain the subparser " + sub_parser.info.app_name};
    }
}

void speq::args::check_in_file(std::filesystem::path & a_path)
{
    if(!a_path.empty())
    {
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
    
}

void speq::args::check_out_file(std::filesystem::path & a_path, const bool is_force)
{
    if(a_path.is_relative())
    {
        a_path = std::filesystem::current_path() / a_path;
        if(std::filesystem::exists(a_path) && !is_force)
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "Cowardly refusing to use an existing output file. Use '-f' to overwrite."));
        }
    }

    return;
}