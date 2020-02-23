#include "arg_parse.hpp"

cmd_arguments initialize_argument_parser(const std::string name, int argc, char ** argv)
{
    cmd_arguments args{};
    seqan3::argument_parser parser{"HOPA", argc, argv};
    
    parser.info.app_name = "hopa";
    parser.info.author = "Ian B Harvey";
    parser.info.short_description = "Homogenize Orientation by Pairwise Alignment.";
    parser.info.synopsis = {"./hopa -t 4 ~/Documents/my_sequences.fa"};
    parser.info.version = "0.1";

    parser.add_positional_option(args.in_file_path, "A fast(aq) file to be oriented.");
    parser.add_option(args.out_file_path, 'o', "output", "Output file for oriented sequences.");
    parser.add_flag(args.reverse, 'r', "reverse", "Reverse the global orientation of the output "
                                                    "opposite the first sequence in the file.");
    // TODO: Implement gzip and bzip functionality
    /*
    parser.add_flag(args.use_bzip,'j',"bzip", "BZip the output file.");
    parser.add_flag(args.use_gzip,'z',"gzip", "GZip the output file.");
    */
    parser.add_option(args.threads, 't', "threads", "Number of threads to use.", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1,std::thread::hardware_concurrency()});

    int error_handle = 0;
    if(check_arguments(args, error_handle) < 0) return args;
    try
    {
        parser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext)
    {
        seqan3::debug_stream << "Invalid argument: " << ext.what() << "\n"; // customise your error message
    }
    return args;
}


int check_arguments(cmd_arguments & args, int & err_handle)
{
    // Check the input file
    if(args.in_file_path.is_relative())
    {
        if(std::filesystem::exists(std::filesystem::current_path() / args.in_file_path))
        {
            args.in_file_path = std::filesystem::current_path() / args.in_file_path;
        }
        else
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "The relative-path file ", args.in_file_path, " was not found."));
        }
    }else{
        if(!std::filesystem::exists(args.in_file_path))
        {
            throw seqan3::validation_error(seqan3::detail::to_string( "The full-path file ", args.in_file_path, " was not found."));
        }
    }
    // Set the output file path if not already a good value.
    if(args.out_file_path.is_relative())
    {
        args.out_file_path = std::filesystem::current_path() / args.out_file_path;
    }
    return err_handle;
}