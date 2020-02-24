/*


*/

#include <aligner.hpp>
#include <arg_parse.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/exceptions.hpp>

int main(int argc, char ** argv)
{
    try
    {
        auto args = initialize_argument_parser("SPeQ", argc, argv);
        if(args.in_file_reads_path_2)
        {
            return speq_run_2(args);
        }
        else
        {
            return speq_run_1(args);
        }
    }
    catch (seqan3::validation_error const & ext)
    {
        seqan3::debug_stream << "Invalid argument: " << ext.what() << "\n"; // customise your error message
        return -1;
    }
}
