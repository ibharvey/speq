/*


*/

#include <fm_indexer.h>
#include <fm_scanner.h>
#include <arg_parse.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/exceptions.hpp>

#include <vector>

int main(int argc, char ** argv)
{
    auto args = speq::args::initialize_argument_parser(argc, argv);
    if(args.is_parsed)
    {
        
        if(args.is_indexer)
        {
            speq::fm::index(args);
        }
        if(args.is_scanner)
        {
            // If there is only one file of reads
            if(args.in_file_reads_path_2.empty())
            {
                speq::scan::async_one(args);
            }
            else
            {
                seqan3::debug_stream << "Path2: " << args.in_file_reads_path_2 << "\n";
                //speq::scan::two(args);
            }
            
        }
    }
    
}
