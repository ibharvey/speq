/*


*/

#include <fm_indexer.h>
#include <arg_parse.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/std/filesystem>

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>
#include <fstream>

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
        // if(args.is_scanner)
        // {
        //     // If there is only one file of reads
        //     if(std::filesystem::path::empty(args.in_file_reads_path_2))
        //     {
        //         speq::scan_1(args);
        //     }
        //     else
        //     {
        //         //speq::scan_2(args);
        //     }
            
        // }
    }
    
}
