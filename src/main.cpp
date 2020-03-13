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
    {
        // Make a cereal archive
        std::vector<std::size_t> faks{10,20,30};
        std::ofstream os("output.i", std::ios::binary);
        cereal::BinaryOutputArchive oar(os);
        oar(faks);
        std::vector<std::size_t> faks2{15,22,39};
        oar(faks2);
    }
    {
        std::ifstream is1("output.i", std::ios::binary);
        cereal::BinaryInputArchive iar1(is1);
        for(size_t i = 0; i < 2; ++i)
        {
            std::vector<std::size_t> in1;
            iar1(in1);
            std::ifstream is2("output.i", std::ios::binary);
            cereal::BinaryInputArchive iar2(is2);
            for(size_t j = 0; j < 2; ++j)
            {
                std::vector<std::size_t> in2;
                iar2(in2);
                for(auto it = in1.begin(); it != in1.end(); ++it)
                {
                    for(auto jt = in2.begin(); jt != in2.end(); ++jt)
                    {
                        std::cout << *it + *jt << ", ";
                    }
                    std::cout << "\t";
                }
                std::cout << std::endl;
            }
        }
    }
    
    
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
