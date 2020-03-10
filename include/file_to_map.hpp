
#include <string>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>

#include <seqan3/std/filesystem>
#include <range/v3/view/split.hpp>

void file_to_map(   std::filesystem::path a_path, 
                    std::vector<std::string> & group_names, 
                    std::vector<std::size_t> & group_scaffolds);