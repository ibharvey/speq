
#include <map>
#include <string>
#include <fstream>
#include <iterator>
#include <sstream>

#include <seqan3/std/filesystem>
#include <range/v3/view/split.hpp>

std::map<size_t,std::string> file_to_map(std::filesystem::path a_path);