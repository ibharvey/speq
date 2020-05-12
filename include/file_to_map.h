#pragma once

#include <string>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>
#include <iostream>

#include <seqan3/std/filesystem>
#include <range/v3/view/split.hpp>

namespace speq
{
    namespace utils
    {
        void string_strip(std::string & s);

        inline bool is_integer(    const std::string & s);
    }
    void file_to_map(   const std::filesystem::path a_path, 
                    std::vector<std::string> & group_names, 
                    std::vector<int> & group_scaffolds,
                    std::vector<int> & group_count);
}


