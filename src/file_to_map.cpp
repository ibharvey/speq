#include <file_to_map.hpp>

void file_to_map(   std::filesystem::path a_path,
                    std::vector<std::string> & group_names, 
                    std::vector<int> & group_scaffolds)
{
    std::ifstream inf(a_path);
    std::string line,variant;
    size_t isolate_index; size_t variant_index = -1;
    if(inf.is_open())
    {
        while(std::getline(inf,line))
        {
            // Ignore in-line comments in the file
            if(line.find("#") != std::string::npos) line = line.substr(0,line.find("#"));
            // Ignore full-line comments in the file
            if(line.find(":") != std::string::npos)
            {
                // Get the strain/species group name
                auto col_ind = line.find(":");
                variant = line.substr(0,col_ind);
                group_names.push_back(variant);
                variant_index++;
                // Get the scaffold indices that are assigned to this strain/species
                auto gis = line.substr(col_ind + 1);
                std::stringstream ss;
                ss << gis;
                while(ss.rdbuf()->in_avail() > 0)
                {
                    ss >> isolate_index;
                    if(isolate_index + 1 > group_scaffolds.size())
                        group_scaffolds.resize(isolate_index+1,-1);
                    group_scaffolds[isolate_index] = variant_index;
                }
            }
        }
    }
    return;
}