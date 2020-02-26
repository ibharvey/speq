#include <file_to_map.hpp>

std::map<size_t,std::string> file_to_map(std::filesystem::path a_path)
{
    std::map<size_t,std::string> map_out;
    std::ifstream inf(a_path);
    std::string line,variant;
    size_t genome_index;
    if(inf.is_open())
    {
        while(std::getline(inf,line))
        {
            auto col_ind = line.find(":");
            variant = line.substr(0,col_ind);
            auto gis = line.substr(col_ind + 1);
            std::stringstream ss;
            ss << gis;
            while(ss.rdbuf()->in_avail() > 0)
            {
                ss >> genome_index;
                map_out[genome_index] = variant;
            }
        }
    }
    return map_out;
}