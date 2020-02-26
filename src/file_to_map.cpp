#include <file_to_map.hpp>

std::map<size_t,std::string> file_to_map(std::filesystem::path a_path)
{
    std::map<size_t,std::string> map_out;
    std::ifstream inf(a_path);
    std::string line,variant;
    size_t genome_index;
    std::stringstream ss;
    if(inf.is_open())
    {
        while(std::getline(inf,line))
        {
            auto temp = ranges::views::split(line,':');
            auto temp_it = std::begin(temp);
            variant = *temp_it | ranges::to<std::string>;
            temp_it++;
            auto inds = ranges::views::split(*temp_it,',');
            for(auto it = std::begin(inds); it != std::end(inds); it++)
            {
                ss << *it;
                ss >> genome_index;
                map_out[genome_index] = variant;
            }
        }
    }
    return map_out;
}