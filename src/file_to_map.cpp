#include <file_to_map.h>

void speq::utils::string_strip(std::string & s)
{
    s.erase(std::remove(s.begin(), s.end(), '\t'), s.end());
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
    return;
}

inline bool speq::utils::is_integer(const std::string & s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;

   char * p;
   strtol(s.c_str(), &p, 10);

   return (*p == 0);
}

void speq::file_to_map(   std::filesystem::path a_path,
                    std::vector<std::string> & group_names, 
                    std::vector<int> & group_scaffolds)
{
    std::ifstream inf(a_path);
    std::string line,variant,token;
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
                //while(ss.rdbuf()->in_avail() > 0)
                while(std::getline(ss, token, ','))
                {
                    auto hyphen_loc = token.find('-');
                    std::size_t start_index, end_index;
                    bool is_parsable = true;
                    if(hyphen_loc != std::string::npos)
                    {
                        assert(std::count(token.begin(), token.end(), '-') == 1);
                        // Check that the first number is an integer
                        auto start_s = token.substr(0,hyphen_loc);
                        speq::utils::string_strip(start_s);
                        if(speq::utils::is_integer(start_s))
                        {
                            start_index = std::atoi(start_s.c_str());
                        }
                        else
                        {
                            is_parsable = false;
                        }
                                                // Check that the first number is an integer
                        auto end_s = token.substr(hyphen_loc+1);
                        speq::utils::string_strip(end_s);
                        if(speq::utils::is_integer(end_s))
                        {
                            end_index = std::atoi(end_s.c_str());
                        }
                        else
                        {
                            is_parsable = false;
                        }
                        if(is_parsable)
                        {
                            if (end_index + 1 > group_scaffolds.size()) 
                                group_scaffolds.resize(end_index+1,-1);
                            for(size_t isolate_index = start_index; isolate_index <= end_index; ++isolate_index)
                                group_scaffolds[isolate_index] = variant_index;
                        }
                        else
                        {
                            speq::utils::string_strip(token);
                            std::cerr << "Error in parsing groupings line: " << line
                                      << "\nA non-integer range was detected and ignored at: " << token << std::endl;
                        }
                        

                    }
                    else
                    {
                        speq::utils::string_strip(token);
                        if(speq::utils::is_integer(token))
                        {
                            size_t isolate_index = std::atoi(token.c_str());
                            if(isolate_index + 1 > group_scaffolds.size())
                                group_scaffolds.resize(isolate_index + 1, -1);
                            group_scaffolds[isolate_index] = variant_index;
                        }
                        else
                        {
                            std::cerr << "Error in parsing groupings line: " << line
                                      << "\nA non-integer index was detected and ignored at: " << token << std::endl;
                        }
                        
                    }
                }
            }
        }
    }
    return;
}