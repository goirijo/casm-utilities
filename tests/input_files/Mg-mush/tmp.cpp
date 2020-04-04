#include <filesystem>
#include <iostream>
#include <string>

std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

int main()
{
    namespace fs = std::filesystem;

    for (const auto& d : fs::recursive_directory_iterator("."))
    {
        if (!d.is_directory())
        {
            continue;
        }

        fs::path dir(d);
        if (dir.filename() != "cleave_0.000000")
        {
            continue;
        }

        int x,y;
        auto cleave_dir=dir;

        std::string shift_dirname=cleave_dir.parent_path().filename();
        std::vector<std::string> split_by_period=split(shift_dirname,'.');
        x=std::stoi(split_by_period.back());

        std::vector<std::string> split_by_underscore=split(split_by_period[0],'_');
        y=std::stoi(split_by_underscore.back());

        std::cout<<cleave_dir<<std::endl;
        std::cout<<x<<"    "<<y<<std::endl;

    }
    return 0;
}
