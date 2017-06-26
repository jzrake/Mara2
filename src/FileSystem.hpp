#ifndef FileSystem_hpp
#define FileSystem_hpp

#include <string>
#include <vector>



class FileSystem
{
public:
    static std::vector<std::string> splitPath (std::string pathName);
    static std::string fileExtension (std::string pathName);
    static std::string getParentDirectory (std::string pathName);
    static void ensureDirectoryExists (std::string pathName);
    static void ensureParentDirectoryExists (std::string dirName);
    static std::string makeFilename (
        std::string directory,
        std::string base,
        std::string extension,
        int number,
        int rank=-1);
};

#endif
