#ifndef FileSystem_hpp
#define FileSystem_hpp

#include <string>




class FileSystem
{
public:
    static std::string getParentDirectory (std::string pathName);
    static void ensureDirectoryExists (std::string pathName);
    static void ensureParentDirectoryExists (std::string dirName);
    static std::string makeFilename (std::string directory, std::string base, std::string extension, int number);
};

#endif
