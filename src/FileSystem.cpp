#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include "FileSystem.hpp"


std::string FileSystem::getParentDirectory (std::string pathName)
{
	std::string::size_type lastSlash = pathName.find_last_of ("/");
    return pathName.substr (0, lastSlash);
}

void FileSystem::ensureDirectoryExists (std::string dirName)
{
    mkdir (dirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

void FileSystem::ensureParentDirectoryExists (std::string pathName)
{
	std::string parentDir = getParentDirectory (pathName);
    mkdir (parentDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

std::string FileSystem::makeFilename (std::string directory, std::string base, std::string extension, int number)
{
    std::ostringstream filenameStream;
    filenameStream << directory << "/" << base;

    if (number >= 0)
    {
        filenameStream << "." << std::setfill ('0') << std::setw (6) << number;
    }
    filenameStream << extension;
    return filenameStream.str();
}
