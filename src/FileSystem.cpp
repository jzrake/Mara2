#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include "FileSystem.hpp"




std::vector<std::string> FileSystem::splitPath (std::string pathName)
{
    auto remaining = pathName;
    auto dirs = std::vector<std::string>();

    while (true)
    {
        auto slash = remaining.find ('/');

        if (slash == std::string::npos)
        {
            dirs.push_back (remaining);
            break;
        }
        dirs.push_back (remaining.substr (0, slash));
        remaining = remaining.substr (slash + 1);
    }
    return dirs;
}

std::string FileSystem::fileExtension (std::string pathName)
{
    auto dot = pathName.rfind ('.');

    if (dot != std::string::npos)
    {
        return pathName.substr (dot);
    }
    return "";
}

std::string FileSystem::getParentDirectory (std::string pathName)
{
	std::string::size_type lastSlash = pathName.find_last_of ("/");
    return pathName.substr (0, lastSlash);
}

void FileSystem::ensureDirectoryExists (std::string dirName)
{
    auto path = std::string (".");

    for (auto dir : splitPath (dirName))
    {
        path += "/" + dir;
        mkdir (path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
}

void FileSystem::ensureParentDirectoryExists (std::string pathName)
{
	std::string parentDir = getParentDirectory (pathName);
    mkdir (parentDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

std::string FileSystem::makeFilename (
    std::string directory,
    std::string base,
    std::string extension,
    int number,
    int rank)
{
    std::stringstream filenameStream;
    filenameStream << directory << "/" << base;

    if (number >= 0)
    {
        filenameStream << "." << std::setfill ('0') << std::setw (4) << number;
    }
    if (rank != -1)
    {
        filenameStream << "." << std::setfill ('0') << std::setw (4) << rank;
    }

    filenameStream << extension;
    return filenameStream.str();
}
