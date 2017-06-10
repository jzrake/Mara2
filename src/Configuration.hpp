#ifndef Configuration_hpp
#define Configuration_hpp

#include "Mara.hpp"




class Configuration
{
public:
    Configuration (int argc, const char* argv[]);
    ~Configuration();
    SimulationSetup fromCheckpoint (std::string filename);
    SimulationSetup fromLuaFile (std::string filename);
    int launchFromScript (MaraSession& mara, std::string filename);

private:
    class LuaState;
    std::shared_ptr<LuaState> luaState;
    std::string luaCommandLine;
};


#endif
