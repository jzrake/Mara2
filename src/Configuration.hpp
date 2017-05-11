#ifndef Configuration_hpp
#define Configuration_hpp

#include "Mara.hpp"




class Configuration
{
public:
    class LuaState;

    Configuration();
    ~Configuration();
    SimulationSetup fromLuaFile (std::string filename);

    std::shared_ptr<LuaState> luaState;
};


#endif
