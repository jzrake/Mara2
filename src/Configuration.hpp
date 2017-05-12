#ifndef Configuration_hpp
#define Configuration_hpp

#include "Mara.hpp"




class Configuration
{
public:
    Configuration();
    ~Configuration();
    SimulationSetup fromLuaFile (std::string filename);

private:
    class LuaState;
    std::shared_ptr<LuaState> luaState;
};


#endif
