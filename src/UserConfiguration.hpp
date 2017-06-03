#ifndef UserConfiguration_hpp
#define UserConfiguration_hpp

#include "Variant.hpp"




class UserConfiguration
{
public:
    UserConfiguration (int argc=0, const char* argv[]=nullptr);
    const Variant::NamedValues& getMembers() const;
    void describe() const;
private:
    Variant::NamedValues members;
};

#endif
