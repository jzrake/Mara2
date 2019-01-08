#pragma once
#include "../Mara.hpp"




class BinaryTorque : public SubProgram
{
public:
    int run (int argc, const char* argv[]) override;
};



class BinaryTorqueStressCalculation : public SubProgram
{
public:
    int run (int argc, const char* argv[]) override;
private:
    void processFile(std::string fname) const;
};
