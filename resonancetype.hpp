#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP
#include "particletype.hpp"

class ResonanceType : public ParticleType{
    public:
    ResonanceType(std::string name, double mass, int charge, double width);
    double GetWidth() const {return width_;}
    void Print() const;
    private: 
    double const width_;
};


#endif