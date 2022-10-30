#include "resonancetype.hpp"
#include <iostream>
ResonanceType::ResonanceType(std::string name, double mass, int charge, double width): ParticleType(name,mass,charge), width_(width) {}

void ResonanceType::Print() const {
    ParticleType::Print();
    std::cout<<"Width of resonance: " << width_ << '\n';
}
