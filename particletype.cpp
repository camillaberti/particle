#include "particletype.hpp"
#include <iostream>

ParticleType::ParticleType(std::string name, double mass, int charge): name_{name}, mass_{mass}, charge_{charge} {}

void ParticleType::Print() const {
    std::cout << "Name: " << name_ <<'\n'<< + "Mass: " << mass_ << '\n' << + "Charge: " << charge_ << '\n';
}