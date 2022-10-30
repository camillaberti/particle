#include "particletype.hpp"
#include "resonancetype.hpp"
#include "particle.hpp"
#include <array>
int main(){
    /*ParticleType proton ("Protone",2,1);
    ResonanceType k ("k*",0.5,3,10);
    std::array<ParticleType*,2> arr{&proton,&k};
    for(int i = 0; i != 2; ++i){
        arr[i]->Print();
        
    }*/
    Particle::AddParticleType("Protone",1,2,0);
    Particle::PrintData();

}