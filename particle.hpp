#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include "resonancetype.hpp"
#include <vector>
#include <string>

class Particle{
    public:
    //da mettere default constructor 
    Particle();
    Particle(std::string name, double p_x = 0, double p_y = 0, double p_z = 0);
    int GetIndex() const {return Index_;}
    double GetPx() const {return P_x_;}
    double GetPy() const {return P_y_;}
    double GetPz() const {return P_z_;}
    double GetParticleMass() const{return ParticleType_[Index_]->GetMass();}
    double GetEnergy() const;
    void SetIndex(int i);
    void SetIndex(std::string name);
    void SetP(double px,double py,double pz);
    double InvMass(Particle& p) const;
    static void AddParticleType(std::string name, double mass, double charge, double width = 0);
    static void PrintData();
    void PrintParticleInformation() const;
    int Decay2body(Particle &dau1,Particle &dau2) const;
    private:
    static std::vector <ParticleType*> ParticleType_;
    static const int MaxNumParticleType_ = 10; //(da mettere eccezione per proteggersi dall'utente pazzo)
    //static int NParticleType_; non so se serva
    int Index_;
    double P_x_;
    double P_y_;
    double P_z_;
    static int FindParticle(std::string name);
    void Boost(double bx, double by, double bz);
};
#endif