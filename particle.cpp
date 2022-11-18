#include "particle.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

std::vector <ParticleType*> Particle::ParticleType_{}; //initialization of static vector

int Particle::FindParticle(std::string name){
    auto it = std::find_if(ParticleType_.begin(),ParticleType_.end(),[&](ParticleType* p){return p->GetName() == name;});
    auto i = std::distance(ParticleType_.begin(),it);
    return i;
}

Particle::Particle(){
    Index_ = -1;
}

Particle::Particle(std::string name, double p_x, double p_y, double p_z): P_x_{p_x}, P_y_{p_y}, P_z_{p_z} {
    if (FindParticle(name) == static_cast<int>(ParticleType_.size())){
        std::cout << "No matches found" << '\n';
    }
    else {
        Index_ = FindParticle(name);
    }

}

void Particle::AddParticleType(std::string name, double mass, double charge, double width){
    ParticleType* p;
    (width == 0)? p = new ParticleType(name,mass,charge) : p = new ResonanceType(name,mass,charge,width);

    if (ParticleType_.size() == MaxNumParticleType_){
        std::cout << "Cannot add more particle types" << '\n';} //eventualmente lanciare un'eccezione?    
    else if(FindParticle(p->GetName()) != static_cast<int>(ParticleType_.size())){
        std::cout << "The type has been already added" << '\n';
    }
    else {
        ParticleType_.push_back(p);}
}

void Particle::SetIndex(int i){
    if(i != MaxNumParticleType_){
        Index_ = i;
    }
    else{
        std::cout<< "Can't reassign index: invalid value" << '\n';
    }
}

void Particle::SetIndex(std::string name){
    if(FindParticle(name) == static_cast<int>(ParticleType_.size())){
        std::cout<< "Can't reassign index: invalid value" << '\n';
    }
    else {
        Index_ = FindParticle(name);
    }

}

void Particle::PrintData(){
    std::for_each(ParticleType_.begin(),ParticleType_.end(),[](ParticleType* particle){particle->Print();});
}

void Particle::PrintParticleInformation() const{
    std::cout << "Index: " << Index_ <<'\n';
    std::cout << "Name: " << ParticleType_[Index_]->GetName() << '\n';
    std::cout <<"P_x, P_y, P_z: " << P_x_ << '\t' << P_y_ << '\t' << P_z_;
}

double Particle::GetEnergy() const{
    return sqrt(ParticleType_[Index_]->GetMass()*ParticleType_[Index_]->GetMass()+P_x_ * P_x_ + P_y_* P_y_ + P_z_ * P_z_);
}

double Particle::InvMass(Particle& p) const{
    return sqrt((GetEnergy() + p.GetEnergy())*(GetEnergy() + p.GetEnergy())- s(P_x_ + p.GetPx())*(P_x_ + p.GetPx()) - (P_y_ + p.GetPy())*(P_y_ + p.GetPy())-(P_z_ + p.GetPz())*(P_z_ + p.GetPz()));
}

void Particle::SetP(double px,double py,double pz){
    P_x_ = px;
    P_y_ = py;
    P_z_ = pz;
}