#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP
#include<string>
class ParticleType{
    public:
    ParticleType(std::string name, double mass, int charge);
    std::string GetName() const{return name_;}
    double GetMass() const {return mass_;}
    int GetCharge() const {return charge_;}
    virtual void Print() const;
    virtual double GetWidth() const {return 0.0;}
    private:
    std::string const name_;
    double const mass_;
    int const charge_;
};
#endif