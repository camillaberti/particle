#include "particletype.hpp"
#include "resonancetype.hpp"
#include "particle.hpp"
#include "TRandom.h"
#include "TH1F.h"

int main(){
    R__LOAD_LIBRARY(particletype_cpp.so);
    R__LOAD_LIBRARY(resonancetype_cpp.so);
    R__LOAD_LIBRARY(particle_cpp.so);
    R__LOAD_LIBRARY(extra_cpp.so);
    gRandom->SetSeed();
    Particle::AddParticleType("Pion+",0.13957,1);
    Particle::AddParticleType("Pion-",0.13957,-1);
    Particle::AddParticleType("Kaon+",0.49367,1);
    Particle::AddParticleType("Kaon-",0.49367,-1);
    Particle::AddParticleType("Proton+",0.93827,1);
    Particle::AddParticleType("Proton-",0.93827,-1);
    Particle::AddParticleType("K*",0.89166,0,0.05);

    TH1F* h1 = new TH1F("h1","Types of particles",7,0,7);
    TH1F* h2 = new TH1F("h2","Distribution of #theta",1000,0,M_PI);
    TH1F* h3 = new TH1F("h3","Distribution of #varphi",1000,0,2*M_PI);
    TH1F* h4 = new TH1F("h4","Distribution of P",500,0,10);
    TH1F* h5 = new TH1F("h5","Distribution of transverse P",500,0,10);
    TH1F* h6 = new TH1F("h6","Distribution of energy",500,0,10);





    std::vector<Particle> EventParticles;
    for(int i = 0; i != 1E5; ++i){
        for(int n = 0; n != 100; ++i){
            Particle particle{};
            double phi = gRandom->Rndm()*2*M_PI;
            double theta = gRandom->Rndm()*M_PI;
            double P = gRandom->Exp(1);
            particle.SetP(P*std::sin(theta)*std::cos(phi),P*std::sin(theta)*std::sin(phi),P*std::cos(theta));
            auto x = gRandom->Rndm();
            if(x<0.4){particle.SetIndex(0);}
            else if(x<0.8){particle.SetIndex(1);}
            else if(x<0.85){particle.SetIndex(2);}
            else if(x<0.9){particle.SetIndex(3);}
            else if(x<0.945){particle.SetIndex(4);}
            else if(x<0.99){particle.SetIndex(5);}
            else if(x<1){particle.SetIndex(6);}
            EventParticles.push_back(particle); 
            //filling histograms, non bisogna includere le figlie dei decadimenti  
            h1->Fill(particle.GetIndex()); 
            h2->Fill(theta);
            h3->Fill(phi); 
            h4->Fill(P);
            h5->Fill(sqrt(particle.GetPx()*particle.GetPx()+particle.GetPy()*particle.GetPy()));
            h6->Fill(particle.GetEnergy());    
        }
        for(auto& p : EventParticles){
            if(p.GetIndex() == 6){
                auto j = gRandom->Rndm();
                if(j < 0.5){
                    Particle p1{"Pion+"};
                    Particle p2{"Kaon-"};
                    p.Decay2body(p1,p2);
                    EventParticles.push_back(p1);
                    EventParticles.push_back(p2);
                }
                else {
                    Particle p1{"Pion-"};
                    Particle p2 ={"Kaon+"};
                    p.Decay2body(p1,p2);
                    EventParticles.push_back(p1);
                    EventParticles.push_back(p2);
                }
            }  
        }
        EventParticles.clear();
    }

}