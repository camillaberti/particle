#include "particletype.hpp"
#include "resonancetype.hpp"
#include "particle.hpp"
#include "TRandom.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <algorithm>

void mymain(){
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

    TH1F* h0 = new TH1F("h0","Types of particles",7,0,7);
    TH1F* h1 = new TH1F("h1","Distribution of #theta",1000,0,M_PI);
    TH1F* h2 = new TH1F("h2","Distribution of #varphi",1000,0,2*M_PI);
    TH1F* h3 = new TH1F("h3","Distribution of P",500,0,10);
    TH1F* h4 = new TH1F("h4","Distribution of transverse P",500,0,10);
    TH1F* h5 = new TH1F("h5","Distribution of energy",500,0,10);
    TH1F* h6 = new TH1F("h6","Total invariant mass",160,0,2);
    TH1F* h7 = new TH1F("h7","Invariant mass with opposite charge",160,0,2);
    TH1F* h8 = new TH1F("h8","Invariant mass with same charge",160,0,2);
    TH1F* h9 = new TH1F("h9","Invariant mass with pion+/kaon- and pion-/kaon+",160,0,2);
    TH1F* h10 = new TH1F("h10","Invariant mass with pion+/kaon+ and pion-/kaon-",160,0,2);
    TH1F* h11 = new TH1F("h11","Invariant mass between decayed particles",80,0,2);
    TH1F* histograms[12]{h0,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11};
    
    h7->Sumw2();
    h8->Sumw2();
    h9->Sumw2();
    h10->Sumw2();

    std::vector<Particle> EventParticles;
    for(int i = 0; i != 1E5; ++i){
        for(int n = 0; n != 100; ++n){
            Particle particle{};
            double phi = gRandom->Uniform(0, 2*M_PI);
            double theta = gRandom->Uniform(0,M_PI);
            double P = gRandom->Exp(1);
            particle.SetP(P*std::sin(theta)*std::cos(phi),P*std::sin(theta)*std::sin(phi),P*std::cos(theta));
            auto const x = gRandom->Rndm();
            if(x<0.4){particle.SetIndex(0);}
            else if(x<0.8){particle.SetIndex(1);}
            else if(x<0.85){particle.SetIndex(2);}
            else if(x<0.9){particle.SetIndex(3);}
            else if(x<0.945){particle.SetIndex(4);}
            else if(x<0.99){particle.SetIndex(5);}
            else {particle.SetIndex(6);}
            EventParticles.push_back(particle); 
            //filling histograms, non bisogna includere le figlie dei decadimenti  
            h0->Fill(particle.GetIndex()); 
            h1->Fill(theta);
            h2->Fill(phi); 
            h3->Fill(P);
            h4->Fill(sqrt(particle.GetPx()*particle.GetPx()+particle.GetPy()*particle.GetPy()));
            h5->Fill(particle.GetEnergy());
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
                    h11->Fill(p1.InvMass(p2)); //solo per l'ultimo istogramma di minv bisogna considerare solo le "nuove" particelle
                }
                else {
                    Particle p1{"Pion-"};
                    Particle p2{"Kaon+"};
                    p.Decay2body(p1,p2);
                    EventParticles.push_back(p1);
                    EventParticles.push_back(p2);
                    h11->Fill(p1.InvMass(p2));
                }
            }  
        }
        //nella massa invariante bisogna includere le figlie dei decadimenti?
        auto it = EventParticles.begin();
        for(; it != EventParticles.end(); ++it){
            auto next = std::next(it);
            std::for_each(next,EventParticles.end(),[&](Particle p){ 
                if(it->GetParticleCharge() != 0 && next->GetParticleCharge() != 0){
                h6->Fill(it->InvMass(p));} //massa invariante fra tutte le particelle tranne le k
                if(it->GetParticleCharge() * p.GetParticleCharge() < 0){
                    h7->Fill(it->InvMass(p));
                }
                else {h8->Fill(it->InvMass(p));}
                if((it->GetIndex() == 0 && p.GetIndex() == 3) || (it->GetIndex() == 3 && p.GetIndex() == 0) 
                    || (it->GetIndex() == 1 && p.GetIndex() == 2) || (it->GetIndex() == 2 && p.GetIndex() == 1)){ //pion+/kaon- or pion-/kaon+ 
                    h9->Fill(it->InvMass(p));}
                if((it->GetIndex() == 0 && p.GetIndex() == 2) || (it->GetIndex() == 2 && p.GetIndex() == 0) 
                    || (it->GetIndex() == 1 && p.GetIndex() == 3) || (it->GetIndex() == 3 && p.GetIndex() == 1)){ //pion+/kaon+ or pion-/kaon- 
                    h10->Fill(it->InvMass(p));}        
                }); 

            }
        
        EventParticles.clear();
    } 
    TFile* file = new TFile("histograms.root","RECREATE");
    file->cd();
    for(auto& h : histograms){
        h->Write(); //uno per uno così sennò non funziona
    }
    
    file->Close();
    //drawing, questa parte poi va spostata
    TCanvas* c1 = new TCanvas("c1","test1",900,600);
    c1->Divide(2,3);
    for(int i = 0; i != 6; ++i){
        c1->cd(i+1);
        histograms[i]->Draw();
    }
    TCanvas* c2 = new TCanvas("c2","test inv mass",900,600);
    c2->Divide(2,3);
    for(int i = 1; i != 7; ++i){
        c2->cd(i);
        histograms[i+5]->Draw("HISTO");
    }
}