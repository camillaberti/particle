#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>

void setStyle(){
  //gROOT->SetStyle("Plain");
  //gStyle->SetPalette(57);
  //gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1110);
}

void analysis(){
    double const k_mass{0.89166};
    double const k_width{0.050};
    TH1::AddDirectory(kFALSE); 
    TFile* file = new TFile("histograms.root");
    file->cd(); 
    TH1F* h[12];
    TString histname = "h";
    for(int i = 0; i != 12; ++i){
        h[i] = (TH1F*)file->Get(histname + i);
        std::cout<< "Entries of " + histname + i+":" +'\t' << h[i]->Integral() << '\n';        
    }
    for(int i = 1; i != 8; ++i){
        std::cout << "Bin " << i << " with error: " << h[0]->GetBinContent(i) << "+/-" << h[0]->GetBinError(i) << '\n'; 
    }
    TF1* fpolar = new TF1("fpolar", "pol0", 0, M_PI);
    fpolar->SetLineColor(kRed);
    TF1* fazimutal = new TF1("fazimutal","pol0",0, 2*M_PI);
    fazimutal->SetLineColor(kRed); 
    h[1]->Fit("fpolar","BQ");
    h[2]->Fit("fazimutal","BQ");
    std::cout << '\n' << "Fitting informations of \u03b8 distribution: " << '\n';
    std::cout << "Parameter of fit function: " << fpolar->GetParameter(0) << '\n';
    std::cout << "Chi-square/NDF: " << fpolar->GetChisquare()/fpolar->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << fpolar->GetProb() << '\n'<<'\n';
    std::cout << "Fitting informations of \u03c6 distribution: " << '\n';
    std::cout << "Parameter of fit function: " << fazimutal->GetParameter(0) << '\n';
    std::cout << "Chi-square/NDF: " << fazimutal->GetChisquare()/fazimutal->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << fazimutal->GetProb() << '\n'<<'\n';
    TF1* fexp = new TF1("fexp","expo",0,10);
    fexp->SetParameters(0,-1);
    h[3]->Fit("fexp","BQ");
    std::cout << "Fitting informations of P distribution: " << '\n';
    std::cout << "Parameters of fit function: " << fexp->GetParameter(0) + '\t' + fexp->GetParameter(1) << '\n';
    std::cout << "Chi-square/NDF: " << fexp->GetChisquare()/fexp->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << fexp->GetProb() << '\n'<<'\n';

    TH1F* hdif1 = new TH1F(*h[11]); //copy constructor of TH1F
    TH1F* hdif2 = new TH1F(*h[11]);
    hdif1->Add(h[7],h[8],1,-1);
    hdif2->Add(h[9],h[10],1,-1);
    hdif1->SetNameTitle("hdif1","Mass invariant: opposite charge minus same charge");
    hdif2->SetNameTitle("hdif2","Mass invariant: (p+/k- and p-/k+) minus (p+/k+ + p-/k-)");
    hdif1->SetFillColor(kBlue);
    hdif2->SetFillColor(kBlue);
    TF1* f1 = new TF1("f1","gaus",0,4);
    f1->SetParameter(1,k_mass);
    f1->SetParameter(2,k_width);

    h[11]->Fit("f1","BQ");
    std::cout << "Fitting informations on mass invariant between decayed particles: " << '\n';
    std::cout << "Mean: " << f1->GetParameter(1) << '\n';
    std::cout << "\u03C3: " << f1->GetParameter(2) << '\n';
    std::cout << "Chi-square/NDF: " << f1->GetChisquare()/f1->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f1->GetProb() << '\n'<<'\n';
    
    
    std::cout<< "Entries of hdif1: " << hdif1->Integral() << '\n'; //non so se lasciare o meno
    std::cout<< "Entries of hdif2: " << hdif2->Integral() << '\n';
    TF1* f2 = new TF1("f2","gaus",0,4);
    f2->SetParameter(1,k_mass);
    f2->SetParameter(2,k_width);
    TF1* f3 = new TF1("f3","gaus",0,4);
    f3->SetParameter(1,k_mass);
    f3->SetParameter(2,k_width);
    hdif1->Fit("f2","BQ");
    hdif2->Fit("f3","BQ");
    std::cout << "Fitting informations of mass invariant: opposite charge minus same charge " << '\n';
    std::cout << "Mean: " << f2->GetParameter(1) << '\n';
    std::cout << "\u03C3: " << f2->GetParameter(2) << '\n';
    std::cout << "Chi-square/NDF: " << f2->GetChisquare()/f2->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f2->GetProb() << '\n'<<'\n';
    
    std::cout << "Fitting informations of mass invariant: (p+/k- and p-/k+) minus (p+/k+ + p-/k-)" << '\n';
    std::cout << "Mean: " << f3->GetParameter(1) << '\n';
    std::cout << "\u03C3: " << f3->GetParameter(2) << '\n';
    std::cout << "Chi-square/NDF: " << f3->GetChisquare()/f3->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f3->GetProb() << '\n'<<'\n';
    

    //cosmetics
    for(auto const& hist : h){
        hist->SetFillColor(kBlue);
        hist->GetYaxis()->SetTitle("Entries");
    }
    hdif1->SetFillColor(kBlue);
    hdif2->SetFillColor(kBlue);
    //drawing
    TCanvas* c1 = new TCanvas("c1","Distribution of particles, angles, impulse and energy",900,600);
    c1->Divide(2,3);
    for(int i = 0; i != 6; ++i){
        c1->cd(i+1);
        h[i]->Draw();
    }
    TCanvas* c2 = new TCanvas("c2","Invariant mass 1",900,600);
    c2->Divide(2,3);
    for(int i = 1; i != 6; ++i){
        c2->cd(i);
        h[i+5]->Draw("HISTO");
    }
    TCanvas* c3 = new TCanvas("c3","Invariant mass 2",900,600);
    c3->Divide(2,2);
    c3->cd(1);
    h[11]->Draw("H");
    h[11]->Draw("E,P,SAME");
    c3->cd(3);
    hdif1->Draw("H");
    hdif1->Draw("E,P,SAME");
    c3->cd(4);
    //hdif2->Draw("HISTO");
    hdif2->Draw("H");
    hdif2->Draw("E,P,SAME");


    file->Close();
}