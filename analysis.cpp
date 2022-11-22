#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>

void setStyle(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);
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
        std::cout<< "Entries of " + histname + i+":" +'\t' << h[i]->GetEntries() << '\n';      
    }
    for(int i = 1; i != 8; ++i){
        std::cout << "Bin " << i << " with error: " << h[0]->GetBinContent(i) << " +/- " << h[0]->GetBinError(i) << '\n'; 
    }
    TF1* fpolar = new TF1("fpolar", "pol0", 0, M_PI);
    fpolar->SetLineColor(kBlue-2);
    TF1* fazimutal = new TF1("fazimutal","pol0",0, 2*M_PI);
    fazimutal->SetLineColor(kBlue-2); 
    h[1]->Fit("fpolar","BQ");
    h[2]->Fit("fazimutal","BQ");
    std::cout << '\n' << "Fitting informations of \u03b8 distribution: " << '\n';
    std::cout << "Parameter of fit function: " << fpolar->GetParameter(0) << " +/- " << fpolar->GetParError(0) << '\n';
    std::cout << "Chi-square/NDF: " << fpolar->GetChisquare()/fpolar->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << fpolar->GetProb() << '\n'<<'\n';
    std::cout << "Fitting informations of \u03c6 distribution: " << '\n';
    std::cout << "Parameter of fit function: " << fazimutal->GetParameter(0) << " +/- " << fazimutal->GetParError(0) << '\n';
    std::cout << "Chi-square/NDF: " << fazimutal->GetChisquare()/fazimutal->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << fazimutal->GetProb() << '\n'<<'\n';
    TF1* fexp = new TF1("fexp","expo",0,10);
    fexp->SetParameters(0,-1);
    fexp->SetLineColor(kBlue-2);
    h[3]->Fit("fexp","BQ");
    std::cout << "Fitting informations of P distribution: " << '\n';
    std::cout << "Parameter of fit function: " << fexp->GetParameter(1) << " +/- " << fexp->GetParError(1) << '\n';
    std::cout << "Chi-square/NDF: " << fexp->GetChisquare()/fexp->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << fexp->GetProb() << '\n'<<'\n';

    TH1F* hdif1 = new TH1F(*h[11]); //copy constructor of TH1F
    TH1F* hdif2 = new TH1F(*h[11]);
    hdif1->Add(h[7],h[8],1,-1);
    hdif2->Add(h[9],h[10],1,-1);
    hdif1->SetNameTitle("hdif1","Mass invariant: opposite charge minus same charge");
    hdif2->SetNameTitle("hdif2","Mass invariant: (p+/k- and p-/k+) minus (p+/k+ + p-/k-)");
    hdif1->SetFillColor(kPink-3);
    hdif2->SetFillColor(kPink-3);
    TF1* f1 = new TF1("f1","gaus",0,4);
    f1->SetParameter(1,k_mass);
    f1->SetParameter(2,k_width);
    f1->SetLineColor(kBlue-2);

    h[11]->Fit("f1","BQ");
    std::cout << "Fitting informations on mass invariant between decayed particles: " << '\n';
    std::cout << "Mean: " << f1->GetParameter(1) << " +/- " << f1->GetParError(1)<< '\n';
    std::cout << "\u03C3: " << f1->GetParameter(2) << " +/- " << f1->GetParError(2) << '\n';
    std::cout << "Width: " << f1->GetParameter(0) << " +/- " << f1->GetParError(0) << '\n';
    std::cout << "Chi-square/NDF: " << f1->GetChisquare()/f1->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f1->GetProb() << '\n'<<'\n';
    
    
    std::cout<< "Entries of hdif1: " << hdif1->Integral() << '\n'; 
    std::cout<< "Entries of hdif2: " << hdif2->Integral() << '\n';
    TF1* f2 = new TF1("f2","gaus",0,4);
    f2->SetParameter(1,k_mass);
    f2->SetParameter(2,k_width);
    f2->SetLineColor(kBlue-2);
    TF1* f3 = new TF1("f3","gaus",0,4);
    f3->SetParameter(1,k_mass);
    f3->SetParameter(2,k_width);
    f3->SetLineColor(kBlue-2);
    hdif1->Fit("f2","BQ");
    hdif2->Fit("f3","BQ");
    std::cout << "Fitting informations of mass invariant: opposite charge minus same charge " << '\n';
    std::cout << "Mean: " << f2->GetParameter(1) << " +/- " << f2->GetParError(1)<< '\n';
    std::cout << "\u03C3: " << f2->GetParameter(2) << " +/- " << f2->GetParError(2)<< '\n';
    std::cout << "Width: " << f2->GetParameter(0) << " +/- " << f2->GetParError(0) << '\n';
    std::cout << "Chi-square/NDF: " << f2->GetChisquare()/f2->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f2->GetProb() << '\n'<<'\n';
    
    std::cout << "Fitting informations of mass invariant: (p+/k- and p-/k+) minus (p+/k+ + p-/k-)" << '\n';
    std::cout << "Mean: " << f3->GetParameter(1) << " +/- " << f3->GetParError(1) << '\n';
    std::cout << "\u03C3: " << f3->GetParameter(2) << " +/- " << f3->GetParError(2) << '\n';
    std::cout << "Width: " << f3->GetParameter(0) << " +/- " << f3->GetParError(0) << '\n';
    std::cout << "Chi-square/NDF: " << f3->GetChisquare()/f3->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f3->GetProb() << '\n'<<'\n';
    

    //cosmetics
    for(auto const& hist : h){
        hist->SetFillColor(kPink-3);
        hist->GetYaxis()->SetTitle("Entries");
    }

    h[0]->GetXaxis()->SetBinLabel(1, "#pi +");
    h[0]->GetXaxis()->SetBinLabel(2, "#pi -");
    h[0]->GetXaxis()->SetBinLabel(3, "k+");
    h[0]->GetXaxis()->SetBinLabel(4, "k -");
    h[0]->GetXaxis()->SetBinLabel(5, "p +");
    h[0]->GetXaxis()->SetBinLabel(6, "p -");
    h[0]->GetXaxis()->SetBinLabel(7, "k*");
    h[0]->GetXaxis()->SetLabelSize(0.06);
    h[0]->GetXaxis()->SetTitle("Particle types");
    h[1]->GetXaxis()->SetTitle("#theta (rad)");
    h[2]->GetXaxis()->SetTitle("#Phi (rad)");
    h[3]->GetXaxis()->SetTitle("Impulse (GeV)");
    h[4]->GetXaxis()->SetTitle("Transverse impulse (GeV)");
    h[5]->GetXaxis()->SetTitle("Energy (GeV)");
    for(int i = 6; i != 12; ++i){
        h[i]->GetXaxis()->SetTitle("Invariant mass (GeV)");
    }
    h[1]->SetMinimum(5000);
    h[2]->SetMinimum(5000);
    h[1]->SetMaximum(15000);
    h[2]->SetMaximum(15000);
    hdif1->SetFillColor(kPink-3);
    hdif2->SetFillColor(kPink-3);
    hdif1->GetXaxis()->SetTitle("Invariant mass (GeV)");
    hdif1->GetYaxis()->SetTitle("Entries");
    hdif2->GetXaxis()->SetTitle("Invariant mass (GeV)");
    hdif2->GetYaxis()->SetTitle("Entries");
    //drawing
    TCanvas* c1 = new TCanvas("c1","Distribution of particles, angles, impulse",900,600);
    c1->Divide(2,2);
    for(int i = 0; i != 4; ++i){
        c1->cd(i+1);
        h[i]->Draw();
    }

    TCanvas* c2 = new TCanvas("c2","Transverse impulse and energy",900,600);
    c2->Divide(1,2);
    c2->cd(1);
    h[4]->Draw();
    c2->cd(2);
    h[5]->Draw();

    TCanvas* c3 = new TCanvas("c3","Invariant mass 1",900,600);
    c3->Divide(2,3);
    for(int i = 1; i != 6; ++i){
        c3->cd(i);
        h[i+5]->Draw("HISTO");
    }

    TCanvas* c4 = new TCanvas("c4","Invariant mass 2",900,600);
    c4->Divide(2,2);
    c4->cd(1);
    h[11]->Draw("H");
    h[11]->Draw("E,P,SAME");
    c4->cd(3);
    hdif1->Draw("H");
    hdif1->Draw("E,P,SAME");
    c4->cd(4);
    hdif2->Draw("H");
    hdif2->Draw("E,P,SAME");
    //saving canvas
    c1->Print("canvas1.root");
    c1->Print("canvas1.C");
    c1->Print("canvas1.pdf");
    c2->Print("canvas2.root");
    c2->Print("canvas2.C");
    c2->Print("canvas2.pdf");
    c3->Print("canvas3.root");
    c3->Print("canvas3.C");
    c3->Print("canvas3.pdf");
    c4->Print("canvas4.root");
    c4->Print("canvas4.C");
    c4->Print("canvas4.pdf");
    file->Close();
}