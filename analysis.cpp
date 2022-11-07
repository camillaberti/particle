#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include <iostream>

void analysis(){
    double const k_mass{0.89166};
    double const k_width{0.050};
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
    TF1* f1 = new TF1("f1", "pol0", 0, M_PI);
    f1->SetParameter(0, 1/M_PI);
    TF1* f2 = new TF1("f2","pol0",0, 2*M_PI);
    f2->SetParameter(0, 1/(2*M_PI));
    h[1]->Fit("f1","BQ");
    h[2]->Fit("f2","BQ");
    std::cout << '\n' << "Fitting informations of \u03b8 distribution: " << '\n';
    std::cout << "Parameter of const function: " << f1->GetParameter(0) << '\n';
    std::cout << "Chi-square/NDF: " << f1->GetChisquare()/f1->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f1->GetProb() << '\n'<<'\n';
    std::cout << "Fitting informations of \u03c6 distribution: " << '\n';
    std::cout << "Parameter of const function: " << f2->GetParameter(0) << '\n';
    std::cout << "Chi-square/NDF: " << f2->GetChisquare()/f1->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f2->GetProb() << '\n'<<'\n';
    TF1* f3 = new TF1("f3","expo",0,10);
    f3->SetParameters(0,-1);
    h[3]->Fit("f3","BQ");
    std::cout << "Fitting informations of P distribution: " << '\n';
    std::cout << "Parameters of fit function: " << f3->GetParameter(0) + '\t' + f3->GetParameter(1) << '\n';
    std::cout << "Chi-square/NDF: " << f3->GetChisquare()/f3->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f3->GetProb() << '\n'<<'\n';
    TH1F* hdif1 = new TH1F(*h[8]); //copy constructor of TH1F
    TH1F* hdif2 = new TH1F(*h[8]);
    hdif1->Add(h[7],h[8],1,-1);
    hdif2->Add(h[9],h[10],1,-1);
    hdif1->SetNameTitle("hdif1","Difference between histograms of mass invariant 1");
    hdif2->SetNameTitle("hdif2","Difference between histograms of mass invariant 2");
    TF1* f4 = new TF1("f4","gaus",0,2);
    f4->SetParameter(1,k_mass);
    f4->SetParameter(2,k_width);
    hdif1->Fit("f4","BQ");
    //hdif2->Fit("f4","BQ");
    std::cout << "Fitting informations of hdif1: " << '\n';
    std::cout << "Parameters of fit function: " << f4->GetParameter(0) + '\t' + f4->GetParameter(1) << '\n';
    std::cout << "Chi-square/NDF: " << f4->GetChisquare()/f4->GetNDF() << '\n';
    std::cout << "Chi-square prob: " << f4->GetProb() << '\n';


    file->Close();
}