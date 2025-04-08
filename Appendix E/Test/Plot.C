#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include <TH1F.h>
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void StyleBase()
{
cout << "StyleBase started" << endl;
gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetNdivisions(210, "x");
gStyle->SetNdivisions(10, "y");
gStyle->SetLabelOffset(0.005, "x");
gStyle->SetLabelOffset(0.005, "y");
gStyle->SetLabelSize(0.05, "x");
gStyle->SetLabelSize(0.05, "y");
gStyle->SetTitleOffset(0.9, "x");
gStyle->SetTitleOffset(1.4, "y");
gStyle->SetTitleSize(0.05, "x");
gStyle->SetTitleSize(0.05, "y");
gStyle->SetPadBottomMargin(0.15);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadRightMargin(0.1);
gStyle->SetTickLength(0.05, "x");
gStyle->SetTickLength(0.03, "y");
gStyle->SetPadBorderMode(0);
}

// Open the ROOT file
TFile *file = new TFile("muon_models.root");

void Plot()
{
    StyleBase();

    char hist[11][3] = {"h0","h1","h2","h3","h4","h5","h6","h7", "h8", "h9"};

    for (int i=0; i<10; i++) {
        TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
        TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
        TH1F* histo;
        histo = (TH1F*)file->Get(hist[i]);
        histo->SetLineColor(kRed);
        histo->SetLineWidth(2);
        histo->SetXTitle("log10(Energy/[MeV])");

        // cross section ratio graph
        TLine *line = new TLine(-3, 1.0, 4, 1.0);
        line->SetLineColor(kBlue);
        line->SetLineStyle(2); // 2 for dotted line
        line->SetLineWidth(2);

        switch(i) {
            case 0:
                canvas->SetLogy();
                leg = new TLegend(0.7, 0.7, 0.9, 0.8);
                histo->Draw("HIST SAME L");
                histo->SetTitle("Cross section of#mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}/#mu^{+}Ae^{-}e^{+}");
                histo->SetYTitle("#sigma [g/cm2]");
                histo->GetYaxis()->SetTitleOffset(1.5);
                histo->GetYaxis()->SetRangeUser(1e-12,1e-2);
                histo->GetXaxis()->SetRangeUser(2,7);
                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}", "l");
                leg->Draw("SAME");

                histo = (TH1F*)file->Get(hist[i+1]);
                histo->SetLineColor(kGreen);
                histo->SetLineWidth(2);
                histo->Draw("HIST SAME L");

                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}Ae^{-}e^{+}", "l");
                leg->Draw("SAME");
                canvas->Print("muon_cross_section.png");

                break;
            case 1:
                break;
            case 2:
                break;
            case 3:
                leg = new TLegend(0.6, 0.8, 0.9, 0.9);
                histo->Draw("HIST SAME L");
                histo->SetTitle("#sigma ratio of (#mu^{+}A#rightarrow#mu^{+}Ae^{-}e^{+})/(#mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+})");
                histo->SetYTitle("#sigma ratio");
                histo->GetYaxis()->SetTitleOffset(1.5);
                histo->GetXaxis()->SetRangeUser(2,5);
                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}/#mu^{+}Ae^{-}e^{+}", "l");
                leg->Draw();
                //line->Draw();

                canvas->Print("muon_cross_section_ratio.png");
                break;
            case 4:
                histo->SetFillColorAlpha(kRed, 0.2);
                histo->Draw("HIST SAME");
                histo->SetTitle("Muon energy w.r.t. beam energy of #mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}/#mu^{+}Ae^{-}e^{+}");
                histo->SetYTitle("Events");
                //histo->GetYaxis()->SetTitleOffset(1.5);
                histo->SetXTitle("log10(#mu/E)");
                histo->GetXaxis()->SetRangeUser(-6,0);
                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}", "f");
                leg->Draw("SAME");

                histo = (TH1F*)file->Get(hist[i+1]);
                histo->SetFillColorAlpha(kGreen, 0.2);
                histo->Draw("HIST SAME");
                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}Ae^{-}e^{+}", "f");
                leg->Draw("SAME");

                canvas->Print("muon_energy_ratio.png");

                break;
            case 5:
                break;
            case 6:
                canvas->SetLogy();
                histo->SetFillColorAlpha(kRed, 0.2);
                histo->Draw("HIST SAME");
                histo->SetTitle("Muon cos#theta of #mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}/#mu^{+}Ae^{-}e^{+}");
                histo->SetYTitle("Events");
                histo->GetYaxis()->SetTitleOffset(1.5);
                histo->SetXTitle("cos#theta");
                histo->GetXaxis()->SetRangeUser(0,1);
                //histo->GetYaxis()->SetRangeUser(0,10000);
                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}", "f");
                leg->Draw("SAME");

                histo = (TH1F*)file->Get(hist[i+1]);
                histo->SetFillColorAlpha(kGreen, 0.2);
                histo->Draw("HIST SAME");
                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}Ae^{-}e^{+}", "f");
                leg->Draw("SAME");

                canvas->Print("muon_cos_theta.png");
                break;
            case 7:
                break;
            case 8:
                canvas->SetLogy();
                histo->SetFillColorAlpha(kRed, 0.2);
                histo->Draw("HIST SAME");
                histo->SetTitle("Relative angle of #mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}/#mu^{+}Ae^{-}e^{+}");
                histo->SetYTitle("Events");
                histo->GetYaxis()->SetTitleOffset(1.5);
                histo->SetXTitle("cos#theta");
                histo->GetXaxis()->SetRangeUser(0.5,1);
                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}A#mu^{-}#mu^{+}", "f");
                leg->Draw("SAME");

                histo = (TH1F*)file->Get(hist[i+1]);
                histo->SetFillColorAlpha(kGreen, 0.2);
                histo->Draw("HIST SAME");
                leg->AddEntry(histo, "#mu^{+}A#rightarrow#mu^{+}Ae^{-}e^{+}", "f");
                leg->Draw("SAME");

                canvas->Print("muon_relative_theta.png");
                break;
            case 9:
                break;
            }
    }

};

// execute macro file
// root -b -q Plot.C
