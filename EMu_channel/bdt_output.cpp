#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>
#include <THStack.h>
#include <TApplication.h>
#include <TAxis.h>
#include "iostream"

int main() {

    TApplication app("app",NULL,NULL);

    // Abrir os arquivos contendo os histogramas
    TFile *fileDY = TFile::Open("TMVApp_DY.root");
    TFile *fileSinal = TFile::Open("TMVApp_sinal.root");
    TFile *fileTTJets = TFile::Open("TMVApp_ttjets.root");
    TFile *fileDados = TFile::Open("TMVApp_dados.root");

    // Extrair os histogramas dos arquivos
    TH1F *histDY = (TH1F*)fileDY->Get("MVA_BDT");
    TH1F *histSinal = (TH1F*)fileSinal->Get("MVA_BDT");
    TH1F *histTTJets = (TH1F*)fileTTJets->Get("MVA_BDT");
    TH1F *histDados = (TH1F*)fileDados->Get("MVA_BDT");

    TH1F sum_bkg;

    histDY->Rebin(4); // Reagrupa em 20 bins (80/4 = 20)
    //histQCD->Rebin(4);
    histTTJets->Rebin(4);
    histSinal->Rebin(4);
    histDados->Rebin(4);

    histDY->Scale(46.8/histDY->Integral());
    histTTJets->Scale(137./histTTJets->Integral());
    histSinal->Scale(3.20e-2*5000/histSinal->Integral());
    histDados->Scale(183./histDados->Integral());

    std::cout << "Data integral = " << histDados->Integral() << std::endl;
    
    sum_bkg=*histDY;
    sum_bkg.Add(histTTJets);
    sum_bkg.SetFillColor(kGray+1);
    sum_bkg.SetFillStyle(3354);

        for (int i = 1; i <= sum_bkg.GetNbinsX(); ++i) {
           double content = sum_bkg.GetBinContent(i);
           double error = sqrt(content);
           sum_bkg.SetBinError(i, error);
        }


    // Criar um TCanvas para desenhar o histograma
    TCanvas *canvas = new TCanvas("canvas", "Stacked Histogram", 800, 600);

    // Criar um TStack para sobrepor os histogramas
    THStack *stack = new THStack("stack", "Stacked Histogram");

    // Definir cores para os histogramas
    histDY->SetFillColor(kYellow);
    //histQCD->SetFillColor(kRed);
    histTTJets->SetFillColor(kGreen);
    histDY->SetLineColor(kYellow);
    histSinal->SetLineColor(kBlack);
    histTTJets->SetLineColor(kGreen);
    histSinal->SetLineWidth(3);
    histDados->SetLineWidth(2);
    histDados->SetLineColor(kBlack);

    //histSinal->Scale(3000);

    // Adicionar os histogramas ao TStack
    stack->Add(histDY);
    //stack->Add(histQCD);
    stack->Add(histTTJets);

    // Desenhar o TStack
    //    stack->SetMinimum(0.01); stack->SetMaximum(5000.0);
    stack->Draw("histo"); // "nostack" para sobrepor
    histSinal->Draw("same && histo");
    histDados->Draw(" same && e");
    sum_bkg.Draw("same && e2");

    // Criar a legenda
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    //legend->AddEntry(histQCD, "QCD contribution", "f");
    legend->AddEntry(histDY, "DY contribution", "f");
    legend->AddEntry(histTTJets, "TTJets contribution", "f");
    legend->AddEntry(histSinal, "Signal (x5000)", "l");
    legend->AddEntry(histDados, "Data", "l");
    legend->SetBorderSize(1);
    legend->SetFillColor(0);
    legend->Draw();
    stack->GetXaxis()->SetTitle("BDT output");
    stack->GetYaxis()->SetTitle("Arbitrary units");

    // Mostrar o canvas
    canvas->Draw();

    app.Run("true");

}

