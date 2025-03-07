#include <stdio.h>
#include <iostream>
#include <cmath>
#include "TApplication.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>

using namespace std;


//ficheiro para as corre    es e sistem  ticas

ifstream muon_neg18_lowTrig; 
ifstream muon_pos18_lowTrig; 
ifstream muon_18_highTrig; 
ifstream muon_IdISO;
ifstream ele_trig;
ifstream ele_idISO;
ifstream ele_RECO;

void openFiles(){
 muon_neg18_lowTrig.open("muon_neg18_lowTrig.txt");
 muon_pos18_lowTrig.open("muon_pos18_lowTrig.txt");
 muon_18_highTrig.open("muon_18_highTrig.txt");
 muon_IdISO.open("muon_IdISO.txt");
 ele_trig.open("ele_trig.txt");
 ele_idISO.open("ele_iso_id.txt");
 ele_RECO.open("ele_reco.txt");
}


//fun    es para calcular as corre    es

 double MuonTrigSF(double pt,double eta,double carga,double p_tot, int flag=0){


        double sf=0;

        muon_neg18_lowTrig.clear();
        muon_pos18_lowTrig.clear();
        muon_18_highTrig.clear();
        muon_neg18_lowTrig.seekg(0, std::ios::beg);
        muon_pos18_lowTrig.seekg(0, std::ios::beg);
        muon_18_highTrig.seekg(0, std::ios::beg);

        double eta_min, eta_max, pt_min, pt_max, corr, syst;

        if(pt<200. && carga<0.){

                while(muon_neg18_lowTrig >> eta_min >> eta_max >> pt_min >> pt_max >> corr >> syst){

                        if(abs(eta)>eta_min && abs(eta)<=eta_max && pt>pt_min && pt<=pt_max) {if(flag==0) sf=corr; if(flag==1) sf=syst;};

                }

        }

        muon_neg18_lowTrig.clear();
        muon_pos18_lowTrig.clear();
        muon_18_highTrig.clear();
        muon_neg18_lowTrig.seekg(0, std::ios::beg);
        muon_pos18_lowTrig.seekg(0, std::ios::beg);
        muon_18_highTrig.seekg(0, std::ios::beg);

        if(pt<200. && carga>0.){

                while(muon_pos18_lowTrig >> eta_min >> eta_max >> pt_min >> pt_max >> corr >> syst){

                        if(abs(eta)>eta_min && abs(eta)<=eta_max && pt>pt_min && pt<=pt_max) {if(flag==0) sf=corr; if(flag==1) sf=syst;};

                }

        }

	muon_neg18_lowTrig.clear();
        muon_pos18_lowTrig.clear();
        muon_18_highTrig.clear();
        muon_neg18_lowTrig.seekg(0, std::ios::beg);
        muon_pos18_lowTrig.seekg(0, std::ios::beg);
        muon_18_highTrig.seekg(0, std::ios::beg);

        if(pt>200.){

                while(muon_18_highTrig >> eta_min >> eta_max >> pt_min >> pt_max >> corr >> syst){

                        if(abs(eta)>eta_min && abs(eta)<=eta_max && pt>pt_min && pt<=pt_max) {if(flag==0) sf=corr; if(flag==1) sf=syst;};

                }

        }

	muon_neg18_lowTrig.seekg(0, std::ios::beg);
        muon_pos18_lowTrig.seekg(0, std::ios::beg);
        muon_18_highTrig.seekg(0, std::ios::beg);

       if(abs(eta)>2.4 || pt<35.) sf=1.;

        return sf;
}

double MuonIDISOSF(double pt,double eta,double carga,double p_tot, int flag=0){

        double sf=0;

        double eta_min, eta_max, pt_min, pt_max, corr, syst;

        muon_IdISO.clear();
        muon_IdISO.seekg(0, std::ios::beg);


        while(muon_IdISO >> eta_min >> eta_max >> pt_min >> pt_max >> corr >> syst){

              if(abs(eta)>eta_min && abs(eta)<=eta_max && pt>pt_min && pt<=pt_max) {if(flag==0) sf=corr; if(flag==1) sf=syst;}

        }


	if(abs(eta>2.4) || pt<35.) sf=1;

        return sf;


}

double muRecoSF(double pt, double eta, int flag=0){

        double sf=1;

        if(abs(eta)>2.1 && pt<120) sf=0.999;
        if(pt>=120){

                double p_tot = pt/sin(2*atan(exp(-eta)));
                if(p_tot>100 && p_tot < 150 && abs(eta)<1.6) sf=0.9948;
                if(p_tot>100 && p_tot < 150 && abs(eta)>1.6) sf=0.993;
                if(p_tot>150 && p_tot < 200 && abs(eta)<1.6) sf=0.9950;
                if(p_tot>150 && p_tot < 200 && abs(eta)>1.6) sf=0.990;
                if(p_tot>200 && p_tot < 300 && abs(eta)<1.6) sf=0.994;
                if(p_tot>200 && p_tot < 300 && abs(eta)>1.6) sf=0.988;
                if(p_tot>300 && p_tot < 400 && abs(eta)<1.6) sf=0.9914;
                if(p_tot>300 && p_tot < 400 && abs(eta)>1.6) sf=0.981;
                if(p_tot>400 && p_tot < 600 && abs(eta)<1.6) sf=0.993;
                if(p_tot>400 && p_tot < 600 && abs(eta)>1.6) sf=0.983;
                if(p_tot>600 && p_tot < 1500 && abs(eta)<1.6) sf=0.991;
                if(p_tot>600 && p_tot < 1500 && abs(eta)>1.6) sf=0.978;
                if(p_tot>1500 && p_tot < 3500 && abs(eta)<1.6) sf=1.;
                if(p_tot>1500 && p_tot < 3500 && abs(eta)>1.6) sf=0.98;

        }

	if(pt<200. && flag==1){

                if(abs(eta)<0.9) sf =0.000382;
                if(abs(eta)>=0.9 && abs(eta)<1.2) sf=0.000522;
                if(abs(eta)>=1.2 && abs(eta)<2.1) sf=0.000831;
                if(abs(eta)>=2.1) sf = 0.001724;

        }

	if(pt>200. && flag==1){

                double p_tot = pt/sin(2*atan(exp(-eta)));

                if(abs(eta)<1.6 && p_tot<300.) sf =0.001;
                if(abs(eta)<1.6 && p_tot>=300 && p_tot<400) sf=0.0009;
                if(abs(eta)<1.6 && p_tot>=400 && p_tot<600) sf=0.002;
                if(abs(eta)<1.6 && p_tot>=600 && p_tot<1500) sf = 0.004;
                if(abs(eta)<1.6 && p_tot>=1500) sf = 0.1;

                if(abs(eta)>1.6 && p_tot<300.) sf =0.001;
                if(abs(eta)>1.6 && p_tot>=300 && p_tot<400) sf=0.002;
                if(abs(eta)>1.6 && p_tot>=400 && p_tot<600) sf=0.003;
                if(abs(eta)>1.6 && p_tot>=600 && p_tot<1500) sf = 0.006;
                if(abs(eta)>1.6 && p_tot>=1500) sf = 0.03;

        }

	return sf;
}

double EleTrig(double pt,double eta, int flag=0){

        double sf=0;

        double eta_min, eta_max, pt_min, pt_max, corr, syst;

        ele_trig.clear();
        ele_trig.seekg(0, std::ios::beg);


        while(ele_trig >> eta_min >> eta_max >> pt_min >> pt_max >> corr >> syst){

              if(eta>eta_min && eta<=eta_max && pt>pt_min && pt<=pt_max) {if(flag==0) sf=corr; if(flag==1) sf=syst;}

        }


        if(abs(eta>2.4) || pt<35.) sf=1;

        return sf;


}


double EleIDISO(double pt,double eta, int flag =0){

        double sf=0;

        double eta_min, eta_max, pt_min, pt_max, corr, syst;

        ele_idISO.clear();
        ele_idISO.seekg(0, std::ios::beg);


        while(ele_idISO >> eta_min >> eta_max >> pt_min >> pt_max >> corr >> syst){

              if(eta>eta_min && eta<=eta_max && pt>pt_min && pt<=pt_max) {if(flag==0) sf=corr; if(flag==1) sf=syst;}

        }


	if(abs(eta>2.4) || pt<35.) sf=1;

        return sf;


}

double EleReco(double pt, double eta, int flag=0){


        double sf=0;

        double eta_min, eta_max, pt_min, pt_max, corr, syst;

        ele_RECO.clear();
        ele_RECO.seekg(0, std::ios::beg);


        while(ele_RECO >> eta_min >> eta_max >> pt_min >> pt_max >> corr >> syst){

              if(eta>eta_min && eta<=eta_max && pt>pt_min && pt<=pt_max) {if(flag==0) sf=corr; if(flag==1) sf=syst;}

        }


        if(abs(eta>2.4) || pt<35.) sf=1;

        return sf;

}



int main(){

	openFiles();

  //string prefix = "root:://cms-xrd-global.cern.ch//";
	string input;

	string prefix_output = "DY_fase1_PIC_skimmed_EMu_2018_";

	int k=0;

	double weight = 1.81;

	int numero_linha;
	cin >> numero_linha;

	stringstream ss;
	ss << numero_linha;

	string out_put;
	ss>>out_put;

	string output_tot = "/eos/user/m/mpisano/samples_2018_emu/fase1/" + prefix_output + out_put + ".root";

	ifstream ifile;
	ifile.open("DY_fase1_emu.txt");
	while(k<numero_linha) {
		ifile>>input;
		k++;
	}

	double entry_1,entry_2,entry_3,entry_4,entry_5,entry_6,entry_7,entry_8,entry_9,entry_10,entry_11,entry_12,entry_13,entry_14,entry_15,entry_16,entry_17,entry_18,entry_19,entry_20,entry_21,entry_22
	  ,entry_23,entry_24,entry_25,entry_26,entry_27,entry_28, entry_30, entry_31, entry_32, entry_33, entry_34, entry_35;
	int entry_29;
	//contadores
	double nid=0;
	double ncharge=0;
	double npt=0;

	////////////////////////////////////////
	TLorentzVector e;
	TLorentzVector mu;
	TLorentzVector sistema;

	string total = input;
	cout <<"Input file: " << total << endl;

	TApplication app("app", NULL, NULL);
	cout<<"joaonunes"<<endl;
	//TFile *f = TFile::Open("root:://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2520000/028FE21B-A107-4347-92C3-B533907C13DE.root");

	TFile *f = TFile::Open(total.c_str());
	cout<<"joaotavares"<<endl;

	TTree *tree = (TTree*) f->Get("tree");

	///A siguente linha tem de ser considerada apenas pelo fundo QCD
	TTree *tree_lumi = (TTree*) f->Get("LuminosityBlocks");
	////////////////////////////////////////////////////////////////////


	TFile output (output_tot.c_str(), "RECREATE", "");
	TTree out("tree","");


	out.Branch("e_id", &entry_1 ,"e_id/D");
	out.Branch("mu_id", &entry_2 ,"mu_id/D");
	out.Branch("e_pt", &entry_5 ,"e_pt/D");
	out.Branch("mu_pt", &entry_6, "mu_pt/D");
	out.Branch("e_charge", &entry_7, "e_charge/D");
	out.Branch("mu_charge", &entry_8, "mu_charge/D");
	out.Branch("e_eta", &entry_9, "e_eta/D");
	out.Branch("mu_eta", &entry_10, "mu_eta/D");
	out.Branch("e_phi", &entry_13, "e_phi/D");
	out.Branch("mu_phi", &entry_14, "mu_phi/D");
	out.Branch("e_mass", &entry_15, "e_mass/D");
	out.Branch("mu_mass", &entry_16, "mu_mass/D");
	out.Branch("sist_mass", &entry_17, "sist_mass/D");
	out.Branch("sist_acop", &entry_18, "sist_acop/D");
	out.Branch("sist_pt", &entry_19, "sist_pt/D");
	out.Branch("sist_rap", &entry_20, "sist_rap/D");
	out.Branch("met_pt", &entry_21, "met_pt/D");
	out.Branch("met_phi", &entry_22, "met_phi/D");
	out.Branch("jet_pt", &entry_23, "jet_pt/D");
	out.Branch("jet_eta", &entry_24, "jet_eta/D");
	out.Branch("jet_phi", &entry_25, "jet_phi/D");
	out.Branch("jet_mass", &entry_26, "jet_mass/D");
	out.Branch("jet_btag", &entry_27, "jet_btag/D");
	out.Branch("weight", &entry_28, "weight/D");
	out.Branch("n_b_jet",&entry_29, "n_b_jet/I");
	out.Branch("generator_weight", &entry_30, "generator_weight/D");
        out.Branch("syst_mu_trig", &entry_31, "syst_mu_trig/D");
        out.Branch("syst_mu_idiso", &entry_32, "syst_mu_idiso/D");
        out.Branch("syst_e_trig", &entry_33, "syst_e_trig/D");
        out.Branch("syst_e_idiso", &entry_34, "syst_e_idiso/D");
	out.Branch("syst_mu_reco", &entry_35, "syst_mu_reco/D");

	TH1D histo("histo","histo", 1000, 0, 1000);

	//TFile id_sf("TauID_SF_pt_DeepTau2017v2p1VSjet_UL2018.root");
        //TF1 *f_id_sf = (TF1*) id_sf.Get("VTight_cent");

        //TFile tauenergy("TauES_dm_DeepTau2017v2p1VSjet_2018ReReco.root");
        //TH1F *histtauenergy = (TH1F*) tauenergy.Get("tes");


	for(int i=0; i<tree->GetEntries(); i++){


		int eventos=tree->GetEvent(i);

		double k_e = 1;
                //double tau_id_sf=f_id_sf->Eval(tree->GetLeaf("tau0_pt")->GetValue(0))*f_id_sf->Eval(tree->GetLeaf("tau1_pt")->GetValue(0));

		e.SetPtEtaPhiM(tree->GetLeaf("e_pt")->GetValue(0),tree->GetLeaf("e_eta")->GetValue(0),tree->GetLeaf("electron_phi")->GetValue(0),tree->GetLeaf("electron_mass")->GetValue(0));
		mu.SetPtEtaPhiM(tree->GetLeaf("mu_pt")->GetValue(0),tree->GetLeaf("mu_eta")->GetValue(0),tree->GetLeaf("mu_phi")->GetValue(0),tree->GetLeaf("mu_mass")->GetValue(0));



		if(tree->GetLeaf("e_id")->GetValue(0)==1 && tree->GetLeaf("mu_id")->GetValue(0)==1  && e.DeltaR(mu)>0.4 && abs(e.Eta())<2.4 && abs(mu.Eta())<2.4){

		  double mu_trig_sf = MuonTrigSF(tree->GetLeaf("mu_pt")->GetValue(0), tree->GetLeaf("mu_pt")->GetValue(0), tree->GetLeaf("mu_charge")->GetValue(0), tree->GetLeaf("mu_pt")->GetValue(0)/sin(2.*atan(exp(-tree->GetLeaf("mu_eta")->GetValue(0)))));
                  double mu_idiso_sf = MuonIDISOSF(mu.Pt(), mu.Eta(), tree->GetLeaf("mu_charge")->GetValue(0), mu.Pt()/sin(2.*atan(exp(-mu.Eta()))));
                  double mu_reco_sf = muRecoSF(mu.Pt(), mu.Eta());
                  double ele_trig_sf= EleTrig(e.Pt(), e.Eta());
                  double ele_idiso_sf= EleIDISO(e.Pt(), e.Eta());
		  double ele_reco_sf = EleReco(e.Pt(), e.Eta());

		  nid=nid+weight*mu_trig_sf*mu_idiso_sf*mu_reco_sf*ele_trig_sf*ele_idiso_sf*ele_reco_sf;

		  if(tree->GetLeaf("e_pt")->GetValue(0)>35. && tree->GetLeaf("mu_pt")->GetValue(0)>35.){
		    npt=npt+weight*mu_trig_sf*mu_idiso_sf*mu_reco_sf*ele_trig_sf*ele_idiso_sf*ele_reco_sf;
		    
		    if(tree->GetLeaf("e_charge")->GetValue(0)*tree->GetLeaf("mu_charge")->GetValue(0)<0.){
		      ncharge=ncharge+weight*mu_trig_sf*mu_idiso_sf*mu_reco_sf*ele_trig_sf*ele_idiso_sf*ele_reco_sf;
				
										
		      double mu_trig_syst = MuonTrigSF(mu.Pt(), mu.Eta(), tree->GetLeaf("mu_charge")->GetValue(0), mu.Pt()/sin(2.*atan(exp(-mu.Eta()))),1);
                      double mu_idiso_syst = MuonIDISOSF(mu.Pt(), mu.Eta(), tree->GetLeaf("mu_charge")->GetValue(0), mu.Pt()/sin(2.*atan(exp(-mu.Eta()))),1);
                      double ele_trig_syst= EleTrig(e.Pt(), e.Eta(),1);
                      double ele_idiso_syst= EleIDISO(e.Pt(), e.Eta(),1);
		      double mu_reco_syst = muRecoSF(mu.Pt(), mu.Eta(),1);

                      entry_31 = mu_trig_syst;
                      entry_32 = mu_idiso_syst;
                      entry_33 = ele_trig_syst;
                      entry_34 = ele_idiso_syst;
		      entry_35 = mu_reco_syst;

		      entry_1=tree->GetLeaf("e_id")->GetValue(0);
		      entry_2=tree->GetLeaf("mu_id")->GetValue(0);
		      entry_5=tree->GetLeaf("e_pt")->GetValue(0);
		      entry_6=tree->GetLeaf("mu_pt")->GetValue(0);
		      entry_7=tree->GetLeaf("e_charge")->GetValue(0);
		      entry_8=tree->GetLeaf("mu_charge")->GetValue(0);
		      entry_9=tree->GetLeaf("e_eta")->GetValue(0);
		      entry_10=tree->GetLeaf("mu_eta")->GetValue(0);	  
		      entry_13=tree->GetLeaf("electron_phi")->GetValue(0);
		      entry_14=tree->GetLeaf("mu_phi")->GetValue(0);
		      entry_15=tree->GetLeaf("electron_mass")->GetValue(0);
		      entry_16=tree->GetLeaf("mu_mass")->GetValue(0);
										
		      e.SetPtEtaPhiM(entry_5,entry_9,entry_13,entry_15);
		      mu.SetPtEtaPhiM(entry_6,entry_10,entry_14,entry_16);

		      sistema=e+mu;
		      double Acoplanarity,deltaphi;
		      deltaphi=fabs(entry_14-entry_13);

		      if(deltaphi>M_PI)
			{
			  deltaphi=deltaphi-2*M_PI;
			}

		      Acoplanarity=fabs(deltaphi)/M_PI;
										
		      entry_17=sistema.M();
		      entry_18=Acoplanarity;
		      //cout<<entry_18<<endl;
		      entry_19=sistema.Pt();
		      entry_20=sistema.Rapidity();

		      entry_21=tree->GetLeaf("met_phi")->GetValue(0); // Energia dos neutrinos, não medida pelo detetor
		      entry_22=tree->GetLeaf("met_pt")->GetValue(0);
										
		      entry_23=tree->GetLeaf("jet_pt")->GetValue(0);
		      entry_24=tree->GetLeaf("jet_eta")->GetValue(0);
		      entry_25=tree->GetLeaf("jet_phi")->GetValue(0);
		      entry_26=tree->GetLeaf("jet_mass")->GetValue(0);
		      entry_27=tree->GetLeaf("jet_btag")->GetValue(0);
		      entry_28=weight*mu_trig_sf*mu_idiso_sf*mu_reco_sf*ele_trig_sf*ele_idiso_sf*ele_reco_sf;
		      entry_29=tree->GetLeaf("n_b_jet")->GetValue(0);
		      entry_30=tree->GetLeaf("generator_weight")->GetValue(0);

		      out.Fill();
		    }
		  }					
		}
		
				

		//if(i%1000==0) cout << "progress: " << double (i)/tree->GetEntries()*100 << endl; 	
					
	}
	cout<<"Eventos com Taus e Taus: "<<nid<<endl;
	cout<<"Eventos com carga oposta tautau: "<<ncharge<<endl;
	cout<<"Eventos apos corte no momento dos taus e taus: "<<npt<<endl;
	//histo.Draw("");
	cout << nid << endl;
	cout << ncharge << endl;
	cout << npt << endl;
	cout << -11111111111111 << endl;

	output.Write();
	//app.Run(true);

	return 0;


	}


