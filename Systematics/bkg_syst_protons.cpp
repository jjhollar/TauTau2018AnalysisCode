#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TGraphErrors.h"

using namespace std;

int main(){


	//pego no ficheiro de input
	TFile f("Fundo_total-protons_2018_tautau.root");
	TTree* tree = (TTree*) f.Get("tree");

	int n_evt = tree->GetEntries();

	//abro o ficheiro de output
	TFile out("Fundo_total-protons_2018_tautau_syst.root","RECREATE","");
	TTree tree_out("tree","");

        //Ficheiros para calcular as sistematicas

	TFile proton_sist_1("/afs/cern.ch/user/m/mpisano/PIC_joao/reco_charactersitics_version1.root");
        TGraphErrors* xi_sist_1 = (TGraphErrors*)proton_sist_1.Get("2018_TS1_TS2/multi rp-0/xi/g_systematics_vs_xi");
        TF1 xi_sist_inter_1("xi_sist_inter_1","pol20" ,0,10);
        xi_sist_inter_1.SetParLimits(0,-.8,0.1);
        xi_sist_inter_1.SetParLimits(1,0.1,1.);
        xi_sist_1->Fit("xi_sist_inter_1");
        TGraphErrors* xi_sist_2 = (TGraphErrors*)proton_sist_1.Get("2018_TS1_TS2/multi rp-1/xi/g_systematics_vs_xi");
        TF1 xi_sist_inter_2("xi_sist_inter_2","pol20",0,10);
        xi_sist_inter_2.SetParLimits(0,-1,0.1);
        xi_sist_inter_2.SetParLimits(1,0.1,1);
        xi_sist_2->Fit("xi_sist_inter_2");

	//declaração folhas da árvore de output

	double entry_1, entry_2, entry_3, entry_4, entry_5, entry_6, entry_7, entry_8, entry_9, entry_10, entry_11, entry_12, entry_13, entry_14, entry_15, entry_16;
	double entry_17, entry_18, entry_19, entry_20, entry_21, entry_22, entry_23, entry_24, entry_25, entry_26, entry_27, entry_28, entry_29, entry_30, entry_31;
	double entry_32, entry_33, entry_34, entry_35, entry_36, entry_37, entry_38;
	int entry_39;

	tree_out.Branch("tau0_id1", &entry_1, "tau0_id1/D");
	tree_out.Branch("tau0_id2", &entry_2, "tau0_id2/D");
	tree_out.Branch("tau0_id3", &entry_3, "tau0_id3/D");
	tree_out.Branch("tau1_id1", &entry_4, "tau1_id1/D");
        tree_out.Branch("tau1_id2", &entry_5, "tau1_id2/D");
        tree_out.Branch("tau1_id3", &entry_6, "tau1_id3/D");
	tree_out.Branch("tau0_pt", &entry_7, "tau0_pt/D");
	tree_out.Branch("tau1_pt", &entry_8, "tau1_pt/D");
	tree_out.Branch("tau0_charge", &entry_9, "tau0_charge/D");
        tree_out.Branch("tau1_charge", &entry_10, "tau1_charge/D");
	tree_out.Branch("tau0_eta", &entry_11, "tau0_eta/D");
        tree_out.Branch("tau1_eta", &entry_12, "tau1_eta/D");
	tree_out.Branch("tau_n", &entry_13, "tau_n/D");
	tree_out.Branch("tau0_phi", &entry_14, "tau0_phi/D");
        tree_out.Branch("tau1_phi", &entry_15, "tau1_phi/D");
	tree_out.Branch("tau0_mass", &entry_16, "tau0_mass/D");
        tree_out.Branch("tau1_mass", &entry_17, "tau1_mass/D");
	tree_out.Branch("sist_mass", &entry_18, "sist_mass/D");
	tree_out.Branch("sist_acop", &entry_19, "sist_acop/D");
	tree_out.Branch("sist_pt", &entry_20, "sist_pt/D");
	tree_out.Branch("sist_rap", &entry_21, "sist_rap/D");
	tree_out.Branch("met_pt", &entry_22, "met_pt/D");
	tree_out.Branch("met_phi", &entry_23, "met_phi/D");
	tree_out.Branch("jet_pt", &entry_24, "jet_pt/D");
	tree_out.Branch("jet_eta", &entry_25, "jet_eta/D");
	tree_out.Branch("jet_phi", &entry_26, "jet_phi/D");
	tree_out.Branch("jet_mass", &entry_27, "jet_mass/D");
	tree_out.Branch("jet_btag", &entry_28, "jet_btag/D");
	tree_out.Branch("weight", &entry_38, "weight/D");
	tree_out.Branch("N_b_jet", &entry_39, "N_b_jet/I");	
	tree_out.Branch("generator_weight", &entry_29, "generator_weight/D");
	tree_out.Branch("xi_arm1_1", &entry_30, "xi_arm1_1/D");
	tree_out.Branch("xi_arm1_2", &entry_31, "xi_arm1_2/D");
	tree_out.Branch("xi_arm2_1", &entry_32, "xi_arm2_1/D");
        tree_out.Branch("xi_arm2_2", &entry_33, "xi_arm2_2/D");
	tree_out.Branch("xi_arm1_1_up", &entry_34, "xi_arm1_1_up/D");
        tree_out.Branch("xi_arm1_1_dw", &entry_35, "xi_arm1_1_dw/D");
        tree_out.Branch("xi_arm2_1_up", &entry_36, "xi_arm2_1_up/D");
        tree_out.Branch("xi_arm2_1_dw", &entry_37, "xi_arm2_1_dw/D");

	//preencho a árvore

	for(int i=0; i<n_evt; i++){

		tree->GetEvent(i);

		entry_1=tree->GetLeaf("tau0_id1")->GetValue(0);
		entry_2=tree->GetLeaf("tau0_id2")->GetValue(0);
		entry_3=tree->GetLeaf("tau1_id3")->GetValue(0);
		entry_4=tree->GetLeaf("tau1_id1")->GetValue(0);
                entry_5=tree->GetLeaf("tau1_id2")->GetValue(0);
                entry_6=tree->GetLeaf("tau1_id3")->GetValue(0);
		entry_7=tree->GetLeaf("tau0_pt")->GetValue(0);
		entry_8=tree->GetLeaf("tau1_pt")->GetValue(0);
		entry_9=tree->GetLeaf("tau0_charge")->GetValue(0);
                entry_10=tree->GetLeaf("tau1_charge")->GetValue(0);
		entry_11=tree->GetLeaf("tau0_eta")->GetValue(0);
                entry_12=tree->GetLeaf("tau1_eta")->GetValue(0);
		entry_13=tree->GetLeaf("tau_n")->GetValue(0);
		entry_14=tree->GetLeaf("tau0_phi")->GetValue(0);
                entry_15=tree->GetLeaf("tau1_phi")->GetValue(0);
		entry_16=tree->GetLeaf("tau0_mass")->GetValue(0);
                entry_17=tree->GetLeaf("tau1_mass")->GetValue(0);
		entry_18=tree->GetLeaf("sist_mass")->GetValue(0);
		entry_19=tree->GetLeaf("sist_acop")->GetValue(0);
		entry_20=tree->GetLeaf("sist_pt")->GetValue(0);
		entry_21=tree->GetLeaf("sist_rap")->GetValue(0);
		entry_22=tree->GetLeaf("met_pt")->GetValue(0);
		entry_23=tree->GetLeaf("met_phi")->GetValue(0);
		entry_24=tree->GetLeaf("jet_pt")->GetValue(0);
		entry_25=tree->GetLeaf("jet_eta")->GetValue(0);
		entry_26=tree->GetLeaf("jet_phi")->GetValue(0);
		entry_27=tree->GetLeaf("jet_mass")->GetValue(0);
		entry_28=tree->GetLeaf("jet_btag")->GetValue(0);
		entry_38=tree->GetLeaf("weight")->GetValue(0);
		entry_39=tree->GetLeaf("N_b_jet")->GetValue(0);
		entry_29=tree->GetLeaf("generator_weight")->GetValue(0);
		entry_30=tree->GetLeaf("xi_arm1_1")->GetValue(0);
		entry_31=tree->GetLeaf("xi_arm1_2")->GetValue(0);
		entry_32=tree->GetLeaf("xi_arm2_1")->GetValue(0);
		entry_33=tree->GetLeaf("xi_arm2_2")->GetValue(0);
		entry_34=entry_30+xi_sist_inter_1.Eval(entry_30);
		entry_35=entry_30-xi_sist_inter_1.Eval(entry_30);
		entry_36=entry_32+xi_sist_inter_1.Eval(entry_32);
                entry_37=entry_32-xi_sist_inter_1.Eval(entry_32);

		
		tree_out.Fill();

	}

	out.Write();

	return 0;

}
