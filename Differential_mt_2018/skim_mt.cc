#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "makeHisto.h"
#include "mutau_Tree_mt.h"

int main(int argc, char** argv) {

    using namespace std;
    myMap1 = new map<string, TH1F*>();
    myMap2 = new map<string, TH2F*>();
    string status_sample = *(argv + 1);
    bool isMC = false;
    bool isEmbedded = false;
    bool isData = false;
    if (status_sample.compare("mc") == 0) isMC = true;
    if (status_sample.compare("embedded") == 0) isEmbedded = true;
    if (status_sample.compare("data") == 0) isData = true;
    string out = *(argv + 2);
    string outname= out;
    TFile *fout = TFile::Open(outname.c_str(), "RECREATE");

    string in = *(argv + 3);
    string inname= in;
    TFile *fIn = TFile::Open(inname.c_str());

    int recoil=0;
    string recoilType = *(argv + 4);
    if (recoilType.compare("W") == 0)  recoil=1;
    if (recoilType.compare("Z") == 0)  recoil=2;

    TTree* treePtr = (TTree*) fIn->Get("mt/final/Ntuple");
    TH1F *evCounter = (TH1F*) fIn->Get("mt/eventCount");
    TH1F *evCounterW = (TH1F*) fIn->Get("mt/summedWeights");
    HTauTauTree_mt* tree = new HTauTauTree_mt (treePtr);

    TTree *Run_Tree = new TTree("mutau_tree", "mutau_tree");
    Run_Tree->SetDirectory(0);
    Run_Tree->Branch("run", &run, "run/I");
    Run_Tree->Branch("lumi", &lumi, "lumi/I");
    Run_Tree->Branch("evt", &evt, "evt/I");

    Run_Tree->Branch("lheweight_muR0p5_muF0p5", &lheweight_muR0p5_muF0p5, "lheweight_muR0p5_muF0p5/F");
    Run_Tree->Branch("lheweight_muR0p5_muF1", &lheweight_muR0p5_muF1, "lheweight_muR0p5_muF1/F");
    Run_Tree->Branch("lheweight_muR0p5_muF2", &lheweight_muR0p5_muF2, "lheweight_muR0p5_muF2/F");
    Run_Tree->Branch("lheweight_muR1_muF0p5", &lheweight_muR1_muF0p5, "lheweight_muR1_muF0p5/F");
    Run_Tree->Branch("lheweight_muR1_muF2", &lheweight_muR1_muF2, "lheweight_muR1_muF2/F");
    Run_Tree->Branch("lheweight_muR2_muF0p5", &lheweight_muR2_muF0p5, "lheweight_muR2_muF0p5/F");
    Run_Tree->Branch("lheweight_muR2_muF1", &lheweight_muR2_muF1, "lheweight_muR2_muF1/F");
    Run_Tree->Branch("lheweight_muR2_muF2", &lheweight_muR2_muF2, "lheweight_muR2_muF2/F");

    Run_Tree->Branch("PythiaWeight_fsr_muR0p25", &PythiaWeight_fsr_muR0p25, "PythiaWeight_fsr_muR0p25/F");
    Run_Tree->Branch("PythiaWeight_fsr_muR0p5", &PythiaWeight_fsr_muR0p5, "PythiaWeight_fsr_muR0p5/F");
    Run_Tree->Branch("PythiaWeight_fsr_muR2", &PythiaWeight_fsr_muR2, "PythiaWeight_fsr_muR2/F");
    Run_Tree->Branch("PythiaWeight_fsr_muR4", &PythiaWeight_fsr_muR4, "PythiaWeight_fsr_muR4/F");
    Run_Tree->Branch("PythiaWeight_isr_muR0p25", &PythiaWeight_isr_muR0p25, "PythiaWeight_isr_muR0p25/F");
    Run_Tree->Branch("PythiaWeight_isr_muR0p5", &PythiaWeight_isr_muR0p5, "PythiaWeight_isr_muR0p5/F");
    Run_Tree->Branch("PythiaWeight_isr_muR2", &PythiaWeight_isr_muR2, "PythiaWeight_isr_muR2/F");
    Run_Tree->Branch("PythiaWeight_isr_muR4", &PythiaWeight_isr_muR4, "PythiaWeight_isr_muR4/F");

    Run_Tree->Branch("HTTgenfinalstate", &HTTgenfinalstate, "HTTgenfinalstate/F");
    Run_Tree->Branch("gen_met_pt", &gen_met_pt, "gen_met_pt/F");
    Run_Tree->Branch("gen_met_phi", &gen_met_phi, "gen_met_phi/F");
    Run_Tree->Branch("gen_mu_pt", &gen_mu_pt, "gen_mu_pt/F");
    Run_Tree->Branch("gen_mu_eta", &gen_mu_eta, "gen_mu_eta/F");
    Run_Tree->Branch("gen_mu_phi", &gen_mu_phi, "gen_mu_phi/F");
    Run_Tree->Branch("gen_tauh_pt", &gen_tauh_pt, "gen_tauh_pt/F");
    Run_Tree->Branch("gen_tauh_eta", &gen_tauh_eta, "gen_tauh_eta/F");
    Run_Tree->Branch("gen_tauh_phi", &gen_tauh_phi, "gen_tauh_phi/F");

    Run_Tree->Branch("matchEmbFilter_Mu20Tau27_1", &matchEmbFilter_Mu20Tau27_1, "matchEmbFilter_Mu20Tau27_1/F");
    Run_Tree->Branch("matchEmbFilter_Mu24_1", &matchEmbFilter_Mu24_1, "matchEmbFilter_Mu24_1/F");
    Run_Tree->Branch("matchEmbFilter_Mu27_1", &matchEmbFilter_Mu27_1, "matchEmbFilter_Mu27_1/F");
    Run_Tree->Branch("matchEmbFilter_Mu20Tau27_2", &matchEmbFilter_Mu20Tau27_2, "matchEmbFilter_Mu20Tau27_2/F");
    Run_Tree->Branch("matchEmbFilter_Mu20HPSTau27_2", &matchEmbFilter_Mu20HPSTau27_2, "matchEmbFilter_Mu20HPSTau27_2/F");

    Run_Tree->Branch("genpX", &genpX, "genpX/F");
    Run_Tree->Branch("genpY", &genpY, "genpY/F");
    Run_Tree->Branch("genM", &genM, "genM/F");
    Run_Tree->Branch("genpT", &genpT, "genpT/F");
    Run_Tree->Branch("vispX", &vispX, "vispX/F");
    Run_Tree->Branch("vispY", &vispY, "vispY/F");

    Run_Tree->Branch("Rivet_VEta", &Rivet_VEta, "Rivet_VEta/F");
    Run_Tree->Branch("Rivet_VPt", &Rivet_VPt, "Rivet_VPt/F");
    Run_Tree->Branch("Rivet_errorCode", &Rivet_errorCode, "Rivet_errorCode/F");
    Run_Tree->Branch("Rivet_higgsEta", &Rivet_higgsEta, "Rivet_higgsEta/F");
    Run_Tree->Branch("Rivet_higgsPt", &Rivet_higgsPt, "Rivet_higgsPt/F");
    Run_Tree->Branch("Rivet_nJets25", &Rivet_nJets25, "Rivet_nJets25/F");
    Run_Tree->Branch("Rivet_nJets30", &Rivet_nJets30, "Rivet_nJets30/F");
    Run_Tree->Branch("Rivet_p4decay_VEta", &Rivet_p4decay_VEta, "Rivet_p4decay_VEta/F");
    Run_Tree->Branch("Rivet_p4decay_VPt", &Rivet_p4decay_VPt, "Rivet_p4decay_VPt/F");
    Run_Tree->Branch("Rivet_prodMode", &Rivet_prodMode, "Rivet_prodMode/F");
    Run_Tree->Branch("Rivet_stage0_cat", &Rivet_stage0_cat, "Rivet_stage0_cat/F");
    Run_Tree->Branch("Rivet_stage1_1_fine_cat_pTjet30GeV", &Rivet_stage1_1_fine_cat_pTjet30GeV, "Rivet_stage1_1_fine_cat_pTjet30GeV/F");
    Run_Tree->Branch("Rivet_stage1_1_cat_pTjet30GeV", &Rivet_stage1_1_cat_pTjet30GeV, "Rivet_stage1_1_cat_pTjet30GeV/F");
    Run_Tree->Branch("Rivet_stage1_cat_pTjet30GeV", &Rivet_stage1_cat_pTjet30GeV, "Rivet_stage1_cat_pTjet30GeV/F");
    Run_Tree->Branch("Rivet_j1pt", &Rivet_j1pt, "Rivet_j1pt/F");
    Run_Tree->Branch("Rivet_j2pt", &Rivet_j2pt, "Rivet_j2pt/F");
    Run_Tree->Branch("Rivet_j1eta", &Rivet_j1eta, "Rivet_j1eta/F");
    Run_Tree->Branch("Rivet_j2eta", &Rivet_j2eta, "Rivet_j2eta/F");
    Run_Tree->Branch("Rivet_j1phi", &Rivet_j1phi, "Rivet_j1phi/F");
    Run_Tree->Branch("Rivet_j2phi", &Rivet_j2phi, "Rivet_j2phi/F");
    Run_Tree->Branch("Rivet_j1m", &Rivet_j1m, "Rivet_j1m/F");
    Run_Tree->Branch("Rivet_j2m", &Rivet_j2m, "Rivet_j2m/F");
    Run_Tree->Branch("Rivet_higgsRapidity", &Rivet_higgsRapidity, "Rivet_higgsRapidity/F");
    Run_Tree->Branch("gentau1_eta", &gentau1_eta, "gentau1_eta/F");
    Run_Tree->Branch("gentau2_eta", &gentau2_eta, "gentau2_eta/F");
    Run_Tree->Branch("gentau1_pt", &gentau1_pt, "gentau1_pt/F");
    Run_Tree->Branch("gentau2_pt", &gentau2_pt, "gentau2_pt/F");

    Run_Tree->Branch("npv", &npv, "npv/F");
    Run_Tree->Branch("npu", &npu, "npu/F");
    Run_Tree->Branch("L1iso", &L1iso, "L1iso/F");
    Run_Tree->Branch("L1pt", &L1pt, "L1pt/F");

    Run_Tree->Branch("pt_1", &pt_1, "pt_1/F");
    Run_Tree->Branch("phi_1", &phi_1, "phi_1/F");
    Run_Tree->Branch("eta_1", &eta_1, "eta_1/F");
    Run_Tree->Branch("m_1", &m_1, "m_1/F");
    Run_Tree->Branch("e_1", &e_1, "e_1/F");
    Run_Tree->Branch("q_1", &q_1, "q_1/F");
    Run_Tree->Branch("iso_1", &iso_1, "iso_1/F");

    Run_Tree->Branch("pt_2", &pt_2, "pt_2/F");
    Run_Tree->Branch("phi_2", &phi_2, "phi_2/F");
    Run_Tree->Branch("eta_2", &eta_2, "eta_2/F");
    Run_Tree->Branch("m_2", &m_2, "m_2/F");
    Run_Tree->Branch("e_2", &e_2, "e_2/F");
    Run_Tree->Branch("q_2", &q_2, "q_2/F");
    Run_Tree->Branch("l2_decayMode", &l2_decayMode, "l2_decayMode/F");
    Run_Tree->Branch("decayModeFinding_2", &decayModeFinding_2, "decayModeFinding_2/F");
    Run_Tree->Branch("byVVVLooseDeepVSjet_2", &byVVVLooseDeepVSjet_2, "byVVVLooseDeepVSjet_2/F");
    Run_Tree->Branch("byVVLooseDeepVSjet_2", &byVVLooseDeepVSjet_2, "byVVLooseDeepVSjet_2/F");
    Run_Tree->Branch("byVLooseDeepVSjet_2", &byVLooseDeepVSjet_2, "byVLooseDeepVSjet_2/F");
    Run_Tree->Branch("byLooseDeepVSjet_2", &byLooseDeepVSjet_2, "byLooseDeepVSjet_2/F");
    Run_Tree->Branch("byMediumDeepVSjet_2", &byMediumDeepVSjet_2, "byMediumDeepVSjet_2/F");
    Run_Tree->Branch("byTightDeepVSjet_2", &byTightDeepVSjet_2, "byTightDeepVSjet_2/F");
    Run_Tree->Branch("byVTightDeepVSjet_2", &byVTightDeepVSjet_2, "byVTightDeepVSjet_2/F");
    Run_Tree->Branch("byVVTightDeepVSjet_2", &byVVTightDeepVSjet_2, "byVVTightDeepVSjet_2/F");
    Run_Tree->Branch("byVLooseDeepVSmu_2", &byVLooseDeepVSmu_2, "byVLooseDeepVSmu_2/F");
    Run_Tree->Branch("byLooseDeepVSmu_2", &byLooseDeepVSmu_2, "byLooseDeepVSmu_2/F");
    Run_Tree->Branch("byMediumDeepVSmu_2", &byMediumDeepVSmu_2, "byMediumDeepVSmu_2/F");
    Run_Tree->Branch("byTightDeepVSmu_2", &byTightDeepVSmu_2, "byTightDeepVSmu_2/F");
    Run_Tree->Branch("byVTightDeepVSmu_2", &byVTightDeepVSmu_2, "byVTightDeepVSmu_2/F");
    Run_Tree->Branch("byVVTightDeepVSmu_2", &byVVTightDeepVSmu_2, "byVVTightDeepVSmu_2/F");
    Run_Tree->Branch("byVVVLooseDeepVSe_2", &byVVVLooseDeepVSe_2, "byVVVLooseDeepVSe_2/F");
    Run_Tree->Branch("byVVLooseDeepVSe_2", &byVVLooseDeepVSe_2, "byVVLooseDeepVSe_2/F");
    Run_Tree->Branch("byVLooseDeepVSe_2", &byVLooseDeepVSe_2, "byVLooseDeepVSe_2/F");
    Run_Tree->Branch("byLooseDeepVSe_2", &byLooseDeepVSe_2, "byLooseDeepVSe_2/F");
    Run_Tree->Branch("byMediumDeepVSe_2", &byMediumDeepVSe_2, "byMediumDeepVSe_2/F");
    Run_Tree->Branch("byTightDeepVSe_2", &byTightDeepVSe_2, "byTightDeepVSe_2/F");
    Run_Tree->Branch("byVTightDeepVSe_2", &byVTightDeepVSe_2, "byVTightDeepVSe_2/F");
    Run_Tree->Branch("byVVTightDeepVSe_2", &byVVTightDeepVSe_2, "byVVTightDeepVSe_2/F");

    Run_Tree->Branch("numGenJets", &numGenJets, "numGenJets/F");
    Run_Tree->Branch("jetPt_2", &jetPt_2, "jetPt_2/F");

    Run_Tree->Branch("Flag_ecalBadCalibReducedMINIAODFilter", &Flag_ecalBadCalibReducedMINIAODFilter, "Flag_ecalBadCalibReducedMINIAODFilter/F");
    Run_Tree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/F");
    Run_Tree->Branch("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, "Flag_globalSuperTightHalo2016Filter/F");
    Run_Tree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/F");
    Run_Tree->Branch("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, "Flag_ecalBadCalibFilter/F");
    Run_Tree->Branch("Flag_badMuons", &Flag_badMuons, "Flag_badMuons/F");
    Run_Tree->Branch("Flag_duplicateMuons", &Flag_duplicateMuons, "Flag_duplicateMuons/F");
    Run_Tree->Branch("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, "Flag_HBHENoiseIsoFilter/F");
    Run_Tree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/F");
    Run_Tree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/F");
    Run_Tree->Branch("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, "Flag_BadPFMuonFilter/F");
    Run_Tree->Branch("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, "Flag_BadChargedCandidateFilter/F");

    Run_Tree->Branch("met", &met, "met/F");
    Run_Tree->Branch("metSig", &metSig, "metSig/F");
    Run_Tree->Branch("metcov00", &metcov00, "metcov00/F");
    Run_Tree->Branch("metcov10", &metcov10, "metcov10/F");
    Run_Tree->Branch("metcov11", &metcov11, "metcov11/F");
    Run_Tree->Branch("metcov01", &metcov01, "metcov01/F");
    Run_Tree->Branch("metphi", &metphi, "metphi/F");
    Run_Tree->Branch("met_py", &met_py, "met_py/F");
    Run_Tree->Branch("met_px", &met_px, "met_px/F");
    Run_Tree->Branch("met_UESUp", &met_UESUp, "met_UESUp/F");
    Run_Tree->Branch("metphi_UESUp", &metphi_UESUp, "metphi_UESUp/F");
    Run_Tree->Branch("met_UESDown", &met_UESDown, "met_UESDown/F");
    Run_Tree->Branch("metphi_UESDown", &metphi_UESDown, "metphi_UESDown/F");

    Run_Tree->Branch("met_JetAbsoluteUp", &met_JetAbsoluteUp, "met_JetAbsoluteUp/F");
    Run_Tree->Branch("metphi_JetAbsoluteUp", &metphi_JetAbsoluteUp, "metphi_JetAbsoluteUp/F");
    Run_Tree->Branch("met_JetAbsoluteDown", &met_JetAbsoluteDown, "met_JetAbsoluteDown/F");
    Run_Tree->Branch("metphi_JetAbsoluteDown", &metphi_JetAbsoluteDown, "metphi_JetAbsoluteDown/F");
    Run_Tree->Branch("met_JetAbsoluteyearUp", &met_JetAbsoluteyearUp, "met_JetAbsoluteyearUp/F");
    Run_Tree->Branch("metphi_JetAbsoluteyearUp", &metphi_JetAbsoluteyearUp, "metphi_JetAbsoluteyearUp/F");
    Run_Tree->Branch("met_JetAbsoluteyearDown", &met_JetAbsoluteyearDown, "met_JetAbsoluteyearDown/F");
    Run_Tree->Branch("metphi_JetAbsoluteyearDown", &metphi_JetAbsoluteyearDown, "metphi_JetAbsoluteyearDown/F");
    Run_Tree->Branch("met_JetBBEC1Up", &met_JetBBEC1Up, "met_JetBBEC1Up/F");
    Run_Tree->Branch("metphi_JetBBEC1Up", &metphi_JetBBEC1Up, "metphi_JetBBEC1Up/F");
    Run_Tree->Branch("met_JetBBEC1Down", &met_JetBBEC1Down, "met_JetBBEC1Down/F");
    Run_Tree->Branch("metphi_JetBBEC1Down", &metphi_JetBBEC1Down, "metphi_JetBBEC1Down/F");
    Run_Tree->Branch("met_JetBBEC1yearUp", &met_JetBBEC1yearUp, "met_JetBBEC1yearUp/F");
    Run_Tree->Branch("metphi_JetBBEC1yearUp", &metphi_JetBBEC1yearUp, "metphi_JetBBEC1yearUp/F");
    Run_Tree->Branch("met_JetBBEC1yearDown", &met_JetBBEC1yearDown, "met_JetBBEC1yearDown/F");
    Run_Tree->Branch("metphi_JetBBEC1yearDown", &metphi_JetBBEC1yearDown, "metphi_JetBBEC1yearDown/F");
    Run_Tree->Branch("met_JetEC2Up", &met_JetEC2Up, "met_JetEC2Up/F");
    Run_Tree->Branch("metphi_JetEC2Up", &metphi_JetEC2Up, "metphi_JetEC2Up/F");
    Run_Tree->Branch("met_JetEC2Down", &met_JetEC2Down, "met_JetEC2Down/F");
    Run_Tree->Branch("metphi_JetEC2Down", &metphi_JetEC2Down, "metphi_JetEC2Down/F");
    Run_Tree->Branch("met_JetEC2yearUp", &met_JetEC2yearUp, "met_JetEC2yearUp/F");
    Run_Tree->Branch("metphi_JetEC2yearUp", &metphi_JetEC2yearUp, "metphi_JetEC2yearUp/F");
    Run_Tree->Branch("met_JetEC2yearDown", &met_JetEC2yearDown, "met_JetEC2yearDown/F");
    Run_Tree->Branch("metphi_JetEC2yearDown", &metphi_JetEC2yearDown, "metphi_JetEC2yearDown/F");
    Run_Tree->Branch("met_JetFlavorQCDUp", &met_JetFlavorQCDUp, "met_JetFlavorQCDUp/F");
    Run_Tree->Branch("metphi_JetFlavorQCDUp", &metphi_JetFlavorQCDUp, "metphi_JetFlavorQCDUp/F");
    Run_Tree->Branch("met_JetFlavorQCDDown", &met_JetFlavorQCDDown, "met_JetFlavorQCDDown/F");
    Run_Tree->Branch("metphi_JetFlavorQCDDown", &metphi_JetFlavorQCDDown, "metphi_JetFlavorQCDDown/F");
    Run_Tree->Branch("met_JetHFUp", &met_JetHFUp, "met_JetHFUp/F");
    Run_Tree->Branch("metphi_JetHFUp", &metphi_JetHFUp, "metphi_JetHFUp/F");
    Run_Tree->Branch("met_JetHFDown", &met_JetHFDown, "met_JetHFDown/F");
    Run_Tree->Branch("metphi_JetHFDown", &metphi_JetHFDown, "metphi_JetHFDown/F");
    Run_Tree->Branch("met_JetHFyearUp", &met_JetHFyearUp, "met_JetHFyearUp/F");
    Run_Tree->Branch("metphi_JetHFyearUp", &metphi_JetHFyearUp, "metphi_JetHFyearUp/F");
    Run_Tree->Branch("met_JetHFyearDown", &met_JetHFyearDown, "met_JetHFyearDown/F");
    Run_Tree->Branch("metphi_JetHFyearDown", &metphi_JetHFyearDown, "metphi_JetHFyearDown/F");
    Run_Tree->Branch("met_JetRelativeBalUp", &met_JetRelativeBalUp, "met_JetRelativeBalUp/F");
    Run_Tree->Branch("metphi_JetRelativeBalUp", &metphi_JetRelativeBalUp, "metphi_JetRelativeBalUp/F");
    Run_Tree->Branch("met_JetRelativeBalDown", &met_JetRelativeBalDown, "met_JetRelativeBalDown/F");
    Run_Tree->Branch("metphi_JetRelativeBalDown", &metphi_JetRelativeBalDown, "metphi_JetRelativeBalDown/F");
    Run_Tree->Branch("met_JetRelativeSampleUp", &met_JetRelativeSampleUp, "met_JetRelativeSampleUp/F");
    Run_Tree->Branch("metphi_JetRelativeSampleUp", &metphi_JetRelativeSampleUp, "metphi_JetRelativeSampleUp/F");
    Run_Tree->Branch("met_JetRelativeSampleDown", &met_JetRelativeSampleDown, "met_JetRelativeSampleDown/F");
    Run_Tree->Branch("metphi_JetRelativeSampleDown", &metphi_JetRelativeSampleDown, "metphi_JetRelativeSampleDown/F");
    Run_Tree->Branch("met_JERUp", &met_JERUp, "met_JERUp/F");
    Run_Tree->Branch("metphi_JERUp", &metphi_JERUp, "metphi_JERUp/F");
    Run_Tree->Branch("met_JERDown", &met_JERDown, "met_JERDown/F");
    Run_Tree->Branch("metphi_JERDown", &metphi_JERDown, "metphi_JERDown/F");

    Run_Tree->Branch("met_responseUp", &met_responseUp, "met_responseUp/F");
    Run_Tree->Branch("met_responseDown", &met_responseDown, "met_responseDown/F");
    Run_Tree->Branch("met_resolutionUp", &met_resolutionUp, "met_resolutionUp/F");
    Run_Tree->Branch("met_resolutionDown", &met_resolutionDown, "met_resolutionDown/F");
    Run_Tree->Branch("metphi_responseUp", &metphi_responseUp, "metphi_responseUp/F");
    Run_Tree->Branch("metphi_responseDown", &metphi_responseDown, "metphi_responseDown/F");
    Run_Tree->Branch("metphi_resolutionUp", &metphi_resolutionUp, "metphi_resolutionUp/F");
    Run_Tree->Branch("metphi_resolutionDown", &metphi_resolutionDown, "metphi_resolutionDown/F");

    Run_Tree->Branch("passMu24", &passMu24, "passMu24/F");
    Run_Tree->Branch("passMu27", &passMu27, "passMu27/F");
    Run_Tree->Branch("passMu20Tau27", &passMu20Tau27, "passMu20Tau27/F");
    Run_Tree->Branch("passMu20HPSTau27", &passMu20HPSTau27, "passMu20HPSTau27/F");

    Run_Tree->Branch("matchMu24_1", &matchMu24_1, "matchMu24_1/F");
    Run_Tree->Branch("matchMu27_1", &matchMu27_1, "matchMu27_1/F");
    Run_Tree->Branch("matchMu20Tau27_1", &matchMu20Tau27_1, "matchMu20Tau27_1/F");
    Run_Tree->Branch("matchMu20Tau27_2", &matchMu20Tau27_2, "matchMu20Tau27_2/F");
    Run_Tree->Branch("matchMu20HPSTau27_1", &matchMu20HPSTau27_1, "matchMu20HPSTau27_1/F");
    Run_Tree->Branch("matchMu20HPSTau27_2", &matchMu20HPSTau27_2, "matchMu20HPSTau27_2/F");
    Run_Tree->Branch("filterMu24_1", &filterMu24_1, "filterMu24_1/F");
    Run_Tree->Branch("filterMu27_1", &filterMu27_1, "filterMu27_1/F");
    Run_Tree->Branch("filterMu20Tau27_1", &filterMu20Tau27_1, "filterMu20Tau27_1/F");
    Run_Tree->Branch("filterMu20Tau27_2", &filterMu20Tau27_2, "filterMu20Tau27_2/F");
    Run_Tree->Branch("filterMu20HPSTau27_1", &filterMu20HPSTau27_1, "filterMu20HPSTau27_1/F");
    Run_Tree->Branch("filterMu20HPSTau27_2", &filterMu20HPSTau27_2, "filterMu20HPSTau27_2/F");

    Run_Tree->Branch("mjj", &mjj, "mjj/F");
    Run_Tree->Branch("mjj_JetAbsoluteUp", &mjj_JetAbsoluteUp, "mjj_JetAbsoluteUp/F");
    Run_Tree->Branch("mjj_JetAbsoluteDown", &mjj_JetAbsoluteDown, "mjj_JetAbsoluteDown/F");
    Run_Tree->Branch("mjj_JetAbsoluteyearUp", &mjj_JetAbsoluteyearUp, "mjj_JetAbsoluteyearUp/F");
    Run_Tree->Branch("mjj_JetAbsoluteyearDown", &mjj_JetAbsoluteyearDown, "mjj_JetAbsoluteyearDown/F");
    Run_Tree->Branch("mjj_JetBBEC1Up", &mjj_JetBBEC1Up, "mjj_JetBBEC1Up/F");
    Run_Tree->Branch("mjj_JetBBEC1Down", &mjj_JetBBEC1Down, "mjj_JetBBEC1Down/F");
    Run_Tree->Branch("mjj_JetBBEC1yearUp", &mjj_JetBBEC1yearUp, "mjj_JetBBEC1yearUp/F");
    Run_Tree->Branch("mjj_JetBBEC1yearDown", &mjj_JetBBEC1yearDown, "mjj_JetBBEC1yearDown/F");
    Run_Tree->Branch("mjj_JetEC2Up", &mjj_JetEC2Up, "mjj_JetEC2Up/F");
    Run_Tree->Branch("mjj_JetEC2Down", &mjj_JetEC2Down, "mjj_JetEC2Down/F");
    Run_Tree->Branch("mjj_JetEC2yearUp", &mjj_JetEC2yearUp, "mjj_JetEC2yearUp/F");
    Run_Tree->Branch("mjj_JetEC2yearDown", &mjj_JetEC2yearDown, "mjj_JetEC2yearDown/F");
    Run_Tree->Branch("mjj_JetFlavorQCDUp", &mjj_JetFlavorQCDUp, "mjj_JetFlavorQCDUp/F");
    Run_Tree->Branch("mjj_JetFlavorQCDDown", &mjj_JetFlavorQCDDown, "mjj_JetFlavorQCDDown/F");
    Run_Tree->Branch("mjj_JetHFUp", &mjj_JetHFUp, "mjj_JetHFUp/F");
    Run_Tree->Branch("mjj_JetHFDown", &mjj_JetHFDown, "mjj_JetHFDown/F");
    Run_Tree->Branch("mjj_JetHFyearUp", &mjj_JetHFyearUp, "mjj_JetHFyearUp/F");
    Run_Tree->Branch("mjj_JetHFyearDown", &mjj_JetHFyearDown, "mjj_JetHFyearDown/F");
    Run_Tree->Branch("mjj_JetRelativeBalUp", &mjj_JetRelativeBalUp, "mjj_JetRelativeBalUp/F");
    Run_Tree->Branch("mjj_JetRelativeBalDown", &mjj_JetRelativeBalDown, "mjj_JetRelativeBalDown/F");
    Run_Tree->Branch("mjj_JetRelativeSampleUp", &mjj_JetRelativeSampleUp, "mjj_JetRelativeSampleUp/F");
    Run_Tree->Branch("mjj_JetRelativeSampleDown", &mjj_JetRelativeSampleDown, "mjj_JetRelativeSampleDown/F");
    Run_Tree->Branch("mjj_JERUp", &mjj_JERUp, "mjj_JERUp/F");
    Run_Tree->Branch("mjj_JERDown", &mjj_JERDown, "mjj_JERDown/F");

    Run_Tree->Branch("gen_match_1", &gen_match_1, "gen_match_1/I");
    Run_Tree->Branch("gen_match_2", &gen_match_2, "gen_match_2/I");

    Run_Tree->Branch("bpt_1", &bpt_1, "bpt_1/F");
    Run_Tree->Branch("beta_1", &beta_1, "beta_1/F");
    Run_Tree->Branch("bphi_1", &bphi_1, "bphi_1/F");
    Run_Tree->Branch("bflavor_1", &bflavor_1, "bflavor_1/F");
    Run_Tree->Branch("bcsv_1", &bcsv_1, "bcsv_1/F");

    Run_Tree->Branch("bpt_2", &bpt_2, "bpt_2/F");
    Run_Tree->Branch("beta_2", &beta_2, "beta_2/F");
    Run_Tree->Branch("bphi_2", &bphi_2, "bphi_2/F");
    Run_Tree->Branch("bflavor_2", &bflavor_2, "bflavor_2/F");
    Run_Tree->Branch("bcsv_2", &bcsv_2, "bcsv_2/F");

    Run_Tree->Branch("nbtag", &nbtag, "nbtag/I");
    Run_Tree->Branch("nbtagL", &nbtagL, "nbtagL/I");
    Run_Tree->Branch("njets", &njets, "njets/I");
    Run_Tree->Branch("njets_JetAbsoluteUp", &njets_JetAbsoluteUp, "njets_JetAbsoluteUp/I");
    Run_Tree->Branch("njets_JetAbsoluteDown", &njets_JetAbsoluteDown, "njets_JetAbsoluteDown/I");
    Run_Tree->Branch("njets_JetAbsoluteyearUp", &njets_JetAbsoluteyearUp, "njets_JetAbsoluteyearUp/I");
    Run_Tree->Branch("njets_JetAbsoluteyearDown", &njets_JetAbsoluteyearDown, "njets_JetAbsoluteyearDown/I");
    Run_Tree->Branch("njets_JetBBEC1Up", &njets_JetBBEC1Up, "njets_JetBBEC1Up/I");
    Run_Tree->Branch("njets_JetBBEC1Down", &njets_JetBBEC1Down, "njets_JetBBEC1Down/I");
    Run_Tree->Branch("njets_JetBBEC1yearUp", &njets_JetBBEC1yearUp, "njets_JetBBEC1yearUp/I");
    Run_Tree->Branch("njets_JetBBEC1yearDown", &njets_JetBBEC1yearDown, "njets_JetBBEC1yearDown/I");
    Run_Tree->Branch("njets_JetEC2Up", &njets_JetEC2Up, "njets_JetEC2Up/I");
    Run_Tree->Branch("njets_JetEC2Down", &njets_JetEC2Down, "njets_JetEC2Down/I");
    Run_Tree->Branch("njets_JetEC2yearUp", &njets_JetEC2yearUp, "njets_JetEC2yearUp/I");
    Run_Tree->Branch("njets_JetEC2yearDown", &njets_JetEC2yearDown, "njets_JetEC2yearDown/I");
    Run_Tree->Branch("njets_JetFlavorQCDUp", &njets_JetFlavorQCDUp, "njets_JetFlavorQCDUp/I");
    Run_Tree->Branch("njets_JetFlavorQCDDown", &njets_JetFlavorQCDDown, "njets_JetFlavorQCDDown/I");
    Run_Tree->Branch("njets_JetHFUp", &njets_JetHFUp, "njets_JetHFUp/I");
    Run_Tree->Branch("njets_JetHFDown", &njets_JetHFDown, "njets_JetHFDown/I");
    Run_Tree->Branch("njets_JetHFyearUp", &njets_JetHFyearUp, "njets_JetHFyearUp/I");
    Run_Tree->Branch("njets_JetHFyearDown", &njets_JetHFyearDown, "njets_JetHFyearDown/I");
    Run_Tree->Branch("njets_JetRelativeBalUp", &njets_JetRelativeBalUp, "njets_JetRelativeBalUp/I");
    Run_Tree->Branch("njets_JetRelativeBalDown", &njets_JetRelativeBalDown, "njets_JetRelativeBalDown/I");
    Run_Tree->Branch("njets_JetRelativeSampleUp", &njets_JetRelativeSampleUp, "njets_JetRelativeSampleUp/I");
    Run_Tree->Branch("njets_JetRelativeSampleDown", &njets_JetRelativeSampleDown, "njets_JetRelativeSampleDown/I");
    Run_Tree->Branch("njets_JERUp", &njets_JERUp, "njets_JERUp/I");
    Run_Tree->Branch("njets_JERDown", &njets_JERDown, "njets_JERDown/I");

    Run_Tree->Branch("jpt_1", &jpt_1, "jpt_1/F");
    Run_Tree->Branch("jpt_1_JetAbsoluteUp", &jpt_1_JetAbsoluteUp, "jpt_1_JetAbsoluteUp/F");
    Run_Tree->Branch("jpt_1_JetAbsoluteDown", &jpt_1_JetAbsoluteDown, "jpt_1_JetAbsoluteDown/F");
    Run_Tree->Branch("jpt_1_JetAbsoluteyearUp", &jpt_1_JetAbsoluteyearUp, "jpt_1_JetAbsoluteyearUp/F");
    Run_Tree->Branch("jpt_1_JetAbsoluteyearDown", &jpt_1_JetAbsoluteyearDown, "jpt_1_JetAbsoluteyearDown/F");
    Run_Tree->Branch("jpt_1_JetBBEC1Up", &jpt_1_JetBBEC1Up, "jpt_1_JetBBEC1Up/F");
    Run_Tree->Branch("jpt_1_JetBBEC1Down", &jpt_1_JetBBEC1Down, "jpt_1_JetBBEC1Down/F");
    Run_Tree->Branch("jpt_1_JetBBEC1yearUp", &jpt_1_JetBBEC1yearUp, "jpt_1_JetBBEC1yearUp/F");
    Run_Tree->Branch("jpt_1_JetBBEC1yearDown", &jpt_1_JetBBEC1yearDown, "jpt_1_JetBBEC1yearDown/F");
    Run_Tree->Branch("jpt_1_JetEC2Up", &jpt_1_JetEC2Up, "jpt_1_JetEC2Up/F");
    Run_Tree->Branch("jpt_1_JetEC2Down", &jpt_1_JetEC2Down, "jpt_1_JetEC2Down/F");
    Run_Tree->Branch("jpt_1_JetEC2yearUp", &jpt_1_JetEC2yearUp, "jpt_1_JetEC2yearUp/F");
    Run_Tree->Branch("jpt_1_JetEC2yearDown", &jpt_1_JetEC2yearDown, "jpt_1_JetEC2yearDown/F");
    Run_Tree->Branch("jpt_1_JetFlavorQCDUp", &jpt_1_JetFlavorQCDUp, "jpt_1_JetFlavorQCDUp/F");
    Run_Tree->Branch("jpt_1_JetFlavorQCDDown", &jpt_1_JetFlavorQCDDown, "jpt_1_JetFlavorQCDDown/F");
    Run_Tree->Branch("jpt_1_JetHFUp", &jpt_1_JetHFUp, "jpt_1_JetHFUp/F");
    Run_Tree->Branch("jpt_1_JetHFDown", &jpt_1_JetHFDown, "jpt_1_JetHFDown/F");
    Run_Tree->Branch("jpt_1_JetHFyearUp", &jpt_1_JetHFyearUp, "jpt_1_JetHFyearUp/F");
    Run_Tree->Branch("jpt_1_JetHFyearDown", &jpt_1_JetHFyearDown, "jpt_1_JetHFyearDown/F");
    Run_Tree->Branch("jpt_1_JetRelativeBalUp", &jpt_1_JetRelativeBalUp, "jpt_1_JetRelativeBalUp/F");
    Run_Tree->Branch("jpt_1_JetRelativeBalDown", &jpt_1_JetRelativeBalDown, "jpt_1_JetRelativeBalDown/F");
    Run_Tree->Branch("jpt_1_JetRelativeSampleUp", &jpt_1_JetRelativeSampleUp, "jpt_1_JetRelativeSampleUp/F");
    Run_Tree->Branch("jpt_1_JetRelativeSampleDown", &jpt_1_JetRelativeSampleDown, "jpt_1_JetRelativeSampleDown/F");
    Run_Tree->Branch("jpt_1_JERUp", &jpt_1_JERUp, "jpt_1_JERUp/F");
    Run_Tree->Branch("jpt_1_JERDown", &jpt_1_JERDown, "jpt_1_JERDown/F");

    Run_Tree->Branch("jeta_1", &jeta_1, "jeta_1/F");
    Run_Tree->Branch("jphi_1", &jphi_1, "jphi_1/F");

    Run_Tree->Branch("jpt_2", &jpt_2, "jpt_2/F");
    Run_Tree->Branch("jeta_2", &jeta_2, "jeta_2/F");
    Run_Tree->Branch("jphi_2", &jphi_2, "jphi_2/F");

    Run_Tree->Branch("pt_top1", &pt_top1, "pt_top1/F");
    Run_Tree->Branch("pt_top2", &pt_top2, "pt_top2/F");
    Run_Tree->Branch("genweight", &genweight, "genweight/F");

    Run_Tree->Branch("gen_Higgs_pt", &gen_Higgs_pt, "gen_Higgs_pt/F");
    Run_Tree->Branch("gen_Higgs_mass", &gen_Higgs_mass, "gen_Higgs_mass/F");

    Run_Tree->Branch("genpt_1", &genpt_1, "genpt_1/F");
    Run_Tree->Branch("geneta_1", &geneta_1, "geneta_1/F");
    Run_Tree->Branch("genpt_2", &genpt_2, "genpt_2/F");
    Run_Tree->Branch("geneta_2", &geneta_2, "geneta_2/F");

    int bestEntry=-1;
    ULong64_t evt_now=0;
    ULong64_t evt_before=-1;
    float muiso_before=100;
    float tauiso_before=-1000;
    float mupt_before=0;
    float taupt_before=0;
    plotFill("nevents",0,9,0,9,evCounter->GetBinContent(1));
    plotFill("nevents",1,9,0,9,evCounterW->GetBinContent(1));
    for (int iEntry = 0; iEntry < tree->GetEntries() ; iEntry++)
    {
        float pu=1.0;
        tree->GetEntry(iEntry);
        bool print=false;
        if (iEntry % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d ", iEntry);
        fflush(stdout);
        plotFill("pileup_mc",tree->nTruePU,80,0,80);
        TLorentzVector dau1;
        TLorentzVector dau2;
        dau1.SetPtEtaPhiM(tree->mPt,tree->mEta,tree->mPhi,tree->mMass);
        dau2.SetPtEtaPhiM(tree->tPt,tree->tEta,tree->tPhi,tree->tMass);
        if (isMC && tree->tZTTGenMatching==5 && tree->tDecayMode==0) dau2=dau2*0.987;
        else if (isMC && tree->tZTTGenMatching==5 && tree->tDecayMode==1) dau2=dau2*0.995;
        else if (isMC && tree->tZTTGenMatching==5 && tree->tDecayMode==10) dau2=dau2*0.988;
        else if (isMC && tree->tZTTGenMatching==5 && tree->tDecayMode==11) dau2=dau2*0.988;

        if (isMC && (tree->tZTTGenMatching==1 or tree->tZTTGenMatching==3) && tree->tDecayMode==0 && fabs(dau2.Eta())<1.5) dau2=dau2*1.01362;
        else if (isMC && (tree->tZTTGenMatching==1 or tree->tZTTGenMatching==3) && tree->tDecayMode==1 && fabs(dau2.Eta())<1.5) dau2=dau2*1.01945;
        else if (isMC && (tree->tZTTGenMatching==1 or tree->tZTTGenMatching==3) && tree->tDecayMode==0 && fabs(dau2.Eta())>1.5) dau2=dau2*0.96903;
        else if (isMC && (tree->tZTTGenMatching==1 or tree->tZTTGenMatching==3) && tree->tDecayMode==1 && fabs(dau2.Eta())>1.5) dau2=dau2*0.985;

        /*if (isMC && (tree->tZTTGenMatching==1 or tree->tZTTGenMatching==3) && tree->tDecayMode==0) dau2=dau2*0.968;
        else if (isMC && (tree->tZTTGenMatching==1 or tree->tZTTGenMatching==3) && tree->tDecayMode==1) dau2=dau2*1.026;
        if (isMC && (tree->tZTTGenMatching==2 or tree->tZTTGenMatching==4) && tree->tDecayMode==0) dau2=dau2*0.998;
        else if (isMC && (tree->tZTTGenMatching==2 or tree->tZTTGenMatching==4) && tree->tDecayMode==1) dau2=dau2*0.990;*/

        if (dau1.DeltaR(dau2)<0.5) continue;
	if (fabs(tree->mPVDXY)>0.045) continue;
        if (fabs(tree->mPVDZ)>0.2) continue;
        if (fabs(tree->tPVDZ)>0.2) continue;
	if (dau1.Pt()<19.5 or dau2.Pt()<29.5) continue;
        if (fabs(dau1.Eta())>2.1 or fabs(dau2.Eta())>2.3) continue;
        if (!tree->tVVVLooseDeepTau2017v2p1VSjet) continue;
        if (!tree->tVLooseDeepTau2017v2p1VSmu) continue;
        if (!tree->tVVVLooseDeepTau2017v2p1VSe) continue;
        if (tree->tDecayMode==5 or tree->tDecayMode==6) continue;
        if (tree->mRelPFIsoDBDefault>0.5) continue;
	if (!tree->mPFIDMedium) continue;
	if (tree->eVetoZTTp001dxyzR0>0) continue;
	if (tree->muVetoZTTp001dxyzR0>1) continue;
	if (tree->dimuonVeto>0) continue;
	evt_now=tree->evt;
	if (evt_now!=evt_before){
	   mupt_before=tree->mPt;
	   muiso_before=tree->mRelPFIsoDBDefault;
	   taupt_before=tree->tPt;
	   tauiso_before=tree->tDeepTau2017v2p1VSjetraw;
	}
        if (evt_now!=evt_before){
           if (bestEntry>-1){
              fillTree(Run_Tree,tree,bestEntry,recoil,isMC,isEmbedded);
           }
	   bestEntry=iEntry;
	}
	if (evt_now==evt_before){
	   if (tree->mRelPFIsoDBDefault<muiso_before or (tree->mRelPFIsoDBDefault==muiso_before && tree->mPt>mupt_before) or (tree->mRelPFIsoDBDefault==muiso_before && tree->mPt==mupt_before && tree->tDeepTau2017v2p1VSjetraw>tauiso_before) or (tree->mRelPFIsoDBDefault==muiso_before && tree->mPt==mupt_before && tree->tDeepTau2017v2p1VSjetraw==tauiso_before && tree->tPt>taupt_before) ){
		bestEntry=iEntry;
	        muiso_before=tree->mRelPFIsoDBDefault;
		mupt_before=tree->mPt;
		tauiso_before=tree->tDeepTau2017v2p1VSjetraw;
		taupt_before=tree->tPt;
	   }
	}
	evt_before=evt_now;
    }
    if (bestEntry>-1)
       fillTree(Run_Tree,tree,bestEntry,recoil,isMC,isEmbedded);
    fout->cd();
    Run_Tree->Write();
    map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
    map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
    for (; iMap1 != jMap1; ++iMap1)
        nplot1(iMap1->first)->Write();
    map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
    map<string, TH2F*>::const_iterator jMap2 = myMap2->end();
    for (; iMap2 != jMap2; ++iMap2)
        nplot2(iMap2->first)->Write();

    fout->Close();
    return 0;
}

