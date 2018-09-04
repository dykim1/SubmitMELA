#include "PhysicsTools/FWLite/interface/CommandLineParser.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "TSystem.h"
#include "TMath.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"

#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "ZZMatrixElement/MELA/interface/TUtil.hh"

using namespace std;

Mela mela(13.0, 125.6, TVar::ERROR);

void calculateME(TLorentzVector tau1, TLorentzVector tau2, TLorentzVector jet1, TLorentzVector jet2,
		 int tauCharge1, int tauCharge2,
		 float& ME_sm, float& ME_sm_ggH, float& ME_sm_WH, float& ME_sm_ZH,
		 float& ME_bkg1, float& ME_bkg2,
		 float& Q2V1, float& Q2V2, float& costheta1, float& costheta2, float& Phi, float& costhetastar, float& Phi1);

void processFile(TDirectory*  dir, optutl::CommandLineParser parser, char treeToUse[], bool trueTau = true, int doES = 0);

// files from Tyler R. to copy directory/root file content exactly
void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);

int main (int argc, char* argv[]) {

  optutl::CommandLineParser parser("Input parameters");
  parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile.root");
  parser.addOption("inputFile",optutl::CommandLineParser::kString,"input File");
  parser.addOption("trueTau",optutl::CommandLineParser::kBool,"use true 4-vectors of tau leptons",true);
  parser.addOption("doES",optutl::CommandLineParser::kDouble,"doES",0.0);
  parser.parseArguments (argc, argv);
  std::cout << "Input commands:" 
	    << "\n -- input file: " << parser.stringValue("inputFile")
	    << "\n -- output file: " << parser.stringValue("newFile")
	    << "\n -- use true tau 4-vectors (1 - yes, 0 - no): " << parser.boolValue("trueTau")
            << "\n -- doES: " << parser.doubleValue("doES") << std::endl;
  char treeToUse[80] = "first";

  TFile *fProduce;
  TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"READ");
  std::cout<<"Creating new outputfile"<<std::endl;
  std::string newFileName = parser.stringValue("newFile");

  fProduce = new TFile(newFileName.c_str(),"RECREATE");
  copyFiles(parser, f, fProduce);//new TFile(parser.stringValue("inputFile").c_str()+"SVFit","UPDATE");
  fProduce = new TFile(newFileName.c_str(),"UPDATE");
  std::cout<<"listing the directories================="<<std::endl;
  fProduce->ls();

  processFile(fProduce, parser,  treeToUse, parser.boolValue("trueTau"), parser.doubleValue("doES"));
  f->Close();

  return 0;
}


void processFile(TDirectory*  dir, optutl::CommandLineParser parser, char treeToUse[], bool trueTau, int doES)
{
  std::cout << "Starting!!!" << std::endl;
  std::cout << "treeToUse: " << treeToUse << std::endl;
  char stringA[80]="first";

  TDirectory *dirsav = gDirectory;
  TIter next(dir->GetListOfKeys());
  TKey *key;
  std::string channel = "x";
  dir->cd();
  std::vector<TString> processedNames;

  while ((key = (TKey*)next())) {
    printf("Found key=%s \n",key->GetName());

    TObject *obj = key->ReadObj();
    if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
      std::cout << "This is a directory, diving in!" << std::endl;
      // zero the processedNames vector, to allow processing trees with duplicate names in separate directories
      processedNames.clear();

      dir->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      sprintf(treeToUse,"%s",key->GetName());
      processFile(subdir, parser, treeToUse, parser.boolValue("trueTau"),parser.doubleValue("doES"));

      dirsav->cd();
    } else if(obj->IsA()->InheritsFrom(TTree::Class())) {
      // check  if this tree was already processed
      std::vector<TString>::const_iterator it = find(processedNames.begin(), processedNames.end(), key->GetName());
      if ( it != processedNames.end() ) {
	std::cout << "This tree was already processed, skipping..." <<  std::endl;
        continue;
      }
      std::cout << "This is the tree! Start processing" << std::endl;
      processedNames.push_back(key->GetName());
      
      // Identify the process
      if ( std::string(key->GetName()).find("tt") != std::string::npos )  {
        channel = "tt";
	std::cout << "Identified channel tt" << std::endl;
      } 
      else if ( std::string(key->GetName()).find("em") != std::string::npos )  {
	channel = "em";
	std::cout << "Identified channel em" << std::endl;
      }
      else if ( std::string(key->GetName()).find("et") != std::string::npos ) {
        channel = "et";
	std::cout << "Identified channel et" << std::endl;
      } 
      else if ( std::string(key->GetName()).find("mt") != std::string::npos ) {
        channel = "mt";
	std::cout << "Identified channel mt" << std::endl;
      } else {
	std::cout<<"Tree "<< key->GetName() <<" does not match ... Skipping!!"<<std::endl;
        return;
      }

      TTree* tree = (TTree*)obj;

      float ME_sm_VBF, ME_sm_ggH, ME_sm_WH, ME_sm_ZH, ME_bkg, ME_bkg1, ME_bkg2;
      float mela_Dbkg_VBF;
      float mela_Dbkg_ggH;
      float mela_Dbkg_WH;
      float mela_Dbkg_ZH;
      
      // angles
      float Q2V1;
      float Q2V2;
      float costheta1;
      float costheta2;
      float Phi;
      float costhetastar;
      float Phi1;

      // Uncertainties
      float ME_sm_VBF_UP, ME_sm_ggH_UP, ME_sm_WH_UP, ME_sm_ZH_UP, ME_bkg_UP, ME_bkg1_UP, ME_bkg2_UP;
      float mela_Dbkg_VBF_UP;
      float mela_Dbkg_ggH_UP;
      float mela_Dbkg_WH_UP;
      float mela_Dbkg_ZH_UP;      
      float Q2V1_UP;
      float Q2V2_UP;
      float costheta1_UP;
      float costheta2_UP;
      float Phi_UP;
      float costhetastar_UP;
      float Phi1_UP;

      float ME_sm_VBF_DM0_UP, ME_sm_ggH_DM0_UP, ME_sm_WH_DM0_UP, ME_sm_ZH_DM0_UP, ME_bkg_DM0_UP, ME_bkg1_DM0_UP, ME_bkg2_DM0_UP;
      float mela_Dbkg_VBF_DM0_UP;
      float mela_Dbkg_ggH_DM0_UP;
      float mela_Dbkg_WH_DM0_UP;
      float mela_Dbkg_ZH_DM0_UP;      
      float Q2V1_DM0_UP;
      float Q2V2_DM0_UP;
      float costheta1_DM0_UP;
      float costheta2_DM0_UP;
      float Phi_DM0_UP;
      float costhetastar_DM0_UP;
      float Phi1_DM0_UP;

      float ME_sm_VBF_DM1_UP, ME_sm_ggH_DM1_UP, ME_sm_WH_DM1_UP, ME_sm_ZH_DM1_UP, ME_bkg_DM1_UP, ME_bkg1_DM1_UP, ME_bkg2_DM1_UP;
      float mela_Dbkg_VBF_DM1_UP;
      float mela_Dbkg_ggH_DM1_UP;
      float mela_Dbkg_WH_DM1_UP;
      float mela_Dbkg_ZH_DM1_UP;      
      float Q2V1_DM1_UP;
      float Q2V2_DM1_UP;
      float costheta1_DM1_UP;
      float costheta2_DM1_UP;
      float Phi_DM1_UP;
      float costhetastar_DM1_UP;
      float Phi1_DM1_UP;

      float ME_sm_VBF_DM10_UP, ME_sm_ggH_DM10_UP, ME_sm_WH_DM10_UP, ME_sm_ZH_DM10_UP, ME_bkg_DM10_UP, ME_bkg1_DM10_UP, ME_bkg2_DM10_UP;
      float mela_Dbkg_VBF_DM10_UP;
      float mela_Dbkg_ggH_DM10_UP;
      float mela_Dbkg_WH_DM10_UP;
      float mela_Dbkg_ZH_DM10_UP;      
      float Q2V1_DM10_UP;
      float Q2V2_DM10_UP;
      float costheta1_DM10_UP;
      float costheta2_DM10_UP;
      float Phi_DM10_UP;
      float costhetastar_DM10_UP;
      float Phi1_DM10_UP;

      float ME_sm_VBF_DOWN, ME_sm_ggH_DOWN, ME_sm_WH_DOWN, ME_sm_ZH_DOWN, ME_bkg_DOWN, ME_bkg1_DOWN, ME_bkg2_DOWN;
      float mela_Dbkg_VBF_DOWN;
      float mela_Dbkg_ggH_DOWN;
      float mela_Dbkg_WH_DOWN;
      float mela_Dbkg_ZH_DOWN;      
      float Q2V1_DOWN;
      float Q2V2_DOWN;
      float costheta1_DOWN;
      float costheta2_DOWN;
      float Phi_DOWN;
      float costhetastar_DOWN;
      float Phi1_DOWN;

      float ME_sm_VBF_DM0_DOWN, ME_sm_ggH_DM0_DOWN, ME_sm_WH_DM0_DOWN, ME_sm_ZH_DM0_DOWN, ME_bkg_DM0_DOWN, ME_bkg1_DM0_DOWN, ME_bkg2_DM0_DOWN;
      float mela_Dbkg_VBF_DM0_DOWN;
      float mela_Dbkg_ggH_DM0_DOWN;
      float mela_Dbkg_WH_DM0_DOWN;
      float mela_Dbkg_ZH_DM0_DOWN;      
      float Q2V1_DM0_DOWN;
      float Q2V2_DM0_DOWN;
      float costheta1_DM0_DOWN;
      float costheta2_DM0_DOWN;
      float Phi_DM0_DOWN;
      float costhetastar_DM0_DOWN;
      float Phi1_DM0_DOWN;

      float ME_sm_VBF_DM1_DOWN, ME_sm_ggH_DM1_DOWN, ME_sm_WH_DM1_DOWN, ME_sm_ZH_DM1_DOWN, ME_bkg_DM1_DOWN, ME_bkg1_DM1_DOWN, ME_bkg2_DM1_DOWN;
      float mela_Dbkg_VBF_DM1_DOWN;
      float mela_Dbkg_ggH_DM1_DOWN;
      float mela_Dbkg_WH_DM1_DOWN;
      float mela_Dbkg_ZH_DM1_DOWN;      
      float Q2V1_DM1_DOWN;
      float Q2V2_DM1_DOWN;
      float costheta1_DM1_DOWN;
      float costheta2_DM1_DOWN;
      float Phi_DM1_DOWN;
      float costhetastar_DM1_DOWN;
      float Phi1_DM1_DOWN;

      float ME_sm_VBF_DM10_DOWN, ME_sm_ggH_DM10_DOWN, ME_sm_WH_DM10_DOWN, ME_sm_ZH_DM10_DOWN, ME_bkg_DM10_DOWN, ME_bkg1_DM10_DOWN, ME_bkg2_DM10_DOWN;
      float mela_Dbkg_VBF_DM10_DOWN;
      float mela_Dbkg_ggH_DM10_DOWN;
      float mela_Dbkg_WH_DM10_DOWN;
      float mela_Dbkg_ZH_DM10_DOWN;      
      float Q2V1_DM10_DOWN;
      float Q2V2_DM10_DOWN;
      float costheta1_DM10_DOWN;
      float costheta2_DM10_DOWN;
      float Phi_DM10_DOWN;
      float costhetastar_DM10_DOWN;
      float Phi1_DM10_DOWN;

      float ME_sm_VBF_UncMet_UP, ME_sm_ggH_UncMet_UP, ME_sm_WH_UncMet_UP, ME_sm_ZH_UncMet_UP, ME_bkg_UncMet_UP, ME_bkg1_UncMet_UP, ME_bkg2_UncMet_UP;
      float mela_Dbkg_VBF_UncMet_UP;
      float mela_Dbkg_ggH_UncMet_UP;
      float mela_Dbkg_WH_UncMet_UP;
      float mela_Dbkg_ZH_UncMet_UP;      
      float Q2V1_UncMet_UP;
      float Q2V2_UncMet_UP;
      float costheta1_UncMet_UP;
      float costheta2_UncMet_UP;
      float Phi_UncMet_UP;
      float costhetastar_UncMet_UP;
      float Phi1_UncMet_UP;

      float ME_sm_VBF_UncMet_DOWN, ME_sm_ggH_UncMet_DOWN, ME_sm_WH_UncMet_DOWN, ME_sm_ZH_UncMet_DOWN, ME_bkg_UncMet_DOWN, ME_bkg1_UncMet_DOWN, ME_bkg2_UncMet_DOWN;
      float mela_Dbkg_VBF_UncMet_DOWN;
      float mela_Dbkg_ggH_UncMet_DOWN;
      float mela_Dbkg_WH_UncMet_DOWN;
      float mela_Dbkg_ZH_UncMet_DOWN;      
      float Q2V1_UncMet_DOWN;
      float Q2V2_UncMet_DOWN;
      float costheta1_UncMet_DOWN;
      float costheta2_UncMet_DOWN;
      float Phi_UncMet_DOWN;
      float costhetastar_UncMet_DOWN;
      float Phi1_UncMet_DOWN;

      float ME_sm_VBF_ClusteredMet_UP, ME_sm_ggH_ClusteredMet_UP, ME_sm_WH_ClusteredMet_UP, ME_sm_ZH_ClusteredMet_UP, ME_bkg_ClusteredMet_UP, ME_bkg1_ClusteredMet_UP, ME_bkg2_ClusteredMet_UP;
      float mela_Dbkg_VBF_ClusteredMet_UP;
      float mela_Dbkg_ggH_ClusteredMet_UP;
      float mela_Dbkg_WH_ClusteredMet_UP;
      float mela_Dbkg_ZH_ClusteredMet_UP;      
      float Q2V1_ClusteredMet_UP;
      float Q2V2_ClusteredMet_UP;
      float costheta1_ClusteredMet_UP;
      float costheta2_ClusteredMet_UP;
      float Phi_ClusteredMet_UP;
      float costhetastar_ClusteredMet_UP;
      float Phi1_ClusteredMet_UP;

      float ME_sm_VBF_ClusteredMet_DOWN, ME_sm_ggH_ClusteredMet_DOWN, ME_sm_WH_ClusteredMet_DOWN, ME_sm_ZH_ClusteredMet_DOWN, ME_bkg_ClusteredMet_DOWN, ME_bkg1_ClusteredMet_DOWN, ME_bkg2_ClusteredMet_DOWN;
      float mela_Dbkg_VBF_ClusteredMet_DOWN;
      float mela_Dbkg_ggH_ClusteredMet_DOWN;
      float mela_Dbkg_WH_ClusteredMet_DOWN;
      float mela_Dbkg_ZH_ClusteredMet_DOWN;      
      float Q2V1_ClusteredMet_DOWN;
      float Q2V2_ClusteredMet_DOWN;
      float costheta1_ClusteredMet_DOWN;
      float costheta2_ClusteredMet_DOWN;
      float Phi_ClusteredMet_DOWN;
      float costhetastar_ClusteredMet_DOWN;
      float Phi1_ClusteredMet_DOWN;
      // unc end
      
      Float_t q_1, q_2;
      Float_t m_sv;
      Float_t pt_sv;
      Float_t eta_sv;
      Float_t phi_sv;
      Int_t   njets;
      Float_t jpt_1, jeta_1, jphi_1;
      Float_t jpt_2, jeta_2, jphi_2;
      Float_t jpt_1_JESUp, jpt_1_JESDown;
      Float_t jpt_2_JESUp, jpt_2_JESDown;
      Float_t pt_2;
      Float_t phi_2;
      Float_t eta_2;
      Float_t m_2;
      Float_t pt_1;
      Float_t phi_1;
      Float_t eta_1;
      Float_t m_1;
      
      Float_t tau1_pt;
      Float_t tau1_eta;
      Float_t tau1_phi;
      Float_t tau1_m;
      
      Float_t tau2_pt;
      Float_t tau2_eta;
      Float_t tau2_phi;
      Float_t tau2_m;

      // Kinematics of objects for MET cluster and uncluster uncertainties
      Float_t tau1_pt_UncMet_UP;
      Float_t tau1_eta_UncMet_UP;
      Float_t tau1_phi_UncMet_UP;
      Float_t tau1_m_UncMet_UP;
      Float_t tau2_pt_UncMet_UP;
      Float_t tau2_eta_UncMet_UP;
      Float_t tau2_phi_UncMet_UP;
      Float_t tau2_m_UncMet_UP;

      Float_t tau1_pt_UncMet_DOWN;
      Float_t tau1_eta_UncMet_DOWN;
      Float_t tau1_phi_UncMet_DOWN;
      Float_t tau1_m_UncMet_DOWN;
      Float_t tau2_pt_UncMet_DOWN;
      Float_t tau2_eta_UncMet_DOWN;
      Float_t tau2_phi_UncMet_DOWN;
      Float_t tau2_m_UncMet_DOWN;

      Float_t tau1_pt_ClusteredMet_UP;
      Float_t tau1_eta_ClusteredMet_UP;
      Float_t tau1_phi_ClusteredMet_UP;
      Float_t tau1_m_ClusteredMet_UP;
      Float_t tau2_pt_ClusteredMet_UP;
      Float_t tau2_eta_ClusteredMet_UP;
      Float_t tau2_phi_ClusteredMet_UP;
      Float_t tau2_m_ClusteredMet_UP;

      Float_t tau1_pt_ClusteredMet_DOWN;
      Float_t tau1_eta_ClusteredMet_DOWN;
      Float_t tau1_phi_ClusteredMet_DOWN;
      Float_t tau1_m_ClusteredMet_DOWN;
      Float_t tau2_pt_ClusteredMet_DOWN;
      Float_t tau2_eta_ClusteredMet_DOWN;
      Float_t tau2_phi_ClusteredMet_DOWN;
      Float_t tau2_m_ClusteredMet_DOWN;

      Int_t gen_match_1;
      Int_t gen_match_2;

      float decayMode=-999.;
      float decayMode2;
      
      TBranch *b_q_1;
      TBranch *b_q_2;
      TBranch *b_njets;   //!
      TBranch *b_jpt_1;   //!
      TBranch *b_jpt_1_JESUp;   //!
      TBranch *b_jpt_1_JESDown;   //!
      TBranch *b_jeta_1;   //!
      TBranch *b_jphi_1;   //!
      TBranch *b_jpt_2;   //!
      TBranch *b_jpt_2_JESUp;   //!
      TBranch *b_jpt_2_JESDown;   //!
      TBranch *b_jeta_2;   //!
      TBranch *b_jphi_2;   //!
      TBranch *b_m_sv;   //!
      TBranch *b_pt_sv;   //!
      TBranch *b_eta_sv;   //!
      TBranch *b_phi_sv;   //!
      TBranch *b_pt_2;   //!
      TBranch *b_phi_2;   //!
      TBranch *b_eta_2;   //!
      TBranch *b_m_2;   //!
      TBranch *b_pt_1;   //!
      TBranch *b_phi_1;   //!
      TBranch *b_eta_1;   //!
      TBranch *b_m_1;   //!
      TBranch *b_tau1_pt;
      TBranch *b_tau1_phi;
      TBranch *b_tau1_eta;
      TBranch *b_tau1_m;
      TBranch *b_tau2_pt;
      TBranch *b_tau2_phi;
      TBranch *b_tau2_eta;
      TBranch *b_tau2_m;
      // Kinematics of objects for MET cluster and uncluster uncertainties
      TBranch *b_tau1_pt_UncMet_UP;
      TBranch *b_tau1_phi_UncMet_UP;
      TBranch *b_tau1_eta_UncMet_UP;
      TBranch *b_tau1_m_UncMet_UP;
      TBranch *b_tau2_pt_UncMet_UP;
      TBranch *b_tau2_phi_UncMet_UP;
      TBranch *b_tau2_eta_UncMet_UP;
      TBranch *b_tau2_m_UncMet_UP;

      TBranch *b_tau1_pt_UncMet_DOWN;
      TBranch *b_tau1_phi_UncMet_DOWN;
      TBranch *b_tau1_eta_UncMet_DOWN;
      TBranch *b_tau1_m_UncMet_DOWN;
      TBranch *b_tau2_pt_UncMet_DOWN;
      TBranch *b_tau2_phi_UncMet_DOWN;
      TBranch *b_tau2_eta_UncMet_DOWN;
      TBranch *b_tau2_m_UncMet_DOWN;

      TBranch *b_tau1_pt_ClusteredMet_UP;
      TBranch *b_tau1_phi_ClusteredMet_UP;
      TBranch *b_tau1_eta_ClusteredMet_UP;
      TBranch *b_tau1_m_ClusteredMet_UP;
      TBranch *b_tau2_pt_ClusteredMet_UP;
      TBranch *b_tau2_phi_ClusteredMet_UP;
      TBranch *b_tau2_eta_ClusteredMet_UP;
      TBranch *b_tau2_m_ClusteredMet_UP;

      TBranch *b_tau1_pt_ClusteredMet_DOWN;
      TBranch *b_tau1_phi_ClusteredMet_DOWN;
      TBranch *b_tau1_eta_ClusteredMet_DOWN;
      TBranch *b_tau1_m_ClusteredMet_DOWN;
      TBranch *b_tau2_pt_ClusteredMet_DOWN;
      TBranch *b_tau2_phi_ClusteredMet_DOWN;
      TBranch *b_tau2_eta_ClusteredMet_DOWN;
      TBranch *b_tau2_m_ClusteredMet_DOWN;

      TBranch *b_gen_match_1;
      TBranch *b_gen_match_2;
      TBranch *b_t1_decayMode;
      TBranch *b_t2_decayMode;
      TBranch *b_l2_decayMode;
      
      tree->SetBranchAddress("q_1", &q_1, &b_q_1);
      tree->SetBranchAddress("q_2", &q_2, &b_q_2);
      
      tree->SetBranchAddress("tau1_pt", &tau1_pt, &b_tau1_pt);
      tree->SetBranchAddress("tau1_eta", &tau1_eta, &b_tau1_eta);
      tree->SetBranchAddress("tau1_phi", &tau1_phi, &b_tau1_phi);
      tree->SetBranchAddress("tau1_m", &tau1_m, &b_tau1_m);
      
      tree->SetBranchAddress("tau2_pt", &tau2_pt, &b_tau2_pt);
      tree->SetBranchAddress("tau2_eta", &tau2_eta, &b_tau2_eta);
      tree->SetBranchAddress("tau2_phi", &tau2_phi, &b_tau2_phi);
      tree->SetBranchAddress("tau2_m", &tau2_m, &b_tau2_m);

      // Kinematics of objects for MET cluster and uncluster uncertainties
      tree->SetBranchAddress("tau1_pt_UncMet_UP", &tau1_pt_UncMet_UP, &b_tau1_pt_UncMet_UP);
      tree->SetBranchAddress("tau1_eta_UncMet_UP", &tau1_eta_UncMet_UP, &b_tau1_eta_UncMet_UP);
      tree->SetBranchAddress("tau1_phi_UncMet_UP", &tau1_phi_UncMet_UP, &b_tau1_phi_UncMet_UP);
      tree->SetBranchAddress("tau1_m_UncMet_UP", &tau1_m_UncMet_UP, &b_tau1_m_UncMet_UP);
      
      tree->SetBranchAddress("tau2_pt_UncMet_UP", &tau2_pt_UncMet_UP, &b_tau2_pt_UncMet_UP);
      tree->SetBranchAddress("tau2_eta_UncMet_UP", &tau2_eta_UncMet_UP, &b_tau2_eta_UncMet_UP);
      tree->SetBranchAddress("tau2_phi_UncMet_UP", &tau2_phi_UncMet_UP, &b_tau2_phi_UncMet_UP);
      tree->SetBranchAddress("tau2_m_UncMet_UP", &tau2_m_UncMet_UP, &b_tau2_m_UncMet_UP);

      tree->SetBranchAddress("tau1_pt_UncMet_DOWN", &tau1_pt_UncMet_DOWN, &b_tau1_pt_UncMet_DOWN);
      tree->SetBranchAddress("tau1_eta_UncMet_DOWN", &tau1_eta_UncMet_DOWN, &b_tau1_eta_UncMet_DOWN);
      tree->SetBranchAddress("tau1_phi_UncMet_DOWN", &tau1_phi_UncMet_DOWN, &b_tau1_phi_UncMet_DOWN);
      tree->SetBranchAddress("tau1_m_UncMet_DOWN", &tau1_m_UncMet_DOWN, &b_tau1_m_UncMet_DOWN);
      
      tree->SetBranchAddress("tau2_pt_UncMet_DOWN", &tau2_pt_UncMet_DOWN, &b_tau2_pt_UncMet_DOWN);
      tree->SetBranchAddress("tau2_eta_UncMet_DOWN", &tau2_eta_UncMet_DOWN, &b_tau2_eta_UncMet_DOWN);
      tree->SetBranchAddress("tau2_phi_UncMet_DOWN", &tau2_phi_UncMet_DOWN, &b_tau2_phi_UncMet_DOWN);
      tree->SetBranchAddress("tau2_m_UncMet_DOWN", &tau2_m_UncMet_DOWN, &b_tau2_m_UncMet_DOWN);

      tree->SetBranchAddress("tau1_pt_ClusteredMet_UP", &tau1_pt_ClusteredMet_UP, &b_tau1_pt_ClusteredMet_UP);
      tree->SetBranchAddress("tau1_eta_ClusteredMet_UP", &tau1_eta_ClusteredMet_UP, &b_tau1_eta_ClusteredMet_UP);
      tree->SetBranchAddress("tau1_phi_ClusteredMet_UP", &tau1_phi_ClusteredMet_UP, &b_tau1_phi_ClusteredMet_UP);
      tree->SetBranchAddress("tau1_m_ClusteredMet_UP", &tau1_m_ClusteredMet_UP, &b_tau1_m_ClusteredMet_UP);
      
      tree->SetBranchAddress("tau2_pt_ClusteredMet_UP", &tau2_pt_ClusteredMet_UP, &b_tau2_pt_ClusteredMet_UP);
      tree->SetBranchAddress("tau2_eta_ClusteredMet_UP", &tau2_eta_ClusteredMet_UP, &b_tau2_eta_ClusteredMet_UP);
      tree->SetBranchAddress("tau2_phi_ClusteredMet_UP", &tau2_phi_ClusteredMet_UP, &b_tau2_phi_ClusteredMet_UP);
      tree->SetBranchAddress("tau2_m_ClusteredMet_UP", &tau2_m_ClusteredMet_UP, &b_tau2_m_ClusteredMet_UP);

      tree->SetBranchAddress("tau1_pt_ClusteredMet_DOWN", &tau1_pt_ClusteredMet_DOWN, &b_tau1_pt_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau1_eta_ClusteredMet_DOWN", &tau1_eta_ClusteredMet_DOWN, &b_tau1_eta_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau1_phi_ClusteredMet_DOWN", &tau1_phi_ClusteredMet_DOWN, &b_tau1_phi_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau1_m_ClusteredMet_DOWN", &tau1_m_ClusteredMet_DOWN, &b_tau1_m_ClusteredMet_DOWN);
      
      tree->SetBranchAddress("tau2_pt_ClusteredMet_DOWN", &tau2_pt_ClusteredMet_DOWN, &b_tau2_pt_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau2_eta_ClusteredMet_DOWN", &tau2_eta_ClusteredMet_DOWN, &b_tau2_eta_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau2_phi_ClusteredMet_DOWN", &tau2_phi_ClusteredMet_DOWN, &b_tau2_phi_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau2_m_ClusteredMet_DOWN", &tau2_m_ClusteredMet_DOWN, &b_tau2_m_ClusteredMet_DOWN);
      
      tree->SetBranchAddress("njets", &njets, &b_njets);
      tree->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
      tree->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
      tree->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
      tree->SetBranchAddress("jpt_2", &jpt_2, &b_jpt_2);
      tree->SetBranchAddress("jeta_2", &jeta_2, &b_jeta_2);
      tree->SetBranchAddress("jphi_2", &jphi_2, &b_jphi_2);
      tree->SetBranchAddress("m_sv", &m_sv, &b_m_sv);
      tree->SetBranchAddress("pt_sv", &pt_sv, &b_pt_sv);
      tree->SetBranchAddress("eta_sv", &eta_sv, &b_eta_sv);
      tree->SetBranchAddress("phi_sv", &phi_sv, &b_phi_sv);
      tree->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
      tree->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
      tree->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
      tree->SetBranchAddress("m_1", &m_1, &b_m_1);
      tree->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
      tree->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
      tree->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
      tree->SetBranchAddress("m_2", &m_2, &b_m_2);

      tree->SetBranchAddress("gen_match_1",&gen_match_1,&b_gen_match_1);
      tree->SetBranchAddress("gen_match_2",&gen_match_2,&b_gen_match_2);
      if ( channel == "tt" ) tree->SetBranchAddress("t1_decayMode", &decayMode, &b_t1_decayMode);
      if ( channel == "tt" ) tree->SetBranchAddress("t2_decayMode", &decayMode2, &b_t2_decayMode);
      if ( channel != "tt" ) tree->SetBranchAddress("l2_decayMode",&decayMode2, &b_l2_decayMode);

      // Kinematics of objects for MET cluster and uncluster uncertainties
      tree->SetBranchAddress("tau1_pt_UncMet_UP", &tau1_pt_UncMet_UP, &b_tau1_pt_UncMet_UP);
      tree->SetBranchAddress("tau1_eta_UncMet_UP", &tau1_eta_UncMet_UP, &b_tau1_eta_UncMet_UP);
      tree->SetBranchAddress("tau1_phi_UncMet_UP", &tau1_phi_UncMet_UP, &b_tau1_phi_UncMet_UP);
      tree->SetBranchAddress("tau1_m_UncMet_UP", &tau1_m_UncMet_UP, &b_tau1_m_UncMet_UP);
      tree->SetBranchAddress("tau2_pt_UncMet_UP", &tau2_pt_UncMet_UP, &b_tau2_pt_UncMet_UP);
      tree->SetBranchAddress("tau2_eta_UncMet_UP", &tau2_eta_UncMet_UP, &b_tau2_eta_UncMet_UP);
      tree->SetBranchAddress("tau2_phi_UncMet_UP", &tau2_phi_UncMet_UP, &b_tau2_phi_UncMet_UP);
      tree->SetBranchAddress("tau2_m_UncMet_UP", &tau2_m_UncMet_UP, &b_tau2_m_UncMet_UP);

      tree->SetBranchAddress("tau1_pt_UncMet_DOWN", &tau1_pt_UncMet_DOWN, &b_tau1_pt_UncMet_DOWN);
      tree->SetBranchAddress("tau1_eta_UncMet_DOWN", &tau1_eta_UncMet_DOWN, &b_tau1_eta_UncMet_DOWN);
      tree->SetBranchAddress("tau1_phi_UncMet_DOWN", &tau1_phi_UncMet_DOWN, &b_tau1_phi_UncMet_DOWN);
      tree->SetBranchAddress("tau1_m_UncMet_DOWN", &tau1_m_UncMet_DOWN, &b_tau1_m_UncMet_DOWN);
      tree->SetBranchAddress("tau2_pt_UncMet_DOWN", &tau2_pt_UncMet_DOWN, &b_tau2_pt_UncMet_DOWN);
      tree->SetBranchAddress("tau2_eta_UncMet_DOWN", &tau2_eta_UncMet_DOWN, &b_tau2_eta_UncMet_DOWN);
      tree->SetBranchAddress("tau2_phi_UncMet_DOWN", &tau2_phi_UncMet_DOWN, &b_tau2_phi_UncMet_DOWN);
      tree->SetBranchAddress("tau2_m_UncMet_DOWN", &tau2_m_UncMet_DOWN, &b_tau2_m_UncMet_DOWN);


      tree->SetBranchAddress("tau1_pt_ClusteredMet_UP", &tau1_pt_ClusteredMet_UP, &b_tau1_pt_ClusteredMet_UP);
      tree->SetBranchAddress("tau1_eta_ClusteredMet_UP", &tau1_eta_ClusteredMet_UP, &b_tau1_eta_ClusteredMet_UP);
      tree->SetBranchAddress("tau1_phi_ClusteredMet_UP", &tau1_phi_ClusteredMet_UP, &b_tau1_phi_ClusteredMet_UP);
      tree->SetBranchAddress("tau1_m_ClusteredMet_UP", &tau1_m_ClusteredMet_UP, &b_tau1_m_ClusteredMet_UP);
      tree->SetBranchAddress("tau2_pt_ClusteredMet_UP", &tau2_pt_ClusteredMet_UP, &b_tau2_pt_ClusteredMet_UP);
      tree->SetBranchAddress("tau2_eta_ClusteredMet_UP", &tau2_eta_ClusteredMet_UP, &b_tau2_eta_ClusteredMet_UP);
      tree->SetBranchAddress("tau2_phi_ClusteredMet_UP", &tau2_phi_ClusteredMet_UP, &b_tau2_phi_ClusteredMet_UP);
      tree->SetBranchAddress("tau2_m_ClusteredMet_UP", &tau2_m_ClusteredMet_UP, &b_tau2_m_ClusteredMet_UP);

      tree->SetBranchAddress("tau1_pt_ClusteredMet_DOWN", &tau1_pt_ClusteredMet_DOWN, &b_tau1_pt_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau1_eta_ClusteredMet_DOWN", &tau1_eta_ClusteredMet_DOWN, &b_tau1_eta_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau1_phi_ClusteredMet_DOWN", &tau1_phi_ClusteredMet_DOWN, &b_tau1_phi_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau1_m_ClusteredMet_DOWN", &tau1_m_ClusteredMet_DOWN, &b_tau1_m_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau2_pt_ClusteredMet_DOWN", &tau2_pt_ClusteredMet_DOWN, &b_tau2_pt_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau2_eta_ClusteredMet_DOWN", &tau2_eta_ClusteredMet_DOWN, &b_tau2_eta_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau2_phi_ClusteredMet_DOWN", &tau2_phi_ClusteredMet_DOWN, &b_tau2_phi_ClusteredMet_DOWN);
      tree->SetBranchAddress("tau2_m_ClusteredMet_DOWN", &tau2_m_ClusteredMet_DOWN, &b_tau2_m_ClusteredMet_DOWN);
      
      tree->SetBranchAddress("jpt_1_JESUp", &jpt_1_JESUp, &b_jpt_1_JESUp);
      tree->SetBranchAddress("jpt_1_JESDown", &jpt_1_JESDown, &b_jpt_1_JESDown);
      tree->SetBranchAddress("jpt_2_JESUp", &jpt_2_JESUp, &b_jpt_2_JESUp);
      tree->SetBranchAddress("jpt_2_JESDown", &jpt_2_JESDown, &b_jpt_2_JESDown);


      // new branches that will need to be filled
      vector<TBranch*> newBranches;
      newBranches.push_back(tree->Branch("Dbkg_VBF", &mela_Dbkg_VBF));
      newBranches.push_back(tree->Branch("Dbkg_ggH", &mela_Dbkg_ggH));
      newBranches.push_back(tree->Branch("Dbkg_WH", &mela_Dbkg_WH));
      newBranches.push_back(tree->Branch("Dbkg_ZH", &mela_Dbkg_ZH));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF", &ME_sm_VBF));
      newBranches.push_back(tree->Branch("ME_sm_ggH", &ME_sm_ggH));
      newBranches.push_back(tree->Branch("ME_sm_WH", &ME_sm_WH));
      newBranches.push_back(tree->Branch("ME_sm_ZH", &ME_sm_ZH));
      newBranches.push_back(tree->Branch("ME_bkg", &ME_bkg));
      newBranches.push_back(tree->Branch("ME_bkg1", &ME_bkg1));
      newBranches.push_back(tree->Branch("ME_bkg2", &ME_bkg2));
      
      // angles
      newBranches.push_back(tree->Branch("Q2V1", &Q2V1));
      newBranches.push_back(tree->Branch("Q2V2", &Q2V2));
      newBranches.push_back(tree->Branch("costheta1", &costheta1));
      newBranches.push_back(tree->Branch("costheta2", &costheta2));
      newBranches.push_back(tree->Branch("Phi", &Phi));
      newBranches.push_back(tree->Branch("costhetastar", &costhetastar));
      newBranches.push_back(tree->Branch("Phi1", &Phi1));



      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_UP", &mela_Dbkg_VBF_UP));
      newBranches.push_back(tree->Branch("Dbkg_ggH_UP", &mela_Dbkg_ggH_UP));
      newBranches.push_back(tree->Branch("Dbkg_WH_UP", &mela_Dbkg_WH_UP));
      newBranches.push_back(tree->Branch("Dbkg_ZH_UP", &mela_Dbkg_ZH_UP));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_UP", &ME_sm_VBF_UP));
      newBranches.push_back(tree->Branch("ME_sm_ggH_UP", &ME_sm_ggH_UP));
      newBranches.push_back(tree->Branch("ME_sm_WH_UP", &ME_sm_WH_UP));
      newBranches.push_back(tree->Branch("ME_sm_ZH_UP", &ME_sm_ZH_UP));
      newBranches.push_back(tree->Branch("ME_bkg_UP", &ME_bkg_UP));
      newBranches.push_back(tree->Branch("ME_bkg1_UP", &ME_bkg1_UP));
      newBranches.push_back(tree->Branch("ME_bkg2_UP", &ME_bkg2_UP));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_UP", &Q2V1_UP));
      newBranches.push_back(tree->Branch("Q2V2_UP", &Q2V2_UP));
      newBranches.push_back(tree->Branch("costheta1_UP", &costheta1_UP));
      newBranches.push_back(tree->Branch("costheta2_UP", &costheta2_UP));
      newBranches.push_back(tree->Branch("Phi_UP", &Phi_UP));
      newBranches.push_back(tree->Branch("costhetastar_UP", &costhetastar_UP));
      newBranches.push_back(tree->Branch("Phi1_UP", &Phi1_UP));


      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_DM0_UP", &mela_Dbkg_VBF_DM0_UP));
      newBranches.push_back(tree->Branch("Dbkg_ggH_DM0_UP", &mela_Dbkg_ggH_DM0_UP));
      newBranches.push_back(tree->Branch("Dbkg_WH_DM0_UP", &mela_Dbkg_WH_DM0_UP));
      newBranches.push_back(tree->Branch("Dbkg_ZH_DM0_UP", &mela_Dbkg_ZH_DM0_UP));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_DM0_UP", &ME_sm_VBF_DM0_UP));
      newBranches.push_back(tree->Branch("ME_sm_ggH_DM0_UP", &ME_sm_ggH_DM0_UP));
      newBranches.push_back(tree->Branch("ME_sm_WH_DM0_UP", &ME_sm_WH_DM0_UP));
      newBranches.push_back(tree->Branch("ME_sm_ZH_DM0_UP", &ME_sm_ZH_DM0_UP));
      newBranches.push_back(tree->Branch("ME_bkg_DM0_UP", &ME_bkg_DM0_UP));
      newBranches.push_back(tree->Branch("ME_bkg1_DM0_UP", &ME_bkg1_DM0_UP));
      newBranches.push_back(tree->Branch("ME_bkg2_DM0_UP", &ME_bkg2_DM0_UP));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_DM0_UP", &Q2V1_DM0_UP));
      newBranches.push_back(tree->Branch("Q2V2_DM0_UP", &Q2V2_DM0_UP));
      newBranches.push_back(tree->Branch("costheta1_DM0_UP", &costheta1_DM0_UP));
      newBranches.push_back(tree->Branch("costheta2_DM0_UP", &costheta2_DM0_UP));
      newBranches.push_back(tree->Branch("Phi_DM0_UP", &Phi_DM0_UP));
      newBranches.push_back(tree->Branch("costhetastar_DM0_UP", &costhetastar_DM0_UP));
      newBranches.push_back(tree->Branch("Phi1_DM0_UP", &Phi1_DM0_UP));

      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_DM1_UP", &mela_Dbkg_VBF_DM1_UP));
      newBranches.push_back(tree->Branch("Dbkg_ggH_DM1_UP", &mela_Dbkg_ggH_DM1_UP));
      newBranches.push_back(tree->Branch("Dbkg_WH_DM1_UP", &mela_Dbkg_WH_DM1_UP));
      newBranches.push_back(tree->Branch("Dbkg_ZH_DM1_UP", &mela_Dbkg_ZH_DM1_UP));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_DM1_UP", &ME_sm_VBF_DM1_UP));
      newBranches.push_back(tree->Branch("ME_sm_ggH_DM1_UP", &ME_sm_ggH_DM1_UP));
      newBranches.push_back(tree->Branch("ME_sm_WH_DM1_UP", &ME_sm_WH_DM1_UP));
      newBranches.push_back(tree->Branch("ME_sm_ZH_DM1_UP", &ME_sm_ZH_DM1_UP));
      newBranches.push_back(tree->Branch("ME_bkg_DM1_UP", &ME_bkg_DM1_UP));
      newBranches.push_back(tree->Branch("ME_bkg1_DM1_UP", &ME_bkg1_DM1_UP));
      newBranches.push_back(tree->Branch("ME_bkg2_DM1_UP", &ME_bkg2_DM1_UP));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_DM1_UP", &Q2V1_DM1_UP));
      newBranches.push_back(tree->Branch("Q2V2_DM1_UP", &Q2V2_DM1_UP));
      newBranches.push_back(tree->Branch("costheta1_DM1_UP", &costheta1_DM1_UP));
      newBranches.push_back(tree->Branch("costheta2_DM1_UP", &costheta2_DM1_UP));
      newBranches.push_back(tree->Branch("Phi_DM1_UP", &Phi_DM1_UP));
      newBranches.push_back(tree->Branch("costhetastar_DM1_UP", &costhetastar_DM1_UP));
      newBranches.push_back(tree->Branch("Phi1_DM1_UP", &Phi1_DM1_UP));

      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_DM10_UP", &mela_Dbkg_VBF_DM10_UP));
      newBranches.push_back(tree->Branch("Dbkg_ggH_DM10_UP", &mela_Dbkg_ggH_DM10_UP));
      newBranches.push_back(tree->Branch("Dbkg_WH_DM10_UP", &mela_Dbkg_WH_DM10_UP));
      newBranches.push_back(tree->Branch("Dbkg_ZH_DM10_UP", &mela_Dbkg_ZH_DM10_UP));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_DM10_UP", &ME_sm_VBF_DM10_UP));
      newBranches.push_back(tree->Branch("ME_sm_ggH_DM10_UP", &ME_sm_ggH_DM10_UP));
      newBranches.push_back(tree->Branch("ME_sm_WH_DM10_UP", &ME_sm_WH_DM10_UP));
      newBranches.push_back(tree->Branch("ME_sm_ZH_DM10_UP", &ME_sm_ZH_DM10_UP));
      newBranches.push_back(tree->Branch("ME_bkg_DM10_UP", &ME_bkg_DM10_UP));
      newBranches.push_back(tree->Branch("ME_bkg1_DM10_UP", &ME_bkg1_DM10_UP));
      newBranches.push_back(tree->Branch("ME_bkg2_DM10_UP", &ME_bkg2_DM10_UP));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_DM10_UP", &Q2V1_DM10_UP));
      newBranches.push_back(tree->Branch("Q2V2_DM10_UP", &Q2V2_DM10_UP));
      newBranches.push_back(tree->Branch("costheta1_DM10_UP", &costheta1_DM10_UP));
      newBranches.push_back(tree->Branch("costheta2_DM10_UP", &costheta2_DM10_UP));
      newBranches.push_back(tree->Branch("Phi_DM10_UP", &Phi_DM10_UP));
      newBranches.push_back(tree->Branch("costhetastar_DM10_UP", &costhetastar_DM10_UP));
      newBranches.push_back(tree->Branch("Phi1_DM10_UP", &Phi1_DM10_UP));


      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_DOWN", &mela_Dbkg_VBF_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ggH_DOWN", &mela_Dbkg_ggH_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_WH_DOWN", &mela_Dbkg_WH_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ZH_DOWN", &mela_Dbkg_ZH_DOWN));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_DOWN", &ME_sm_VBF_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ggH_DOWN", &ME_sm_ggH_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_WH_DOWN", &ME_sm_WH_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ZH_DOWN", &ME_sm_ZH_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg_DOWN", &ME_bkg_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg1_DOWN", &ME_bkg1_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg2_DOWN", &ME_bkg2_DOWN));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_DOWN", &Q2V1_DOWN));
      newBranches.push_back(tree->Branch("Q2V2_DOWN", &Q2V2_DOWN));
      newBranches.push_back(tree->Branch("costheta1_DOWN", &costheta1_DOWN));
      newBranches.push_back(tree->Branch("costheta2_DOWN", &costheta2_DOWN));
      newBranches.push_back(tree->Branch("Phi_DOWN", &Phi_DOWN));
      newBranches.push_back(tree->Branch("costhetastar_DOWN", &costhetastar_DOWN));
      newBranches.push_back(tree->Branch("Phi1_DOWN", &Phi1_DOWN));

      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_DM0_DOWN", &mela_Dbkg_VBF_DM0_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ggH_DM0_DOWN", &mela_Dbkg_ggH_DM0_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_WH_DM0_DOWN", &mela_Dbkg_WH_DM0_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ZH_DM0_DOWN", &mela_Dbkg_ZH_DM0_DOWN));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_DM0_DOWN", &ME_sm_VBF_DM0_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ggH_DM0_DOWN", &ME_sm_ggH_DM0_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_WH_DM0_DOWN", &ME_sm_WH_DM0_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ZH_DM0_DOWN", &ME_sm_ZH_DM0_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg_DM0_DOWN", &ME_bkg_DM0_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg1_DM0_DOWN", &ME_bkg1_DM0_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg2_DM0_DOWN", &ME_bkg2_DM0_DOWN));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_DM0_DOWN", &Q2V1_DM0_DOWN));
      newBranches.push_back(tree->Branch("Q2V2_DM0_DOWN", &Q2V2_DM0_DOWN));
      newBranches.push_back(tree->Branch("costheta1_DM0_DOWN", &costheta1_DM0_DOWN));
      newBranches.push_back(tree->Branch("costheta2_DM0_DOWN", &costheta2_DM0_DOWN));
      newBranches.push_back(tree->Branch("Phi_DM0_DOWN", &Phi_DM0_DOWN));
      newBranches.push_back(tree->Branch("costhetastar_DM0_DOWN", &costhetastar_DM0_DOWN));
      newBranches.push_back(tree->Branch("Phi1_DM0_DOWN", &Phi1_DM0_DOWN));


      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_DM1_DOWN", &mela_Dbkg_VBF_DM1_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ggH_DM1_DOWN", &mela_Dbkg_ggH_DM1_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_WH_DM1_DOWN", &mela_Dbkg_WH_DM1_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ZH_DM1_DOWN", &mela_Dbkg_ZH_DM1_DOWN));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_DM1_DOWN", &ME_sm_VBF_DM1_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ggH_DM1_DOWN", &ME_sm_ggH_DM1_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_WH_DM1_DOWN", &ME_sm_WH_DM1_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ZH_DM1_DOWN", &ME_sm_ZH_DM1_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg_DM1_DOWN", &ME_bkg_DM1_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg1_DM1_DOWN", &ME_bkg1_DM1_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg2_DM1_DOWN", &ME_bkg2_DM1_DOWN));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_DM1_DOWN", &Q2V1_DM1_DOWN));
      newBranches.push_back(tree->Branch("Q2V2_DM1_DOWN", &Q2V2_DM1_DOWN));
      newBranches.push_back(tree->Branch("costheta1_DM1_DOWN", &costheta1_DM1_DOWN));
      newBranches.push_back(tree->Branch("costheta2_DM1_DOWN", &costheta2_DM1_DOWN));
      newBranches.push_back(tree->Branch("Phi_DM1_DOWN", &Phi_DM1_DOWN));
      newBranches.push_back(tree->Branch("costhetastar_DM1_DOWN", &costhetastar_DM1_DOWN));
      newBranches.push_back(tree->Branch("Phi1_DM1_DOWN", &Phi1_DM1_DOWN));

      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_DM10_DOWN", &mela_Dbkg_VBF_DM10_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ggH_DM10_DOWN", &mela_Dbkg_ggH_DM10_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_WH_DM10_DOWN", &mela_Dbkg_WH_DM10_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ZH_DM10_DOWN", &mela_Dbkg_ZH_DM10_DOWN));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_DM10_DOWN", &ME_sm_VBF_DM10_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ggH_DM10_DOWN", &ME_sm_ggH_DM10_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_WH_DM10_DOWN", &ME_sm_WH_DM10_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ZH_DM10_DOWN", &ME_sm_ZH_DM10_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg_DM10_DOWN", &ME_bkg_DM10_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg1_DM10_DOWN", &ME_bkg1_DM10_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg2_DM10_DOWN", &ME_bkg2_DM10_DOWN));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_DM10_DOWN", &Q2V1_DM10_DOWN));
      newBranches.push_back(tree->Branch("Q2V2_DM10_DOWN", &Q2V2_DM10_DOWN));
      newBranches.push_back(tree->Branch("costheta1_DM10_DOWN", &costheta1_DM10_DOWN));
      newBranches.push_back(tree->Branch("costheta2_DM10_DOWN", &costheta2_DM10_DOWN));
      newBranches.push_back(tree->Branch("Phi_DM10_DOWN", &Phi_DM10_DOWN));
      newBranches.push_back(tree->Branch("costhetastar_DM10_DOWN", &costhetastar_DM10_DOWN));
      newBranches.push_back(tree->Branch("Phi1_DM10_DOWN", &Phi1_DM10_DOWN));


      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_UncMet_UP", &mela_Dbkg_VBF_UncMet_UP));
      newBranches.push_back(tree->Branch("Dbkg_ggH_UncMet_UP", &mela_Dbkg_ggH_UncMet_UP));
      newBranches.push_back(tree->Branch("Dbkg_WH_UncMet_UP", &mela_Dbkg_WH_UncMet_UP));
      newBranches.push_back(tree->Branch("Dbkg_ZH_UncMet_UP", &mela_Dbkg_ZH_UncMet_UP));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_UncMet_UP", &ME_sm_VBF_UncMet_UP));
      newBranches.push_back(tree->Branch("ME_sm_ggH_UncMet_UP", &ME_sm_ggH_UncMet_UP));
      newBranches.push_back(tree->Branch("ME_sm_WH_UncMet_UP", &ME_sm_WH_UncMet_UP));
      newBranches.push_back(tree->Branch("ME_sm_ZH_UncMet_UP", &ME_sm_ZH_UncMet_UP));
      newBranches.push_back(tree->Branch("ME_bkg_UncMet_UP", &ME_bkg_UncMet_UP));
      newBranches.push_back(tree->Branch("ME_bkg1_UncMet_UP", &ME_bkg1_UncMet_UP));
      newBranches.push_back(tree->Branch("ME_bkg2_UncMet_UP", &ME_bkg2_UncMet_UP));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_UncMet_UP", &Q2V1_UncMet_UP));
      newBranches.push_back(tree->Branch("Q2V2_UncMet_UP", &Q2V2_UncMet_UP));
      newBranches.push_back(tree->Branch("costheta1_UncMet_UP", &costheta1_UncMet_UP));
      newBranches.push_back(tree->Branch("costheta2_UncMet_UP", &costheta2_UncMet_UP));
      newBranches.push_back(tree->Branch("Phi_UncMet_UP", &Phi_UncMet_UP));
      newBranches.push_back(tree->Branch("costhetastar_UncMet_UP", &costhetastar_UncMet_UP));
      newBranches.push_back(tree->Branch("Phi1_UncMet_UP", &Phi1_UncMet_UP));

      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_UncMet_DOWN", &mela_Dbkg_VBF_UncMet_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ggH_UncMet_DOWN", &mela_Dbkg_ggH_UncMet_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_WH_UncMet_DOWN", &mela_Dbkg_WH_UncMet_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ZH_UncMet_DOWN", &mela_Dbkg_ZH_UncMet_DOWN));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_UncMet_DOWN", &ME_sm_VBF_UncMet_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ggH_UncMet_DOWN", &ME_sm_ggH_UncMet_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_WH_UncMet_DOWN", &ME_sm_WH_UncMet_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ZH_UncMet_DOWN", &ME_sm_ZH_UncMet_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg_UncMet_DOWN", &ME_bkg_UncMet_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg1_UncMet_DOWN", &ME_bkg1_UncMet_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg2_UncMet_DOWN", &ME_bkg2_UncMet_DOWN));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_UncMet_DOWN", &Q2V1_UncMet_DOWN));
      newBranches.push_back(tree->Branch("Q2V2_UncMet_DOWN", &Q2V2_UncMet_DOWN));
      newBranches.push_back(tree->Branch("costheta1_UncMet_DOWN", &costheta1_UncMet_DOWN));
      newBranches.push_back(tree->Branch("costheta2_UncMet_DOWN", &costheta2_UncMet_DOWN));
      newBranches.push_back(tree->Branch("Phi_UncMet_DOWN", &Phi_UncMet_DOWN));
      newBranches.push_back(tree->Branch("costhetastar_UncMet_DOWN", &costhetastar_UncMet_DOWN));
      newBranches.push_back(tree->Branch("Phi1_UncMet_DOWN", &Phi1_UncMet_DOWN));

      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_ClusteredMet_UP", &mela_Dbkg_VBF_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("Dbkg_ggH_ClusteredMet_UP", &mela_Dbkg_ggH_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("Dbkg_WH_ClusteredMet_UP", &mela_Dbkg_WH_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("Dbkg_ZH_ClusteredMet_UP", &mela_Dbkg_ZH_ClusteredMet_UP));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_ClusteredMet_UP", &ME_sm_VBF_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("ME_sm_ggH_ClusteredMet_UP", &ME_sm_ggH_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("ME_sm_WH_ClusteredMet_UP", &ME_sm_WH_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("ME_sm_ZH_ClusteredMet_UP", &ME_sm_ZH_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("ME_bkg_ClusteredMet_UP", &ME_bkg_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("ME_bkg1_ClusteredMet_UP", &ME_bkg1_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("ME_bkg2_ClusteredMet_UP", &ME_bkg2_ClusteredMet_UP));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_ClusteredMet_UP", &Q2V1_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("Q2V2_ClusteredMet_UP", &Q2V2_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("costheta1_ClusteredMet_UP", &costheta1_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("costheta2_ClusteredMet_UP", &costheta2_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("Phi_ClusteredMet_UP", &Phi_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("costhetastar_ClusteredMet_UP", &costhetastar_ClusteredMet_UP));
      newBranches.push_back(tree->Branch("Phi1_ClusteredMet_UP", &Phi1_ClusteredMet_UP));

      // Unc
      newBranches.push_back(tree->Branch("Dbkg_VBF_ClusteredMet_DOWN", &mela_Dbkg_VBF_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ggH_ClusteredMet_DOWN", &mela_Dbkg_ggH_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_WH_ClusteredMet_DOWN", &mela_Dbkg_WH_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("Dbkg_ZH_ClusteredMet_DOWN", &mela_Dbkg_ZH_ClusteredMet_DOWN));
      // ME
      newBranches.push_back(tree->Branch("ME_sm_VBF_ClusteredMet_DOWN", &ME_sm_VBF_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ggH_ClusteredMet_DOWN", &ME_sm_ggH_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_WH_ClusteredMet_DOWN", &ME_sm_WH_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("ME_sm_ZH_ClusteredMet_DOWN", &ME_sm_ZH_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg_ClusteredMet_DOWN", &ME_bkg_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg1_ClusteredMet_DOWN", &ME_bkg1_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("ME_bkg2_ClusteredMet_DOWN", &ME_bkg2_ClusteredMet_DOWN));
      // angles
      newBranches.push_back(tree->Branch("Q2V1_ClusteredMet_DOWN", &Q2V1_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("Q2V2_ClusteredMet_DOWN", &Q2V2_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("costheta1_ClusteredMet_DOWN", &costheta1_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("costheta2_ClusteredMet_DOWN", &costheta2_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("Phi_ClusteredMet_DOWN", &Phi_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("costhetastar_ClusteredMet_DOWN", &costhetastar_ClusteredMet_DOWN));
      newBranches.push_back(tree->Branch("Phi1_ClusteredMet_DOWN", &Phi1_ClusteredMet_DOWN));

      float mjj = 0;
      newBranches.push_back(tree->Branch("mjj", &mjj));
      float mjj_ClusteredMet_UP = 0;
      newBranches.push_back(tree->Branch("mjj_ClusteredMet_UP", &mjj_ClusteredMet_UP));
      float mjj_ClusteredMet_DOWN = 0;
      newBranches.push_back(tree->Branch("mjj_ClusteredMet_DOWN", &mjj_ClusteredMet_DOWN));

      Long64_t nentries = tree->GetEntries();

      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
	if ( jentry % 100000 == 0 ) std::cout << jentry << std::endl;
	nb = tree->GetEntry(jentry);   nbytes += nb;
	
	mjj = 0;
	ME_sm_VBF = -100; // ME for SM process VBF H->tt
	ME_sm_ggH = -100; // ME for ggH + 2 jets
	ME_sm_WH  = -100; // ME for WH (W->jj)
	ME_sm_ZH  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1 = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2 = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH = -100; // same for ggH and Ztt
	mela_Dbkg_WH  = -100; // same for WH and Ztt
	mela_Dbkg_ZH  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1 = -100;
	Q2V2 = -100;
	costheta1    = -100;
	costheta2    = -100;
	Phi          = -100;
	costhetastar = -100;
	Phi1         = -100;


	ME_sm_VBF_UP = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_UP = -100; // ME for ggH + 2 jets
	ME_sm_WH_UP  = -100; // ME for WH (W->jj)
	ME_sm_ZH_UP  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_UP = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_UP = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_UP  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_UP = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_UP = -100; // same for ggH and Ztt
	mela_Dbkg_WH_UP  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_UP  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_UP = -100;
	Q2V2_UP = -100;
	costheta1_UP    = -100;
	costheta2_UP    = -100;
	Phi_UP          = -100;
	costhetastar_UP = -100;
	Phi1_UP         = -100;

	ME_sm_VBF_DM0_UP = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_DM0_UP = -100; // ME for ggH + 2 jets
	ME_sm_WH_DM0_UP  = -100; // ME for WH (W->jj)
	ME_sm_ZH_DM0_UP  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_DM0_UP = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_DM0_UP = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_DM0_UP  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_DM0_UP = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_DM0_UP = -100; // same for ggH and Ztt
	mela_Dbkg_WH_DM0_UP  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_DM0_UP  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_DM0_UP = -100;
	Q2V2_DM0_UP = -100;
	costheta1_DM0_UP    = -100;
	costheta2_DM0_UP    = -100;
	Phi_DM0_UP          = -100;
	costhetastar_DM0_UP = -100;
	Phi1_DM0_UP         = -100;

	ME_sm_VBF_DM1_UP = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_DM1_UP = -100; // ME for ggH + 2 jets
	ME_sm_WH_DM1_UP  = -100; // ME for WH (W->jj)
	ME_sm_ZH_DM1_UP  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_DM1_UP = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_DM1_UP = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_DM1_UP  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_DM1_UP = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_DM1_UP = -100; // same for ggH and Ztt
	mela_Dbkg_WH_DM1_UP  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_DM1_UP  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_DM1_UP = -100;
	Q2V2_DM1_UP = -100;
	costheta1_DM1_UP    = -100;
	costheta2_DM1_UP    = -100;
	Phi_DM1_UP          = -100;
	costhetastar_DM1_UP = -100;
	Phi1_DM1_UP         = -100;

	ME_sm_VBF_DM10_UP = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_DM10_UP = -100; // ME for ggH + 2 jets
	ME_sm_WH_DM10_UP  = -100; // ME for WH (W->jj)
	ME_sm_ZH_DM10_UP  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_DM10_UP = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_DM10_UP = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_DM10_UP  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_DM10_UP = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_DM10_UP = -100; // same for ggH and Ztt
	mela_Dbkg_WH_DM10_UP  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_DM10_UP  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_DM10_UP = -100;
	Q2V2_DM10_UP = -100;
	costheta1_DM10_UP    = -100;
	costheta2_DM10_UP    = -100;
	Phi_DM10_UP          = -100;
	costhetastar_DM10_UP = -100;
	Phi1_DM10_UP         = -100;

	ME_sm_VBF_DOWN = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_DOWN = -100; // ME for ggH + 2 jets
	ME_sm_WH_DOWN  = -100; // ME for WH (W->jj)
	ME_sm_ZH_DOWN  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_DOWN = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_DOWN = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_DOWN  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_DOWN = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_DOWN = -100; // same for ggH and Ztt
	mela_Dbkg_WH_DOWN  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_DOWN  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_DOWN = -100;
	Q2V2_DOWN = -100;
	costheta1_DOWN    = -100;
	costheta2_DOWN    = -100;
	Phi_DOWN          = -100;
	costhetastar_DOWN = -100;
	Phi1_DOWN         = -100;

	ME_sm_VBF_DM0_DOWN = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_DM0_DOWN = -100; // ME for ggH + 2 jets
	ME_sm_WH_DM0_DOWN  = -100; // ME for WH (W->jj)
	ME_sm_ZH_DM0_DOWN  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_DM0_DOWN = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_DM0_DOWN = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_DM0_DOWN  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_DM0_DOWN = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_DM0_DOWN = -100; // same for ggH and Ztt
	mela_Dbkg_WH_DM0_DOWN  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_DM0_DOWN  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_DM0_DOWN = -100;
	Q2V2_DM0_DOWN = -100;
	costheta1_DM0_DOWN    = -100;
	costheta2_DM0_DOWN    = -100;
	Phi_DM0_DOWN          = -100;
	costhetastar_DM0_DOWN = -100;
	Phi1_DM0_DOWN         = -100;

	ME_sm_VBF_DM1_DOWN = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_DM1_DOWN = -100; // ME for ggH + 2 jets
	ME_sm_WH_DM1_DOWN  = -100; // ME for WH (W->jj)
	ME_sm_ZH_DM1_DOWN  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_DM1_DOWN = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_DM1_DOWN = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_DM1_DOWN  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_DM1_DOWN = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_DM1_DOWN = -100; // same for ggH and Ztt
	mela_Dbkg_WH_DM1_DOWN  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_DM1_DOWN  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_DM1_DOWN = -100;
	Q2V2_DM1_DOWN = -100;
	costheta1_DM1_DOWN    = -100;
	costheta2_DM1_DOWN    = -100;
	Phi_DM1_DOWN          = -100;
	costhetastar_DM1_DOWN = -100;
	Phi1_DM1_DOWN         = -100;

	ME_sm_VBF_DM10_DOWN = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_DM10_DOWN = -100; // ME for ggH + 2 jets
	ME_sm_WH_DM10_DOWN  = -100; // ME for WH (W->jj)
	ME_sm_ZH_DM10_DOWN  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_DM10_DOWN = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_DM10_DOWN = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_DM10_DOWN  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_DM10_DOWN = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_DM10_DOWN = -100; // same for ggH and Ztt
	mela_Dbkg_WH_DM10_DOWN  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_DM10_DOWN  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_DM10_DOWN = -100;
	Q2V2_DM10_DOWN = -100;
	costheta1_DM10_DOWN    = -100;
	costheta2_DM10_DOWN    = -100;
	Phi_DM10_DOWN          = -100;
	costhetastar_DM10_DOWN = -100;
	Phi1_DM10_DOWN         = -100;

	ME_sm_VBF_UncMet_UP = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_UncMet_UP = -100; // ME for ggH + 2 jets
	ME_sm_WH_UncMet_UP  = -100; // ME for WH (W->jj)
	ME_sm_ZH_UncMet_UP  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_UncMet_UP = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_UncMet_UP = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_UncMet_UP  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_UncMet_UP = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_UncMet_UP = -100; // same for ggH and Ztt
	mela_Dbkg_WH_UncMet_UP  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_UncMet_UP  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_UncMet_UP = -100;
	Q2V2_UncMet_UP = -100;
	costheta1_UncMet_UP    = -100;
	costheta2_UncMet_UP    = -100;
	Phi_UncMet_UP          = -100;
	costhetastar_UncMet_UP = -100;
	Phi1_UncMet_UP         = -100;

	ME_sm_VBF_UncMet_DOWN = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_UncMet_DOWN = -100; // ME for ggH + 2 jets
	ME_sm_WH_UncMet_DOWN  = -100; // ME for WH (W->jj)
	ME_sm_ZH_UncMet_DOWN  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_UncMet_DOWN = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_UncMet_DOWN = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_UncMet_DOWN  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_UncMet_DOWN = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_UncMet_DOWN = -100; // same for ggH and Ztt
	mela_Dbkg_WH_UncMet_DOWN  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_UncMet_DOWN  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_UncMet_DOWN = -100;
	Q2V2_UncMet_DOWN = -100;
	costheta1_UncMet_DOWN    = -100;
	costheta2_UncMet_DOWN    = -100;
	Phi_UncMet_DOWN          = -100;
	costhetastar_UncMet_DOWN = -100;
	Phi1_UncMet_DOWN         = -100;

	mjj_ClusteredMet_UP = 0;
	ME_sm_VBF_ClusteredMet_UP = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_ClusteredMet_UP = -100; // ME for ggH + 2 jets
	ME_sm_WH_ClusteredMet_UP  = -100; // ME for WH (W->jj)
	ME_sm_ZH_ClusteredMet_UP  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_ClusteredMet_UP = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_ClusteredMet_UP = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_ClusteredMet_UP  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_ClusteredMet_UP = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_ClusteredMet_UP = -100; // same for ggH and Ztt
	mela_Dbkg_WH_ClusteredMet_UP  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_ClusteredMet_UP  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_ClusteredMet_UP = -100;
	Q2V2_ClusteredMet_UP = -100;
	costheta1_ClusteredMet_UP    = -100;
	costheta2_ClusteredMet_UP    = -100;
	Phi_ClusteredMet_UP          = -100;
	costhetastar_ClusteredMet_UP = -100;
	Phi1_ClusteredMet_UP         = -100;


	mjj_ClusteredMet_DOWN = 0;
	ME_sm_VBF_ClusteredMet_DOWN = -100; // ME for SM process VBF H->tt
	ME_sm_ggH_ClusteredMet_DOWN = -100; // ME for ggH + 2 jets
	ME_sm_WH_ClusteredMet_DOWN  = -100; // ME for WH (W->jj)
	ME_sm_ZH_ClusteredMet_DOWN  = -100; // ME for ZH (Z->jj)
	
	ME_bkg1_ClusteredMet_DOWN = -100;   // ME for Z+2jets with leading jet being first, trailing second
	ME_bkg2_ClusteredMet_DOWN = -100;   // ME for Z+2jets with trailing jet being first, leading second
	ME_bkg_ClusteredMet_DOWN  = -100;   // Sum of the two above (what we need to use)
	
	mela_Dbkg_VBF_ClusteredMet_DOWN = -100; // ME_sm_VBF / (ME_sm_VBF + ME_bkg) <- normalized probability to separate H->tt and Z->tt
	mela_Dbkg_ggH_ClusteredMet_DOWN = -100; // same for ggH and Ztt
	mela_Dbkg_WH_ClusteredMet_DOWN  = -100; // same for WH and Ztt
	mela_Dbkg_ZH_ClusteredMet_DOWN  = -100; // same for ZH and Ztt

	// angles inputs to MELA
	Q2V1_ClusteredMet_DOWN = -100;
	Q2V2_ClusteredMet_DOWN = -100;
	costheta1_ClusteredMet_DOWN    = -100;
	costheta2_ClusteredMet_DOWN    = -100;
	Phi_ClusteredMet_DOWN          = -100;
	costhetastar_ClusteredMet_DOWN = -100;
	Phi1_ClusteredMet_DOWN         = -100;
	

	if (njets>=2){
	  if (channel=="tt") {	  
	    TLorentzVector tau1, tau2;
	    tau1.SetPtEtaPhiM(tau1_pt, tau1_eta, tau1_phi, tau1_m);
	    tau2.SetPtEtaPhiM(tau2_pt, tau2_eta, tau2_phi, tau2_m);

	    // jet 4-vectors
	    TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0),
	      blank1(0, 0, 0, 0);
	    jet1.SetPtEtaPhiM(jpt_1, jeta_1, jphi_1, 0);
	    jet2.SetPtEtaPhiM(jpt_2, jeta_2, jphi_2, 0);
	    mjj = (jet1 +  jet2).M();
	  
	    // tau lepton 4-vectors
	    TLorentzVector pDaughters1, pDaughters2;
	  
	    TLorentzVector visTau1, visTau2;
	    visTau1.SetPtEtaPhiM(pt_1, eta_1, phi_1, m_1);
	    visTau2.SetPtEtaPhiM(pt_2, eta_2, phi_2, m_2);
	  
	    if ( !trueTau ) {
	      pDaughters1 = visTau1;
	      pDaughters2 = visTau2;
	    } else {
	      pDaughters1 = tau1;
	      pDaughters2 = tau2;
	    }

	    // Determine the signs of the tau leptons      
	    int tauCharge1 = q_1;
	    int tauCharge2 = q_2;
	  
	    if ( tau1.DeltaR(visTau1) > tau1.DeltaR(visTau2) )
	      tauCharge1 = q_2;
	    if ( tau2.DeltaR(visTau2) > tau2.DeltaR(visTau1) ) 
	      tauCharge2 = q_1;
	  
	    calculateME(pDaughters1, pDaughters2, jet1, jet2, tauCharge1, tauCharge2,
			ME_sm_VBF, ME_sm_ggH, ME_sm_WH, ME_sm_ZH, ME_bkg1, ME_bkg2,
			Q2V1, Q2V2, costheta1, costheta2, Phi, costhetastar, Phi1);
	    
	    ME_bkg = ME_bkg1 + ME_bkg2;
	    
	    mela_Dbkg_VBF = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
	    mela_Dbkg_ggH = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
	    mela_Dbkg_WH  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
	    mela_Dbkg_ZH  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	    
	    //*****************************************************
	    // MET and JETs SYSTEMATICS
	    // Taus Pt have been corrected at skimming level (TEC)
	    // TEC propagated to shifted METs
	    //*****************************************************
	    std::cout << "MET Unclustered Energy Up   ---  " << std::endl;
	    TLorentzVector pDaughters1_UncMet_UP, pDaughters2_UncMet_UP;
	    pDaughters1_UncMet_UP.SetPtEtaPhiM(tau1_pt_UncMet_UP, tau1_eta_UncMet_UP, tau1_phi_UncMet_UP, tau1_m_UncMet_UP);
	    pDaughters2_UncMet_UP.SetPtEtaPhiM(tau2_pt_UncMet_UP, tau2_eta_UncMet_UP, tau2_phi_UncMet_UP, tau2_m_UncMet_UP);

	    calculateME(pDaughters1_UncMet_UP, pDaughters2_UncMet_UP, jet1, jet2, tauCharge1, tauCharge2,
			ME_sm_VBF_UncMet_UP, ME_sm_ggH_UncMet_UP, ME_sm_WH_UncMet_UP, ME_sm_ZH_UncMet_UP, ME_bkg1_UncMet_UP, ME_bkg2_UncMet_UP,
			Q2V1_UncMet_UP, Q2V2_UncMet_UP, costheta1_UncMet_UP, costheta2_UncMet_UP, Phi_UncMet_UP, costhetastar_UncMet_UP, Phi1_UncMet_UP);

	    ME_bkg_UncMet_UP = ME_bkg1_UncMet_UP + ME_bkg2_UncMet_UP;
	    
	    mela_Dbkg_VBF_UncMet_UP = ME_sm_VBF_UncMet_UP / ( ME_sm_VBF_UncMet_UP + ME_bkg_UncMet_UP);
	    mela_Dbkg_ggH_UncMet_UP = ME_sm_ggH_UncMet_UP / ( ME_sm_ggH_UncMet_UP + ME_bkg_UncMet_UP);
	    mela_Dbkg_WH_UncMet_UP  = ME_sm_WH_UncMet_UP / ( ME_sm_WH_UncMet_UP + ME_bkg_UncMet_UP);
	    mela_Dbkg_ZH_UncMet_UP  = ME_sm_ZH_UncMet_UP / ( ME_sm_ZH_UncMet_UP + ME_bkg_UncMet_UP);

	    std::cout << "MET Unclustered Energy Down ---  " << std::endl;
	    TLorentzVector pDaughters1_UncMet_DOWN, pDaughters2_UncMet_DOWN;
	    pDaughters1_UncMet_DOWN.SetPtEtaPhiM(tau1_pt_UncMet_DOWN, tau1_eta_UncMet_DOWN, tau1_phi_UncMet_DOWN, tau1_m_UncMet_DOWN);
	    pDaughters2_UncMet_DOWN.SetPtEtaPhiM(tau2_pt_UncMet_DOWN, tau2_eta_UncMet_DOWN, tau2_phi_UncMet_DOWN, tau2_m_UncMet_DOWN);

	    calculateME(pDaughters1_UncMet_DOWN, pDaughters2_UncMet_DOWN, jet1, jet2, tauCharge1, tauCharge2,
			ME_sm_VBF_UncMet_DOWN, ME_sm_ggH_UncMet_DOWN, ME_sm_WH_UncMet_DOWN, ME_sm_ZH_UncMet_DOWN, ME_bkg1_UncMet_DOWN, ME_bkg2_UncMet_DOWN,
			Q2V1_UncMet_DOWN, Q2V2_UncMet_DOWN, costheta1_UncMet_DOWN, costheta2_UncMet_DOWN, Phi_UncMet_DOWN, costhetastar_UncMet_DOWN, Phi1_UncMet_DOWN);

	    ME_bkg_UncMet_DOWN = ME_bkg1_UncMet_DOWN + ME_bkg2_UncMet_DOWN;
	    
	    mela_Dbkg_VBF_UncMet_DOWN = ME_sm_VBF_UncMet_DOWN / ( ME_sm_VBF_UncMet_DOWN + ME_bkg_UncMet_DOWN);
	    mela_Dbkg_ggH_UncMet_DOWN = ME_sm_ggH_UncMet_DOWN / ( ME_sm_ggH_UncMet_DOWN + ME_bkg_UncMet_DOWN);
	    mela_Dbkg_WH_UncMet_DOWN  = ME_sm_WH_UncMet_DOWN / ( ME_sm_WH_UncMet_DOWN + ME_bkg_UncMet_DOWN);
	    mela_Dbkg_ZH_UncMet_DOWN  = ME_sm_ZH_UncMet_DOWN / ( ME_sm_ZH_UncMet_DOWN + ME_bkg_UncMet_DOWN);


	    std::cout << "MET Clustered Energy Up     ---  " << std::endl;
	    TLorentzVector jet1_ClusteredMet_UP(0, 0, 1e-3, 1e-3), jet2_ClusteredMet_UP(0, 0, 1e-3, 1e-3);
	    jet1_ClusteredMet_UP.SetPtEtaPhiM(jpt_1_JESUp, jeta_1, jphi_1, 0);
	    jet2_ClusteredMet_UP.SetPtEtaPhiM(jpt_2_JESUp, jeta_2, jphi_2, 0);
	    mjj_ClusteredMet_UP = (jet1_ClusteredMet_UP +  jet2_ClusteredMet_UP).M();

	    TLorentzVector pDaughters1_ClusteredMet_UP, pDaughters2_ClusteredMet_UP;
	    pDaughters1_ClusteredMet_UP.SetPtEtaPhiM(tau1_pt_ClusteredMet_UP, tau1_eta_ClusteredMet_UP, tau1_phi_ClusteredMet_UP, tau1_m_ClusteredMet_UP);
	    pDaughters2_ClusteredMet_UP.SetPtEtaPhiM(tau2_pt_ClusteredMet_UP, tau2_eta_ClusteredMet_UP, tau2_phi_ClusteredMet_UP, tau2_m_ClusteredMet_UP);

	    calculateME(pDaughters1_ClusteredMet_UP, pDaughters2_ClusteredMet_UP, jet1_ClusteredMet_UP, jet2_ClusteredMet_UP, tauCharge1, tauCharge2,
			ME_sm_VBF_ClusteredMet_UP, ME_sm_ggH_ClusteredMet_UP, ME_sm_WH_ClusteredMet_UP, ME_sm_ZH_ClusteredMet_UP, ME_bkg1_ClusteredMet_UP, ME_bkg2_ClusteredMet_UP,
			Q2V1_ClusteredMet_UP, Q2V2_ClusteredMet_UP, costheta1_ClusteredMet_UP, costheta2_ClusteredMet_UP, Phi_ClusteredMet_UP, costhetastar_ClusteredMet_UP, Phi1_ClusteredMet_UP);

	    ME_bkg_ClusteredMet_UP = ME_bkg1_ClusteredMet_UP + ME_bkg2_ClusteredMet_UP;
	    
	    mela_Dbkg_VBF_ClusteredMet_UP = ME_sm_VBF_ClusteredMet_UP / ( ME_sm_VBF_ClusteredMet_UP + ME_bkg_ClusteredMet_UP);
	    mela_Dbkg_ggH_ClusteredMet_UP = ME_sm_ggH_ClusteredMet_UP / ( ME_sm_ggH_ClusteredMet_UP + ME_bkg_ClusteredMet_UP);
	    mela_Dbkg_WH_ClusteredMet_UP  = ME_sm_WH_ClusteredMet_UP / ( ME_sm_WH_ClusteredMet_UP + ME_bkg_ClusteredMet_UP);
	    mela_Dbkg_ZH_ClusteredMet_UP  = ME_sm_ZH_ClusteredMet_UP / ( ME_sm_ZH_ClusteredMet_UP + ME_bkg_ClusteredMet_UP);

	    std::cout << "MET Clustered Energy Down   ---  " << std::endl;
	    TLorentzVector jet1_ClusteredMet_DOWN(0, 0, 1e-3, 1e-3), jet2_ClusteredMet_DOWN(0, 0, 1e-3, 1e-3);
	    jet1_ClusteredMet_DOWN.SetPtEtaPhiM(jpt_1_JESDown, jeta_1, jphi_1, 0);
	    jet2_ClusteredMet_DOWN.SetPtEtaPhiM(jpt_2_JESDown, jeta_2, jphi_2, 0);
	    mjj_ClusteredMet_DOWN = (jet1_ClusteredMet_DOWN +  jet2_ClusteredMet_DOWN).M();

	    TLorentzVector pDaughters1_ClusteredMet_DOWN, pDaughters2_ClusteredMet_DOWN;
	    pDaughters1_ClusteredMet_DOWN.SetPtEtaPhiM(tau1_pt_ClusteredMet_DOWN, tau1_eta_ClusteredMet_DOWN, tau1_phi_ClusteredMet_DOWN, tau1_m_ClusteredMet_DOWN);
	    pDaughters2_ClusteredMet_DOWN.SetPtEtaPhiM(tau2_pt_ClusteredMet_DOWN, tau2_eta_ClusteredMet_DOWN, tau2_phi_ClusteredMet_DOWN, tau2_m_ClusteredMet_DOWN);

	    calculateME(pDaughters1_ClusteredMet_DOWN, pDaughters2_ClusteredMet_DOWN, jet1_ClusteredMet_DOWN, jet2_ClusteredMet_DOWN, tauCharge1, tauCharge2,
			ME_sm_VBF_ClusteredMet_DOWN, ME_sm_ggH_ClusteredMet_DOWN, ME_sm_WH_ClusteredMet_DOWN, ME_sm_ZH_ClusteredMet_DOWN, ME_bkg1_ClusteredMet_DOWN, ME_bkg2_ClusteredMet_DOWN,
			Q2V1_ClusteredMet_DOWN, Q2V2_ClusteredMet_DOWN, costheta1_ClusteredMet_DOWN, costheta2_ClusteredMet_DOWN, Phi_ClusteredMet_DOWN, costhetastar_ClusteredMet_DOWN, Phi1_ClusteredMet_DOWN);

	    ME_bkg_ClusteredMet_DOWN = ME_bkg1_ClusteredMet_DOWN + ME_bkg2_ClusteredMet_DOWN;
	    
	    mela_Dbkg_VBF_ClusteredMet_DOWN = ME_sm_VBF_ClusteredMet_DOWN / ( ME_sm_VBF_ClusteredMet_DOWN + ME_bkg_ClusteredMet_DOWN);
	    mela_Dbkg_ggH_ClusteredMet_DOWN = ME_sm_ggH_ClusteredMet_DOWN / ( ME_sm_ggH_ClusteredMet_DOWN + ME_bkg_ClusteredMet_DOWN);
	    mela_Dbkg_WH_ClusteredMet_DOWN  = ME_sm_WH_ClusteredMet_DOWN / ( ME_sm_WH_ClusteredMet_DOWN + ME_bkg_ClusteredMet_DOWN);
	    mela_Dbkg_ZH_ClusteredMet_DOWN  = ME_sm_ZH_ClusteredMet_DOWN / ( ME_sm_ZH_ClusteredMet_DOWN + ME_bkg_ClusteredMet_DOWN);
	    

	    if(doES) {
	      double tesSize = 0.012; // 0.6% uncertainty is considered for each decay mode. AN line1275. Confirm with Abdollah, it should be 1.2%
	      double tesUP = 1.0 + tesSize;
	      double tesDOWN = 1.0 - tesSize;
	      //***************************************************************************
	      //********************* Two taus shifted up *********************************
	      //***************************************************************************
	      // for now, only tt channel
	      if (gen_match_2==5 or gen_match_1==5){
		std::cout << "Two UP    ---  ";
		float ES_UP_scale1 = 1.0;
		float ES_UP_scale2 = 1.0;
		if(gen_match_1==5) ES_UP_scale1 = tesUP;
		if(gen_match_2==5) ES_UP_scale2 = tesUP;
		std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
		std::cout << "   tes1: " << ES_UP_scale1;
		std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
		std::cout << "   tes2: " << ES_UP_scale2 << std::endl;
		
		// AN line 1281 : TES uncertainty is applied by shifting the tau 4-vector up and down 0.6%,
		TLorentzVector pDaughters1_scaled, pDaughters2_scaled; 
		pDaughters1_scaled.SetPtEtaPhiM(pDaughters1.Pt()*ES_UP_scale1,pDaughters1.Eta(),pDaughters1.Phi(),pDaughters1.M());
		pDaughters2_scaled.SetPtEtaPhiM(pDaughters2.Pt()*ES_UP_scale2,pDaughters2.Eta(),pDaughters2.Phi(),pDaughters2.M());
		// recomputin composite variables in the analysis.
		calculateME(pDaughters1_scaled, pDaughters2_scaled, jet1, jet2, tauCharge1, tauCharge2,
			    ME_sm_VBF_UP, ME_sm_ggH_UP, ME_sm_WH_UP, ME_sm_ZH_UP, ME_bkg1_UP, ME_bkg2_UP,
			    Q2V1_UP, Q2V2_UP, costheta1_UP, costheta2_UP, Phi_UP, costhetastar_UP, Phi1_UP);
		
		ME_bkg_UP = ME_bkg1_UP + ME_bkg2_UP;
		
		mela_Dbkg_VBF_UP = ME_sm_VBF_UP / ( ME_sm_VBF_UP + ME_bkg_UP);
		mela_Dbkg_ggH_UP = ME_sm_ggH_UP / ( ME_sm_ggH_UP + ME_bkg_UP);
		mela_Dbkg_WH_UP  = ME_sm_WH_UP / ( ME_sm_WH_UP + ME_bkg_UP);
		mela_Dbkg_ZH_UP  = ME_sm_ZH_UP / ( ME_sm_ZH_UP + ME_bkg_UP);
	      } else {
		ME_bkg_UP = ME_bkg1 + ME_bkg2;
		
		mela_Dbkg_VBF_UP = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
		mela_Dbkg_ggH_UP = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
		mela_Dbkg_WH_UP  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
		mela_Dbkg_ZH_UP  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	      }
	      //***************************************************************************
	      //********************** Tau DM0 shifted up *********************************
	      //***************************************************************************
	      if ((gen_match_2==5 && decayMode2==0) or (gen_match_1==5 && decayMode==0)){
		std::cout << "DM0 UP    ---  ";
		float ES_UP_scale1 = 1.0;
		float ES_UP_scale2 = 1.0;
		if(gen_match_1==5 && decayMode==0) ES_UP_scale1 = tesUP;
		if(gen_match_2==5 && decayMode2==0) ES_UP_scale2 = tesUP;
		std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
		std::cout << "   tes1: " << ES_UP_scale1;
		std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
		std::cout << "   tes2: " << ES_UP_scale2 << std::endl;
		
		// AN line 1281 : TES uncertainty is applied by shifting the tau 4-vector up and down 0.6%,
		TLorentzVector pDaughters1_scaled, pDaughters2_scaled; 
		pDaughters1_scaled.SetPtEtaPhiM(pDaughters1.Pt()*ES_UP_scale1,pDaughters1.Eta(),pDaughters1.Phi(),pDaughters1.M());
		pDaughters2_scaled.SetPtEtaPhiM(pDaughters2.Pt()*ES_UP_scale2,pDaughters2.Eta(),pDaughters2.Phi(),pDaughters2.M());
		// recomputin composite variables in the analysis.
		calculateME(pDaughters1_scaled, pDaughters2_scaled, jet1, jet2, tauCharge1, tauCharge2,
			    ME_sm_VBF_DM0_UP, ME_sm_ggH_DM0_UP, ME_sm_WH_DM0_UP, ME_sm_ZH_DM0_UP, ME_bkg1_DM0_UP, ME_bkg2_DM0_UP,
			    Q2V1_DM0_UP, Q2V2_DM0_UP, costheta1_DM0_UP, costheta2_DM0_UP, Phi_DM0_UP, costhetastar_DM0_UP, Phi1_DM0_UP);
		
		ME_bkg_DM0_UP = ME_bkg1_DM0_UP + ME_bkg2_DM0_UP;
		
		mela_Dbkg_VBF_DM0_UP = ME_sm_VBF_DM0_UP / ( ME_sm_VBF_DM0_UP + ME_bkg_DM0_UP);
		mela_Dbkg_ggH_DM0_UP = ME_sm_ggH_DM0_UP / ( ME_sm_ggH_DM0_UP + ME_bkg_DM0_UP);
		mela_Dbkg_WH_DM0_UP  = ME_sm_WH_DM0_UP / ( ME_sm_WH_DM0_UP + ME_bkg_DM0_UP);
		mela_Dbkg_ZH_DM0_UP  = ME_sm_ZH_DM0_UP / ( ME_sm_ZH_DM0_UP + ME_bkg_DM0_UP);
	      } else {
		ME_bkg_DM0_UP = ME_bkg1 + ME_bkg2;
		
		mela_Dbkg_VBF_DM0_UP = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
		mela_Dbkg_ggH_DM0_UP = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
		mela_Dbkg_WH_DM0_UP  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
		mela_Dbkg_ZH_DM0_UP  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	      }
	      //***************************************************************************
	      //********************** Tau DM1 shifted up *********************************
	      //***************************************************************************
	      if ((decayMode==1 && gen_match_1==5) or (decayMode2==1 && gen_match_2==5)){
		std::cout << "DM1 UP    ---  ";
		float ES_UP_scale1 = 1.0;
		float ES_UP_scale2 = 1.0;
		if(gen_match_1==5 && decayMode==0) ES_UP_scale1 = tesUP;
		if(gen_match_2==5 && decayMode2==0) ES_UP_scale2 = tesUP;
		std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
		std::cout << "   tes1: " << ES_UP_scale1;
		std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
		std::cout << "   tes2: " << ES_UP_scale2 << std::endl;
		
		// AN line 1281 : TES uncertainty is applied by shifting the tau 4-vector up and down 0.6%,
		TLorentzVector pDaughters1_scaled, pDaughters2_scaled; 
		pDaughters1_scaled.SetPtEtaPhiM(pDaughters1.Pt()*ES_UP_scale1,pDaughters1.Eta(),pDaughters1.Phi(),pDaughters1.M());
		pDaughters2_scaled.SetPtEtaPhiM(pDaughters2.Pt()*ES_UP_scale2,pDaughters2.Eta(),pDaughters2.Phi(),pDaughters2.M());
		// recomputin composite variables in the analysis.
		calculateME(pDaughters1_scaled, pDaughters2_scaled, jet1, jet2, tauCharge1, tauCharge2,
			    ME_sm_VBF_DM1_UP, ME_sm_ggH_DM1_UP, ME_sm_WH_DM1_UP, ME_sm_ZH_DM1_UP, ME_bkg1_DM1_UP, ME_bkg2_DM1_UP,
			    Q2V1_DM1_UP, Q2V2_DM1_UP, costheta1_DM1_UP, costheta2_DM1_UP, Phi_DM1_UP, costhetastar_DM1_UP, Phi1_DM1_UP);
		
		ME_bkg_DM1_UP = ME_bkg1_DM1_UP + ME_bkg2_DM1_UP;
		
		mela_Dbkg_VBF_DM1_UP = ME_sm_VBF_DM1_UP / ( ME_sm_VBF_DM1_UP + ME_bkg_DM1_UP);
		mela_Dbkg_ggH_DM1_UP = ME_sm_ggH_DM1_UP / ( ME_sm_ggH_DM1_UP + ME_bkg_DM1_UP);
		mela_Dbkg_WH_DM1_UP  = ME_sm_WH_DM1_UP / ( ME_sm_WH_DM1_UP + ME_bkg_DM1_UP);
		mela_Dbkg_ZH_DM1_UP  = ME_sm_ZH_DM1_UP / ( ME_sm_ZH_DM1_UP + ME_bkg_DM1_UP);
	      } else {
		ME_bkg_DM1_UP = ME_bkg1 + ME_bkg2;
		
		mela_Dbkg_VBF_DM1_UP = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
		mela_Dbkg_ggH_DM1_UP = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
		mela_Dbkg_WH_DM1_UP  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
		mela_Dbkg_ZH_DM1_UP  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	      }
	      //***************************************************************************
	      //********************* Tau DM10 shifted up *********************************
	      //***************************************************************************
	      if ((decayMode2==10 && gen_match_2==5) or (decayMode==10 && gen_match_1==5)){
		std::cout << "DM10 UP    ---  ";
		float ES_UP_scale1 = 1.0;
		float ES_UP_scale2 = 1.0;
		if(gen_match_1==5 && decayMode==0) ES_UP_scale1 = tesUP;
		if(gen_match_2==5 && decayMode2==0) ES_UP_scale2 = tesUP;
		std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
		std::cout << "   tes1: " << ES_UP_scale1;
		std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
		std::cout << "   tes2: " << ES_UP_scale2 << std::endl;
		
		// AN line 1281 : TES uncertainty is applied by shifting the tau 4-vector up and down 0.6%,
		TLorentzVector pDaughters1_scaled, pDaughters2_scaled; 
		pDaughters1_scaled.SetPtEtaPhiM(pDaughters1.Pt()*ES_UP_scale1,pDaughters1.Eta(),pDaughters1.Phi(),pDaughters1.M());
		pDaughters2_scaled.SetPtEtaPhiM(pDaughters2.Pt()*ES_UP_scale2,pDaughters2.Eta(),pDaughters2.Phi(),pDaughters2.M());	   
		// recomputin composite variables in the analysis.
		calculateME(pDaughters1_scaled, pDaughters2_scaled, jet1, jet2, tauCharge1, tauCharge2,
			    ME_sm_VBF_DM10_UP, ME_sm_ggH_DM10_UP, ME_sm_WH_DM10_UP, ME_sm_ZH_DM10_UP, ME_bkg1_DM10_UP, ME_bkg2_DM10_UP,
			    Q2V1_DM10_UP, Q2V2_DM10_UP, costheta1_DM10_UP, costheta2_DM10_UP, Phi_DM10_UP, costhetastar_DM10_UP, Phi1_DM10_UP);
		
		ME_bkg_DM10_UP = ME_bkg1_DM10_UP + ME_bkg2_DM10_UP;
		
		mela_Dbkg_VBF_DM10_UP = ME_sm_VBF_DM10_UP / ( ME_sm_VBF_DM10_UP + ME_bkg_DM10_UP);
		mela_Dbkg_ggH_DM10_UP = ME_sm_ggH_DM10_UP / ( ME_sm_ggH_DM10_UP + ME_bkg_DM10_UP);
		mela_Dbkg_WH_DM10_UP  = ME_sm_WH_DM10_UP / ( ME_sm_WH_DM10_UP + ME_bkg_DM10_UP);
		mela_Dbkg_ZH_DM10_UP  = ME_sm_ZH_DM10_UP / ( ME_sm_ZH_DM10_UP + ME_bkg_DM10_UP);
	      } else {
		ME_bkg_DM10_UP = ME_bkg1 + ME_bkg2;
		
		mela_Dbkg_VBF_DM10_UP = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
		mela_Dbkg_ggH_DM10_UP = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
		mela_Dbkg_WH_DM10_UP  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
		mela_Dbkg_ZH_DM10_UP  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	      }
	      //*****************************************************
	      //************* Two taus shifted down *****************
	      //*****************************************************
	      if (gen_match_1==5 or gen_match_2==5){
		std::cout << "Two DOWN  ---  ";
		float ES_DOWN_scale1 = 1.0;
		float ES_DOWN_scale2 = 1.0;
		if (gen_match_1==5) ES_DOWN_scale1 = tesDOWN;
		if (gen_match_2==5) ES_DOWN_scale2 = tesDOWN;
		std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
		std::cout << "   tes1: " << ES_DOWN_scale1;
		std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
		std::cout << "   tes2: " << ES_DOWN_scale2 << std::endl;
		
		// AN line 1281 : TES uncertainty is applied by shifting the tau 4-vector up and down 0.6%,
		TLorentzVector pDaughters1_scaled, pDaughters2_scaled; 
		pDaughters1_scaled.SetPtEtaPhiM(pDaughters1.Pt()*ES_DOWN_scale1,pDaughters1.Eta(),pDaughters1.Phi(),pDaughters1.M());
		pDaughters2_scaled.SetPtEtaPhiM(pDaughters2.Pt()*ES_DOWN_scale2,pDaughters2.Eta(),pDaughters2.Phi(),pDaughters2.M());
		// recomputin composite variables in the analysis.
		calculateME(pDaughters1_scaled, pDaughters2_scaled, jet1, jet2, tauCharge1, tauCharge2,
			    ME_sm_VBF_DOWN, ME_sm_ggH_DOWN, ME_sm_WH_DOWN, ME_sm_ZH_DOWN, ME_bkg1_DOWN, ME_bkg2_DOWN,
			    Q2V1_DOWN, Q2V2_DOWN, costheta1_DOWN, costheta2_DOWN, Phi_DOWN, costhetastar_DOWN, Phi1_DOWN);
		
		ME_bkg_DOWN = ME_bkg1_DOWN + ME_bkg2_DOWN;
		
		mela_Dbkg_VBF_DOWN = ME_sm_VBF_DOWN / ( ME_sm_VBF_DOWN + ME_bkg_DOWN);
		mela_Dbkg_ggH_DOWN = ME_sm_ggH_DOWN / ( ME_sm_ggH_DOWN + ME_bkg_DOWN);
		mela_Dbkg_WH_DOWN  = ME_sm_WH_DOWN / ( ME_sm_WH_DOWN + ME_bkg_DOWN);
		mela_Dbkg_ZH_DOWN  = ME_sm_ZH_DOWN / ( ME_sm_ZH_DOWN + ME_bkg_DOWN);
	      } else {
		ME_bkg_DOWN = ME_bkg1 + ME_bkg2;
		
		mela_Dbkg_VBF_DOWN = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
		mela_Dbkg_ggH_DOWN = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
		mela_Dbkg_WH_DOWN  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
		mela_Dbkg_ZH_DOWN  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	      }
	      //*****************************************************
	      //************* Tau DM0 shifted down  *****************
	      //*****************************************************
	      if ((decayMode==0 && gen_match_1==5) or (decayMode2==0 && gen_match_2==5)){
		float ES_DOWN_scale1 = 1.0;
		float ES_DOWN_scale2 = 1.0;
		if (gen_match_1==5) ES_DOWN_scale1 = tesDOWN;
		if (gen_match_2==5) ES_DOWN_scale2 = tesDOWN;	  
		std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
		std::cout << "   tes1: " << ES_DOWN_scale1;
		std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
		std::cout << "   tes2: " << ES_DOWN_scale2 << std::endl;
		
		// AN line 1281 : TES uncertainty is applied by shifting the tau 4-vector up and down 0.6%,
		TLorentzVector pDaughters1_scaled, pDaughters2_scaled; 
		pDaughters1_scaled.SetPtEtaPhiM(pDaughters1.Pt()*ES_DOWN_scale1,pDaughters1.Eta(),pDaughters1.Phi(),pDaughters1.M());
		pDaughters2_scaled.SetPtEtaPhiM(pDaughters2.Pt()*ES_DOWN_scale2,pDaughters2.Eta(),pDaughters2.Phi(),pDaughters2.M());
		// recomputin composite variables in the analysis.
		calculateME(pDaughters1_scaled, pDaughters2_scaled, jet1, jet2, tauCharge1, tauCharge2,
			    ME_sm_VBF_DM0_DOWN, ME_sm_ggH_DM0_DOWN, ME_sm_WH_DM0_DOWN, ME_sm_ZH_DM0_DOWN, ME_bkg1_DM0_DOWN, ME_bkg2_DM0_DOWN,
			    Q2V1_DM0_DOWN, Q2V2_DM0_DOWN, costheta1_DM0_DOWN, costheta2_DM0_DOWN, Phi_DM0_DOWN, costhetastar_DM0_DOWN, Phi1_DM0_DOWN);
		
		ME_bkg_DM0_DOWN = ME_bkg1_DM0_DOWN + ME_bkg2_DM0_DOWN;
		
		mela_Dbkg_VBF_DM0_DOWN = ME_sm_VBF_DM0_DOWN / ( ME_sm_VBF_DM0_DOWN + ME_bkg_DM0_DOWN);
		mela_Dbkg_ggH_DM0_DOWN = ME_sm_ggH_DM0_DOWN / ( ME_sm_ggH_DM0_DOWN + ME_bkg_DM0_DOWN);
		mela_Dbkg_WH_DM0_DOWN  = ME_sm_WH_DM0_DOWN / ( ME_sm_WH_DM0_DOWN + ME_bkg_DM0_DOWN);
		mela_Dbkg_ZH_DM0_DOWN  = ME_sm_ZH_DM0_DOWN / ( ME_sm_ZH_DM0_DOWN + ME_bkg_DM0_DOWN);
	      } else {
		ME_bkg_DM0_DOWN = ME_bkg1 + ME_bkg2;
		
		mela_Dbkg_VBF_DM0_DOWN = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
		mela_Dbkg_ggH_DM0_DOWN = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
		mela_Dbkg_WH_DM0_DOWN  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
		mela_Dbkg_ZH_DM0_DOWN  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	      }
	      //*****************************************************
	      //************** Tau DM1 shifted down *****************
	      //*****************************************************
	      if ((decayMode==1 && gen_match_1==5) or (decayMode2==1 && gen_match_2==5)){
		std::cout << "DM1 DOWN  ---  ";
		float ES_DOWN_scale1 = 1.0;
		float ES_DOWN_scale2 = 1.0;
		if (gen_match_1==5) ES_DOWN_scale1 = tesDOWN;
		if (gen_match_2==5) ES_DOWN_scale2 = tesDOWN;	  
		std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
		std::cout << "   tes1: " << ES_DOWN_scale1;
		std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
		std::cout << "   tes2: " << ES_DOWN_scale2 << std::endl;
		
		// AN line 1281 : TES uncertainty is applied by shifting the tau 4-vector up and down 0.6%,
		TLorentzVector pDaughters1_scaled, pDaughters2_scaled; 
		pDaughters1_scaled.SetPtEtaPhiM(pDaughters1.Pt()*ES_DOWN_scale1,pDaughters1.Eta(),pDaughters1.Phi(),pDaughters1.M());
		pDaughters2_scaled.SetPtEtaPhiM(pDaughters2.Pt()*ES_DOWN_scale2,pDaughters2.Eta(),pDaughters2.Phi(),pDaughters2.M());
		// recomputin composite variables in the analysis.
		calculateME(pDaughters1_scaled, pDaughters2_scaled, jet1, jet2, tauCharge1, tauCharge2,
			    ME_sm_VBF_DM1_DOWN, ME_sm_ggH_DM1_DOWN, ME_sm_WH_DM1_DOWN, ME_sm_ZH_DM1_DOWN, ME_bkg1_DM1_DOWN, ME_bkg2_DM1_DOWN,
			    Q2V1_DM1_DOWN, Q2V2_DM1_DOWN, costheta1_DM1_DOWN, costheta2_DM1_DOWN, Phi_DM1_DOWN, costhetastar_DM1_DOWN, Phi1_DM1_DOWN);
		
		ME_bkg_DM1_DOWN = ME_bkg1_DM1_DOWN + ME_bkg2_DM1_DOWN;
		
		mela_Dbkg_VBF_DM1_DOWN = ME_sm_VBF_DM1_DOWN / ( ME_sm_VBF_DM1_DOWN + ME_bkg_DM1_DOWN);
		mela_Dbkg_ggH_DM1_DOWN = ME_sm_ggH_DM1_DOWN / ( ME_sm_ggH_DM1_DOWN + ME_bkg_DM1_DOWN);
		mela_Dbkg_WH_DM1_DOWN  = ME_sm_WH_DM1_DOWN / ( ME_sm_WH_DM1_DOWN + ME_bkg_DM1_DOWN);
		mela_Dbkg_ZH_DM1_DOWN  = ME_sm_ZH_DM1_DOWN / ( ME_sm_ZH_DM1_DOWN + ME_bkg_DM1_DOWN);
	      } else {
		ME_bkg_DM1_DOWN = ME_bkg1 + ME_bkg2;
		
		mela_Dbkg_VBF_DM1_DOWN = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
		mela_Dbkg_ggH_DM1_DOWN = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
		mela_Dbkg_WH_DM1_DOWN  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
		mela_Dbkg_ZH_DM1_DOWN  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	      }
	      //*****************************************************
	      //************* Tau DM10 shifted down *****************
	      //*****************************************************
	      if ((decayMode==10 && gen_match_1==5) or (decayMode2==10 && gen_match_2==5)){
		float ES_DOWN_scale1 = 1.0;
		float ES_DOWN_scale2 = 1.0;
		if (gen_match_1==5) ES_DOWN_scale1 = tesDOWN;
		if (gen_match_2==5) ES_DOWN_scale2 = tesDOWN;	  
		std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
		std::cout << "   tes1: " << ES_DOWN_scale1;
		std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
		std::cout << "   tes2: " << ES_DOWN_scale2 << std::endl;
		
		// AN line 1281 : TES uncertainty is applied by shifting the tau 4-vector up and down 0.6%,
		TLorentzVector pDaughters1_scaled, pDaughters2_scaled; 
		pDaughters1_scaled.SetPtEtaPhiM(pDaughters1.Pt()*ES_DOWN_scale1,pDaughters1.Eta(),pDaughters1.Phi(),pDaughters1.M());
		pDaughters2_scaled.SetPtEtaPhiM(pDaughters2.Pt()*ES_DOWN_scale2,pDaughters2.Eta(),pDaughters2.Phi(),pDaughters2.M());
		// recomputin composite variables in the analysis.
		calculateME(pDaughters1_scaled, pDaughters2_scaled, jet1, jet2, tauCharge1, tauCharge2,
			    ME_sm_VBF_DM10_DOWN, ME_sm_ggH_DM10_DOWN, ME_sm_WH_DM10_DOWN, ME_sm_ZH_DM10_DOWN, ME_bkg1_DM10_DOWN, ME_bkg2_DM10_DOWN,
			    Q2V1_DM10_DOWN, Q2V2_DM10_DOWN, costheta1_DM10_DOWN, costheta2_DM10_DOWN, Phi_DM10_DOWN, costhetastar_DM10_DOWN, Phi1_DM10_DOWN);
		
		ME_bkg_DM10_DOWN = ME_bkg1_DM10_DOWN + ME_bkg2_DM10_DOWN;
		
		mela_Dbkg_VBF_DM10_DOWN = ME_sm_VBF_DM10_DOWN / ( ME_sm_VBF_DM10_DOWN + ME_bkg_DM10_DOWN);
		mela_Dbkg_ggH_DM10_DOWN = ME_sm_ggH_DM10_DOWN / ( ME_sm_ggH_DM10_DOWN + ME_bkg_DM10_DOWN);
		mela_Dbkg_WH_DM10_DOWN  = ME_sm_WH_DM10_DOWN / ( ME_sm_WH_DM10_DOWN + ME_bkg_DM10_DOWN);
		mela_Dbkg_ZH_DM10_DOWN  = ME_sm_ZH_DM10_DOWN / ( ME_sm_ZH_DM10_DOWN + ME_bkg_DM10_DOWN);
	      } else {
		ME_bkg_DM10_DOWN = ME_bkg1 + ME_bkg2;
		
		mela_Dbkg_VBF_DM10_DOWN = ME_sm_VBF / ( ME_sm_VBF + ME_bkg);
		mela_Dbkg_ggH_DM10_DOWN = ME_sm_ggH / ( ME_sm_ggH + ME_bkg);
		mela_Dbkg_WH_DM10_DOWN  = ME_sm_WH / ( ME_sm_WH + ME_bkg);
		mela_Dbkg_ZH_DM10_DOWN  = ME_sm_ZH / ( ME_sm_ZH + ME_bkg);
	      }
	    } // end of doES
	  } // end of tt channel
	}
	
	// Fill new branches
	for(auto branchToFill : newBranches) branchToFill->Fill();

      } // loop over events

      dir->cd();
      tree->Write("", TObject::kOverwrite);
      strcpy(treeToUse, stringA);

    }
  }
}

void calculateME(TLorentzVector tau1, TLorentzVector tau2, TLorentzVector jet1, TLorentzVector jet2,
                 int tauCharge1, int tauCharge2,
                 float& ME_sm_VBF, float& ME_sm_ggH, float& ME_sm_WH, float& ME_sm_ZH,
                 float& ME_bkg1, float& ME_bkg2,
                 float& Q2V1, float& Q2V2, float& costheta1, float& costheta2, float& Phi, float& costhetastar, float& Phi1) {

  /*
  std::cout << "Input " << std::endl;
  std::cout << "\t tau1: " << tau1.Pt() << "\t" << tau1.Eta() << "\t" << tau1.Phi() << "\t" << tau1.M() << std::endl;
  std::cout << "\t tau2: " << tau2.Pt() << "\t" << tau2.Eta() << "\t" << tau2.Phi() << "\t" << tau2.M() << std::endl;
  std::cout << "\t jet1: " << jet1.Pt() << "\t" << jet1.Eta() << "\t" << jet1.Phi() << "\t" << jet1.M() << std::endl;
  std::cout << "\t jet2: " << jet2.Pt() << "\t" << jet2.Eta() << "\t" << jet2.Phi() << "\t" << jet2.M() << std::endl;
  std::cout << "And charges: " << tauCharge1 << "\t" << tauCharge2 << std::endl;
  std::cout << std::endl;
  */

  // tau lepton 4-vectors
  SimpleParticleCollection_t daughters;
  daughters.push_back(SimpleParticle_t(15*tauCharge1, tau1)); 
  daughters.push_back(SimpleParticle_t(15*tauCharge2, tau2));

  SimpleParticleCollection_t associated;
  associated.push_back(SimpleParticle_t(0, jet1));
  associated.push_back(SimpleParticle_t(0, jet2));
  
  // SM Higgs hypothesis ===========================================================
  mela.setCandidateDecayMode(TVar::CandidateDecay_ff);
  mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);
  mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
  mela.computeProdP(ME_sm_VBF, false);

  mela.computeVBFAngles(Q2V1, Q2V2, costheta1, costheta2, Phi, costhetastar, Phi1);
  
  //get the ggH + 2jets hypothesis
  mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJQCD);
  mela.selfDHggcoupl[0][gHIGGS_GG_2][0] = 1;
  mela.computeProdP(ME_sm_ggH, false);
  
  //get the WH + 2 jets hypothesis
  mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_WH);
  mela.selfDHzzcoupl[0][gHIGGS_VV_1][0] = 1;
  mela.computeProdP(ME_sm_WH, false);
  
  //get the ZH + 2 jets hypothesis
  mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_ZH);
  mela.selfDHzzcoupl[0][gHIGGS_VV_1][0] = 1;
  mela.computeProdP(ME_sm_ZH, false);
  
  // Z+2jets
  mela.setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
  mela.computeProdP(ME_bkg1, false);
  
  associated.clear();
  associated.push_back(SimpleParticle_t(0, jet2));
  associated.push_back(SimpleParticle_t(0, jet1));
  mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);
  mela.setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
  mela.computeProdP(ME_bkg2, false);
  mela.resetInputEvent();
  return;
}

//Thank you Renee Brun :)                                                                                                                                                  
void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
  //copy all objects and subdirs of directory source as a subdir of the current directory   
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir;
  if(source->GetName()!=parser.stringValue("inputFile")){
    adir = savdir->mkdir(source->GetName());
    std::cout<<"Source name is not outputfile name"<<std::endl;
    adir->cd();
  }
  else{
    //adir = savdir->mkdir("input");                                                     
    adir->cd();
  }

  //loop on all entries of this directory                                                                                                                                                               
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir,parser);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void CopyFile(const char *fname, optutl::CommandLineParser parser) {
  //Copy all objects and subdirs of file fname as a subdir of the current directory                                                                                                                     
  TDirectory *target = gDirectory;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot copy file: %s\n",fname);
    target->cd();
    return;
  }
  target->cd();
  CopyDir(f,parser);
  delete f;
  target->cd();
}
void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew)
{
  //prepare files to be copied                                                                                                         
  if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
  }

  fNew->cd();
  CopyFile(parser.stringValue("inputFile").c_str(),parser);
  fNew->ls();
  fNew->Close();

}


