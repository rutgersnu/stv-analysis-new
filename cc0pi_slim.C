#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

void cc0pi_slim(TString filename) {

  TChain tc("stv_tree");
  tc.Add(filename);
  assert(tc && tc.IsOpen());

  TObjArray* branches = tc.GetListOfBranches();
  for (int i=0; i<branches->GetEntries(); i++) {
    TString name = branches->At(i)->GetName();
    if (name.Contains("weight_")) {
      std::cout << "DISABLE " << name << std::endl;
      tc.SetBranchStatus(name, 0);
    }
  }

  //Create a new file + a clone of old tree header. Do not copy events
  TString n = "slim.root";
  TFile* newfile = new TFile(n, "recreate");
  TTree* newtree = tc.CloneTree();

  newfile->Write();
  delete newfile;
}

