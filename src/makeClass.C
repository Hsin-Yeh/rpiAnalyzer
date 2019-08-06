#include "TH1.h"
#include "TTree.h"

void makeClass(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  //TApplication *a = new TApplication("a", 0, 0);
  TFile f("module124cosmic/ana_output/Module124cosmic_23-7-2019_22-29_pedestal.root");
  TTree *tt = (TTree *) f.Get("pulseshapeplotter/tree");
  //TTree *tt3 = (TTree *) f.Get("metadatantupler/meta");
    
  tt->MakeClass("tree");    
}
