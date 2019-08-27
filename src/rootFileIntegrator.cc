#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "TFile.h"
#include "TKey.h"

#include "rootFileIntegrator.h"

//Constructor
rootFileIntegrator::rootFileIntegrator( char* fname )
{
    outfile = new TFile(fname,"update");
}

//Destructor
rootFileIntegrator::~rootFileIntegrator()
{
    //delete c; // delete canvas first before deleting app
    //delete app;
}

void rootFileIntegrator::Integrate( TObject* obj, string name )
{
    char title[200];
    //sprintf(title,"root_plot/%s",moduleName);
    sprintf(title,"output.root");
    TFile *outfile = new TFile(title,"UPDATE");
    std::cout << "output file = " <<  title << std::endl << std::endl;

    //obj->cd();
    

    outfile->WriteObject(obj,"test");
    
    outfile->Write();
}

void rootFileIntegrator::CopyDir(TDirectory *source) {
   //copy all objects and subdirs of directory source as a subdir of the current directory   
   source->ls();
   TDirectory *savdir = gDirectory;
   TDirectory *adir = savdir->mkdir(source->GetName());
   adir->cd();
   //loop on all entries of this directory
   TKey *key;
   TIter nextkey(source->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      const char *classname = key->GetClassName();
      TClass *cl = gROOT->GetClass(classname);
      if (!cl) continue;
      if (cl->InheritsFrom("TDirectory")) {
         source->cd(key->GetName());
         TDirectory *subdir = gDirectory;
         adir->cd();
         CopyDir(subdir);
         adir->cd();
      } else if (cl->InheritsFrom("TTree")) {
         TTree *T = (TTree*)source->Get(key->GetName());
         adir->cd();
         TTree *newT = T->CloneTree();
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


void rootFileIntegrator::CopyFile(const char *fname) {
   //Copy all objects and subdirs of file fname as a subdir of the current directory

    char title[200];
    //sprintf(title,"root_plot/%s",moduleName);
    sprintf(title,"output.root");
    TFile *outfile = new TFile(title,"UPDATE");
    std::cout << "output file = " <<  title << std::endl << std::endl;

   TDirectory *target = gDirectory;
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) {
      printf("Cannot copy file: %s\n",fname);
      target->cd();
      return;
   }
   target->cd();
   CopyDir(f);
   delete f;
   target->cd();
}  

