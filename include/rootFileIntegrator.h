// This class opens a root file at default and will append all the analysis plots for the same module to this file.
// Also it has functions to make compare plots and write to that file as well

#ifndef rootFileIntegrator_h
#define rootFileIntegrator_h

#include "makePlots.h"


class rootFileIntegrator{
    
public:
    
    rootFileIntegrator ( char* fname );
    ~rootFileIntegrator();

    
    //public function
    void Integrate( TObject *obj , string name);
    void CopyDir(TDirectory *source);
    void CopyFile(const char *fname);
    void ComparePlots();

    //public parameter
    TFile *outfile;
    string input_fileName;



private:

    
    
};

#endif
	
	
