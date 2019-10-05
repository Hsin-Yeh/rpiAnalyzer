#include "makePlots.h"
//#include "rootFileIntegrator.h"
#include "PlotSetting.h"
#include <fstream>
#include <iostream>

using namespace std;

/// default setup
string main_outpath = "./";
string main_datainput = "./data_input.txt";
string pedfile   = "./pedestal";
string gainfile  = "./src_txtfile/TPro_fittingoutput.txt";
string noisyfile = "./src_txtfile/noisyChannels.txt";

/// utility
void main_makePlots();
void main_gainfactor();
void main_help();

/// global parameters 
bool isNumber(string s);
int  displayChannel = -1;
int  anaType = 0;
int  pulseDisplay_type = 0;
int  lowerR = -1, upperR = -1;
int  startEv = 0;
bool batch_flag = true;
bool subPed_flag = true;
bool oneChannelInjection_flag = false;
bool finalGainFactor_flag = false;
bool makePlots_flag = true;
bool help_flag = false;
bool gain_flag = false;

/// -------------------- main -------------------- ///
int main(int argc, char** argv){
  
    string arg_string;
    vector<string> arg_list;
  
    for(int i = 0 ; i < argc ; ++i){
	arg_string = argv[i];
	arg_list.push_back(arg_string);
    }

    int iarg = 1;
    while ( iarg < argc ) {
	string arg = arg_list[iarg];
	
	if ( arg == "-p" ) {
	    anaType = 1;
	    batch_flag = false;
	    if ( isNumber( arg_list[iarg+1] ) ) {
		displayChannel = atoi(arg_list[iarg+1].c_str());
		cout << "display channel = " << displayChannel << endl;
		iarg+=2;
	    }
	    else iarg++;
	}
	else if ( arg == "-i" ) {
	    cout << "display charge injection channel = " << endl;
	    batch_flag = false;
	    anaType = 1;
	    pulseDisplay_type = 1; 
	    iarg++;
	}
	else if ( arg == "-s" ) {
	    anaType = 1;
	    batch_flag = false;
	    pulseDisplay_type = 2;
	    iarg++;
	}
	else if ( arg == "-r" ) {
	    lowerR =  atoi(arg_list[iarg+1].c_str());
	    upperR =  atoi(arg_list[iarg+2].c_str());
	    cout << "Set Y range = " << lowerR << "-" << upperR << endl;
	    iarg+=3;
	}
	else if ( arg == "-n" || arg == "-noSubPed" ) {
	    cout << "disable pedestal subtraction" << endl;
	    subPed_flag = false;
	    iarg++;
	}
	else if ( arg == "-m" || arg == "-mask" ) {
	    cout << "mask channel" << endl;
	    oneChannelInjection_flag = true;
	    iarg++;
	}
	else if ( arg == "-c" || arg == "-cosmic" ) {
	    anaType = 2;
	    iarg++;
	}
	else if ( arg == "-e" || arg == "-event" ) {
	    startEv = atoi(arg_list[iarg+1].c_str());
	    iarg+=2;
	}
	else if ( arg == "-g" || arg == "-gain" ) {
	    gain_flag = true;
	    makePlots_flag = false;
	    if ( arg_list[iarg+1] == "-f" ) {
		finalGainFactor_flag = true;
		iarg++;
	    }
	    iarg++;
	}
	else if ( arg == "-h" || arg == "-help" ) {
	    help_flag = true;
	    makePlots_flag = false;
	    break;
	}
	else {
	    std::cout << "Unknown option... print usage" << std::endl;
	    makePlots_flag = false;
	    help_flag = true;
	    break;
	}
    }
    if ( makePlots_flag ) { main_makePlots(); }
    else if ( gain_flag ) { main_gainfactor(); }
    else if ( help_flag ) { main_help(); }
    return (0);
}


void main_gainfactor() {

    TChain *chain1 = new TChain("treeproducer/sk2cms");
    TChain *chain2 = new TChain("pulseshapeplotter/tree");
    string filename;
    ifstream infile("data_input.txt");
    infile >> filename;
    infile.close();
    if( filename.length() > 2){
	cout << "inputFile = " << filename << endl << endl;
	chain1->Add(filename.c_str());
	chain2->Add(filename.c_str());
    }
    else
	cout << "There is no input root file written in the input.txt!" << endl;

    ntuplizer* n = new ntuplizer(chain1, chain2);
    n->input_fileName = filename;
    n->Init(pedfile, gainfile, noisyfile);
    if ( n->acquisitionType == "sweep" ) {
	if ( !finalGainFactor_flag ) 
	    n->ntupleProducer();
	else
	    n->readNsort_gainFactor();
    }


}


/// -------------------- main_makePlots -------------------- ///
void main_makePlots() {
    TChain *chain1 = new TChain("treeproducer/sk2cms");
    TChain *chain2 = new TChain("pulseshapeplotter/tree");
    string filename;
    ifstream infile("data_input.txt");
    infile >> filename;
    infile.close();
    if( filename.length() > 2){
	cout << "inputFile = " << filename << endl << endl;
	chain1->Add(filename.c_str());
	chain2->Add(filename.c_str());
    }
    else
	cout << "There is no input root file written in the input.txt!" << endl;

    makePlots* M = new makePlots(chain1, chain2);
    M->input_fileName = filename;
    M->oneChannelInjection_flag = oneChannelInjection_flag;
    M->subPed_flag = subPed_flag;
    M->Init( pedfile, gainfile, noisyfile, batch_flag );
    if ( anaType == 0 ) {
	cout << "Processing PlotProducer " << endl << endl;
	if ( M->acquisitionType == "standard" ) { M->pedestalPlotter(); }
	else if ( M->acquisitionType == "sweep" ) { M->sweepPlotter(); }
	else if ( M->acquisitionType == "const_inj" ) { M->const_injPlotter(); }
    }
    else if ( anaType == 1 ) {
	cout << "Processing Pulse Displayer  " << endl << endl;
	M->Pulse_display( displayChannel, pulseDisplay_type, lowerR, upperR, startEv );
    }
    else if ( anaType == 2 ) {
	cout << "Processing Cosmic Analyzer   " << endl << endl;
	M->cosmicAnalyzer();
    }
	
    delete M;
}


/// -------------------- main_help -------------------- ///
void main_help() {
    cout << endl;
    cout << "Arguments for this code" << endl << endl;
    cout << "./makePlots without arguments will generate a root file containing the analysis plots" << endl << endl;
    cout << "-p : HG pulse display. You can add number after it, like ./makePlots -p 82, and it will display the channel 82 ( chip1, channel 18 ). Otherwise it will loop over every channel" << endl << endl;
    cout << "-i : HG pulse display for charge injection channel. It will search the injection channel in the yaml file and display it " << endl << endl;
    cout << "-n : disable pedestal and common mode subtraction" << endl << endl;
    cout << "-r : set the Y range for pulse display" << endl << "You can type ./makePlots -p 82 -n -r 0 500 to get a display for channel 82 without ped&CM subtraction and with Y range 0-500" << endl << endl;
    cout << "-c : run cosmic analyzer" << endl << endl;
    cout << "-m : mask channel" << endl << endl;
}

/// isNumber 
bool isNumber(string s) { 
    for (int i = 0; i < s.length(); i++) 
	if (isdigit(s[i]) == true) 
	    return true; 

    return false; 
} 


