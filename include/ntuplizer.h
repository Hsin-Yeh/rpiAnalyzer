// Arthor: Hsin-Yeh Wu
// Email : thankyouyou06@gmail.com

// This class is the main class for analyzing the rpi data 

#ifndef ntuplizer_h
#define ntuplizer_h

#include "PlotSetting.h"
#include "TChain.h"
#include "TH2Poly.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TLatex.h"
#include <string>
#include <utility> //std::pair
#include <map>     //std::map

using namespace std;

const int NCHIP = 4;
const int NCH = 64;
const int NSCA = 13;
const int NformatCH = 128;
const int NCHANNEL = 256;
const int NRings = 5;
const int MaxTS = 2;

class ntuplizer{
public:

    ntuplizer (TChain* chain1, TChain* chain2);
    ~ntuplizer();

    
    //public function
    void Init( string pedfile, string gainfile, string noisyfile );
    void ntupleProducer();
    void readNsort_gainFactor();

  
    //public parameter
    string input_fileName;
    bool subPed_flag;
    bool integrate_flag;
    bool oneChannelInjection_flag;

    //yaml parameter
    int injCh;
    string acquisitionType;
    string moduleNumber;
    int injChip;
  
private:

    // private function
    double  CMCalculator( double **sig_SubPed, int *TS );
    double* CMCalculator_v2(double **sig_SubPed, int chip );
    void    Pedestal_CM_Subtractor();
    bool    totFireCheck(int MaxTS_sca);
    bool    totFireCheck_chip(int chip, double lg, double tot);
    void    injectionPlots();
    double  findFitEdge(TF1 *Linear_fit, vector<double> x, vector<double> y);
    

    // src_txtfile reader 
    void    read_P_and_N(string ped_file);
    void    yamlReader();

    // init functions
    void    init_outputFile();
    void    init_analysisParameter();
    void    init_rootBranch();
    void    init_outputBranch();

    // output function
    void    output_gainFactor();


    //private parameter
    PlotSetting    *P;
    TCanvas        *c;
    TApplication   *app;
    TTree          *Chain1;
    TTree          *Chain2;
    TFile          *outfile;
    char           title[200];
    int            TotalEntries;
    int            Nevents;
    char           plot_dir[100];

    // TDirectories
    
    /// Histograms
    
    //pedestal parameter
    float          avg_HG[NCHIP][NCH][NSCA];
    float          sigma_HG[NCHIP][NCH][NSCA];
    float          avg_LG[NCHIP][NCH][NSCA];
    float          sigma_LG[NCHIP][NCH][NSCA];


    //gainFactor parameter
    double LG2DAC[NCHIP][NCH];
    double TOT2DAC[NCHIP][NCH];
    double HGTP[NCHIP][NCH];
    double LGTP[NCHIP][NCH];
    double TOTOffSet[NCHIP][NCH];
    double HG2DAC[NCHIP][NCH];
    double LGTP_default = 900;

    //noisy parameter
    vector<int> noisyChannel;

    //analysis parameter
    double      **hg_sig;
    double      **lg_sig;
    double      *tot_sig;
    double  	**hg_allCh;        
    double  	**lg_allCh;        
    double  	**tot_allCh;      
    double  	**toaf_allCh;      
    double  	**toar_allCh;      
    double  	**toaf_r_allCh;
    int         goodEventCount;
    vector<int> goodEvent;
    vector<double> hg_injCh[NCHIP];
    vector<double> lg_injCh[NCHIP];
    vector<double> tot_injCh[NCHIP];
    vector<double> dac_ctrl;

    //fit resulte function
    TF1* Linear_fit_hg;
    TF1* Linear_fit_lg;
    TF1* Linear_fit_tot;

      
    ///////////////////////////////
    // Declaration of leaf types //
    ///////////////////////////////
  
    //one event consists of 4 entries, every entry consists of one chip
    Int_t         event;          // event number, 
    Int_t         chip;           // chip number 
    Int_t         roll;           // rollposition
    Int_t         dacinj;         // injection dac number 
    Int_t         timesamp[13];   // Time sample [sca]
    Int_t         hg[13][64];     // high gain [ sca ] [ CHANNEL ]
    Int_t         lg[13][64];     // low  gain [ sca ] [ CHANNEL ]
    Int_t         tot_fast[64];   // not used 
    Int_t         tot_slow[64];   // usually use this to be the number for tot 
    Int_t         toa_rise[64];   // Timing 
    Int_t         toa_fall[64];   // Timing


    // Declaration of leaf types
    Int_t           eventID;
    vector<int>     *skirocID;
    vector<int>     *boardID;
    vector<int>     *channelID;
    vector<float>   *HighGainTS3_CM0;
    vector<float>   *HighGainTS3_CM3;
    vector<float>   *HighGainADC;
    vector<float>   *HighGainTmax;
    vector<float>   *HighGainChi2;
    vector<float>   *HighGainErrorADC;
    vector<float>   *HighGainErrorTmax;
    vector<int>     *HighGainStatus;
    vector<int>     *HighGainNCalls;
    vector<float>   *LowGainTS3_CM0;
    vector<float>   *LowGainTS3_CM3;
    vector<float>   *LowGainADC;
    vector<float>   *LowGainTmax;
    vector<float>   *LowGainChi2;
    vector<float>   *LowGainErrorADC;
    vector<float>   *LowGainErrorTmax;
    vector<int>     *LowGainStatus;
    vector<int>     *LowGainNCalls;
    vector<int>     *TotSlow;
    vector<int>     *ToaRise;
    vector<int>     *ToaFall;

    // List of branches
    TBranch        *b_eventID;   //!
    TBranch        *b_skirocID;   //!
    TBranch        *b_boardID;   //!
    TBranch        *b_channelID;   //!
    TBranch        *b_HighGainTS3_CM0;   //!
    TBranch        *b_HighGainTS3_CM3;   //!
    TBranch        *b_HighGainADC;   //!
    TBranch        *b_HighGainTmax;   //!
    TBranch        *b_HighGainChi2;   //!
    TBranch        *b_HighGainErrorADC;   //!
    TBranch        *b_HighGainErrorTmax;   //!
    TBranch        *b_HighGainStatus;   //!
    TBranch        *b_HighGainNCalls;   //!
    TBranch        *b_LowGainTS3_CM0;   //!
    TBranch        *b_LowGainTS3_CM3;   //!
    TBranch        *b_LowGainADC;   //!
    TBranch        *b_LowGainTmax;   //!
    TBranch        *b_LowGainChi2;   //!
    TBranch        *b_LowGainErrorADC;   //!
    TBranch        *b_LowGainErrorTmax;   //!
    TBranch        *b_LowGainStatus;   //!
    TBranch        *b_LowGainNCalls;   //!
    TBranch        *b_TotSlow;   //!
    TBranch        *b_ToaRise;   //!
    TBranch        *b_ToaFall;   //!



    // output TTree
    TTree *outT = new TTree("sk2cms","sk2cms");

    // output root Dir
    TDirectory *cdtree;
    
    // List of output variable
    int event_out;
    int chip_out;
    int roll_out;
    int dacinj_out;
    int timesamp_out[13];
    int hg_out[13][64];
    int lg_out[13][64];
    int tot_slow_out[64];
    int tot_fast_out[64];
    int toa_rise_out[64];
    int toa_fall_out[64];

    // map < key = chip*32+ch/2 , pair <x, y> > 
    map<int,pair < double,double > > CHmap;

    // sensor2hexaboard;
    int sensor2hexaboard[NCHANNEL/2];
};

#endif
