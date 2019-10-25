// Arthor: Hsin-Yeh Wu
// Email : thankyouyou06@gmail.com

// This class is the main class for analyzing the rpi data 

#ifndef makePlots_h
#define makePlots_h

#include "ntuplizer.h"
#include "TChain.h"
#include "TH2Poly.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "PlotSetting.h"
#include "TLatex.h"
#include <string>
#include <utility> //std::pair
#include <map>     //std::map

using namespace std;

class makePlots{
public:

    makePlots (TChain* chain1, TChain* chain2);
    ~makePlots();

    
    //public function
    void Init( string pedfile, string gainfile, string noisyfile, bool batch_flag );
    void pedestalPlotter();
    void sweepPlotter();
    void const_injPlotter();
    void cosmicAnalyzer();
    void Gain_factor_producer();
    void Pulse_display( int displayChannel = -1 , int acq_type = 0, int lowerR = -1, int upperR = -1 , int startEv = 0 );
  
    //public parameter
    string input_fileName;
    bool subPed_flag;
    bool integrate_flag;
    bool oneChannelInjection_flag;
    bool externalChargeInjection_flag;

    //yaml parameter
    int injCh;
    string acquisitionType;
    string moduleNumber;
    int injChip;
  
    //PlotSetting class 
    PlotSetting P;

  
private:

    // private function
    double  mipConverter( double hg_SubPed, double lg_SubPed, double tot , int channel );
    bool    totFireCheck(int eventID);
    bool    totFireCheck_chip(int chip, double lg, double tot);
    int     ringPositionFinder( int inj_channel, int channel);
    double  CMCalculator( double **sig_SubPed, int *TS );
    double* CMCalculator_v2(double **sig_SubPed, int chip );
    void    Pedestal_CM_Subtractor( int chip );
    bool    mipSigCheck( double *sig, int *TS );
    void    Crosstalk(Int_t ch);
    void    deallocate();

    // src_txtfile reader 
    void    yamlReader();
    void    GainFactorReader( string gainfile );
    void    noisyChannelReader( string noisyFileName );
    void    read_P_and_N(string ped_file);
    void    readmap();
    void    readSensor2Hexaboard();

    // init functions
    void    init_outputFile();
    void    init_analysisParameter();
    void    init_rootBranch();
    void    init_rootDir();
    void    init_histo();
    void    InitTH2Poly(TH2Poly& poly);  //Give frame to TH2Poly

    // plot function
    void    pulsePlotter( double *sig, int *TS, int ev, int ichip, int ich, int lowerR, int upperR );
    void    Pedestal_poly();
    void    injectionPlots();
    void    oneChannelInjection_injectionPlots();
    void    injectionPlots_allCh();
    void    XTalk_poly();
    void    Xtalk_1D();
    void    toa_plot();
    void    const_injPlots();

    // fit function
    void    fit_pedestalHisto();
    void    fit_const_inj_Histo();

    // output function
    void    output_xtalkCoupling();
    void    output_xtalkCoupling_all(bool start , bool end, int channel);
	

    //private parameter
    TApplication   *app;
    TCanvas        *c;
    TTree          *Chain1;
    TTree          *Chain2;
    TFile          *outfile;
    int            cross_ch_FirstRing[NCHIP][6];
    double         **hg_sig;
    double         **lg_sig;
    char           title[200];
    int            TotalEntries;
    int            Nevents;
    char           plot_dir[100];
    TLatex         latex;
    int            goodEventCount;
    int            goodEventCount_chip[NCHIP];

    // TDirectories
    TDirectory *cdinjCh;
    TDirectory *cdallCh;
    TDirectory *cdinj;
    TDirectory *cdPedestal_histo;
    TDirectory *cdtot;
    TDirectory *cdmip;
    TDirectory *cdhglgPedestal;
    TDirectory *cdPedestal_poly;
    TDirectory *cdhgPedestal;
    TDirectory *cdlgPedestal;
    TDirectory *cdhgNoise;
    TDirectory *cdlgNoise;
    TDirectory *cdhgChi;
    TDirectory *cdlgChi;

    /// Histograms
    TH1D *h_hgPedestal[NSCA][NCHANNEL];
    TH1D *h_lgPedestal[NSCA][NCHANNEL];
    TH1D *h_tot[NCHANNEL];
    TH1D *h_mip[NCHANNEL];
    TH1D *h_xtalkCoupling[NCHANNEL];
    TH1D *h_toarise[NCHANNEL];
    TH1D *h_toafall[NCHANNEL];
    TH1D *h_toaf_r[NCHANNEL];
    TH1D *h_xtalkCoupling_Ring_1chip;
    
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
    //double ADC2MIP = 0.0227;
    double HG2DAC[NCHIP][NCH];
    double LGTP_default = 900;

    //noisy parameter
    vector<int> noisyChannel;
    bool noisy_flag[NCHANNEL];

    //analysis parameter
    double XTalkCoupling_Ring_1Chip_average;
    double fit_intersept;
    double fit_slope;
    double ringChannelCount;
    double *XTalkCoupling_Average;    
    double *dac_ctrl;                 
    double *hg_NoisyChannel;          
    double **hg_allCh;
    double **lg_allCh;      
    double **tot_allCh;
    double **toaf_allCh;
    double **toar_allCh;
    double **toaf_r_allCh;
    double **mip_allCh;
    double **mip_allCh_goodEvent;
    double **XTalkCoupling;  
    double **hgFitMean;      
    double **hgFitSigma;
    double **hgFitChisquare;
    double **lgFitMean;
    double **lgFitSigma;
    double **lgFitChisquare;
    double *mipFitMean;    
    double *mipFitSigma;
    double *mipFitChisquare;
    double *xtalkCouplingFitMean;    
    double *xtalkCouplingFitSigma;
    double *xtalkCouplingFitChisquare;     
    float **hgMean;
    float **lgMean;
    float **hgSigma;
    float **lgSigma;
    double **hg_SubPed;
    double **lg_SubPed;
    double **hg_SubPedCM;
    double **lg_SubPedCM;
    double ***mip_Ring_4Chip;
    double ***XTalkCoupling_Ring_4Chip;
    double **mip_Ring_1Chip;
    double **XTalkCoupling_Ring_1Chip;

  
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

  
    // map < key = chip*32+ch/2 , pair <x, y> > 
    map<int,pair < double,double > > CHmap;

    // sensor2hexaboard;
    int sensor2hexaboard[NCHANNEL/2];
};

#endif
