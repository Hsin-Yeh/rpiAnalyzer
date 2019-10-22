#include "ntuplizer.h"
#include "PlotSetting.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TGraph.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

//#define DEBUG

ntuplizer::ntuplizer(TChain* chain1, TChain* chain2):Chain1(chain1),Chain2(chain2)
{
    cout << "Constructor of ntuplizer ... \n\n" << endl;   
}

ntuplizer::~ntuplizer()
{
    
}


void ntuplizer::Init(string pedfile, string gainfile, string noisyfile ) {
    
    cout << "----------Init start----------" << endl;
    // read init files and initialize
    yamlReader();
    read_P_and_N( pedfile );

    init_outputFile();
    init_rootBranch();
    init_outputBranch();
    init_analysisParameter(); // always after init_rootBranch();
    
    // init Canvas 
    gROOT->SetBatch("kTRUE");
    app = new TApplication("app",0,0);
    c = new TCanvas();
    cout << "----------Init complete----------" << endl << endl;

    P = new PlotSetting();
    P->root_logon();
}


void ntuplizer::ntupleProducer() {
    

    for(int entry = 0; entry < TotalEntries ; ++entry) {

	if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }
	Chain1 -> GetEntry(entry);
	Chain2 -> GetEntry(entry);

	// Timesample 
	int TS[NSCA];
	int TS0_sca, MaxTS_sca;
	for(int sca = 0 ; sca < NSCA ; sca++) {
	    TS[sca] = timesamp[sca];

	    if (timesamp[sca] == 0) { TS0_sca = sca ; }
	    if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
	}

	/// Passing hg, lg to self define array
	for (int ich = 0; ich < NCH; ich++){
	    for (int sca = 0; sca < NSCA; sca++){
		hg_sig[ chip*64 + ich ][sca] = hg[sca][ich];
		lg_sig[ chip*64 + ich ][sca] = lg[sca][ich];
	    }
	    tot_sig[ chip*64 + ich ] = tot_slow[ich];
	}

	if(subPed_flag) { Pedestal_CM_Subtractor(); } // Ped & CM Subtraction

      
	if ( chip != 3 ) continue;
	
	//if ( totFireCheck(MaxTS_sca) == false ) continue;
	for ( int ichip = 0; ichip < NCHIP; ichip++ ) {
	    int inj_channel = ichip*64 + injCh;
	    hg_injCh[ichip].push_back( hg_sig[inj_channel][MaxTS_sca] );
	    lg_injCh[ichip].push_back( lg_sig[inj_channel][MaxTS_sca] );
	    //cout << event << " " << ichip << " " << hg_sig[inj_channel][MaxTS_sca] << " " << lg_sig[inj_channel][MaxTS_sca] << endl;
	    tot_injCh[ichip].push_back( tot_sig[inj_channel] );
	}
	dac_ctrl.push_back( dacinj );
	goodEventCount++;


#ifdef DEBUG
	cout << "event: " << event << endl;
	for (int ich = 0; ich < NCH; ich++){
	    cout << "chip: " << chip << " ch: " << ich << " == ";
	    for (int sca = 0; sca < NSCA; sca++){
		cout << (int)hg_sig[ich][sca] << " " << (int)lg_sig[ich][sca] << " ";
	    }
	    cout << endl;
	}
	for (int ich = 0; ich < NCH; ich++){
	    cout << "ch: " << ich << " tot_slow: " << tot_slow[ich] << " toa_rise: " << toa_rise[ich] << " toa_fall: " << toa_fall[ich]
		 << " toaf_r: " << toa_fall[ich] - toa_rise[ich] << endl;
	}
#endif
	
    }
    /// End First Loop ///

    injectionPlots();
    /*
    for ( int ichip = 0; ichip < NCHIP; ichip++ ) {
	hg_injCh[ichip].clear();
	lg_injCh[ichip].clear();
	tot_injCh[ichip].clear();
	dac_ctrl.clear();
    }
	
    /// Second Loop --> Exclude the failure injection events ///
    for (int entry = 0; entry < TotalEntries; entry++) {

	if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }
	Chain1 -> GetEntry(entry);
	Chain2 -> GetEntry(entry);

	// Timesample 
	int TS[NSCA];
	int TS0_sca, MaxTS_sca;
	for(int sca = 0 ; sca < NSCA ; sca++) {
	    TS[sca] = timesamp[sca];
	    if (timesamp[sca] == 0) { TS0_sca = sca ; }
	    if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
	}

	/// Passing hg, lg to self define array
	for (int ich = 0; ich < NCH; ich++){
	    for (int sca = 0; sca < NSCA; sca++){
		hg_sig[ chip*64 + ich ][sca] = hg[sca][ich];
		lg_sig[ chip*64 + ich ][sca] = lg[sca][ich];
	    }
	    tot_sig[ chip*64 + ich ] = tot_slow[ich];
	}

	if(subPed_flag) { Pedestal_CM_Subtractor(); } // Ped & CM Subtraction

	if ( chip != 3 ) continue;
	if ( totFireCheck(MaxTS_sca)== false ) continue;
	int hg_residual  = abs ( hg_sig[injCh][MaxTS_sca] - Linear_fit_hg->Eval(dacinj) );
	int lg_residual  = abs ( lg_sig[injCh][MaxTS_sca] - Linear_fit_lg->Eval(dacinj) );
	int tot_residual = abs ( tot_slow[injCh] - Linear_fit_tot->Eval(dacinj) );
	if ( dacinj > 500 ) {
	    if ( tot_residual > 80 ) {
		continue;
	    }
	}
	for ( int ichip = 0; ichip < NCHIP; ichip++ ) {
	    int inj_channel = ichip*64 + injCh;
	    hg_injCh[ichip].push_back( hg_sig[inj_channel][MaxTS_sca] );
	    lg_injCh[ichip].push_back( lg_sig[inj_channel][MaxTS_sca] );
	    tot_injCh[ichip].push_back( tot_sig[inj_channel] );
	}
	dac_ctrl.push_back( dacinj );

	goodEvent.push_back( entry - 3 );
    }

    injectionPlots();

    goodEventCount = 0;
    int ev = 0;
    /// Third Loop --> Fill the ntuples
    for ( int entry = 0; entry < Nevents; entry+=4 ) {
	
	if ( entry != goodEvent.at(goodEventCount) ) continue;

	for ( int ientry = entry; ientry < entry+4; ientry++ ) {
	    Chain1->GetEntry(ientry);
	    event_out = ev;
	    dacinj_out = dacinj;
	    roll_out = roll;
	    chip_out = chip;

	    for ( int sca = 0; sca < NSCA; sca++ ) {
		timesamp_out[sca] = timesamp[sca];
	    }
	
	    for ( int ch = 0; ch < NCH; ch ++ ) {
		for ( int sca = 0; sca < NSCA; sca++ ) {
		    hg_out[sca][ch] = hg[sca][ch];
		    lg_out[sca][ch] = hg[sca][ch];
		}
		tot_slow_out[ch] = tot_slow[ch];
		tot_fast_out[ch] = tot_fast[ch];
		toa_rise_out[ch] = toa_rise[ch];
		toa_fall_out[ch] = toa_fall[ch];
	    }
	    outT->Fill();
	}
	ev++;
	goodEventCount++;
    }
    */
    // Output 
    output_gainFactor();

    outT->Write();
    outfile->Write();
    
}


//void ntuplizer::ntupleProducer_oneChannelInjection() {

    
    
//}

//void ntuplizer::mipProducer() {

//    for ( int ev = 0; ev < 
    
//}


void ntuplizer::init_outputFile() {
    
    /// Set Output Root File
    int start = input_fileName.find_last_of("/");
    int end   = input_fileName.find(".root");
    string outf = input_fileName.substr(start+1,end-start-1);

    sprintf(title,"output_root/%s.root",outf.c_str());
    outfile = new TFile(title,"recreate");
    
    char output_filename[200];
    sprintf(output_filename,"%s",title);
    cout << "output file = " <<  title << endl << endl;

}

void ntuplizer::init_outputBranch() {

    cdtree = outfile->mkdir("treeproducer");
    cdtree->cd();

    // initialize root branch
    outT->Branch("event",&event_out, "event/I");
    outT->Branch("chip",&chip_out, "chip/I");
    outT->Branch("roll",&roll_out, "roll/I");
    outT->Branch("dacinj",&dacinj_out, "dacinj/I");
    outT->Branch("timesamp",timesamp_out, "timesamp[13]/I");
    outT->Branch("hg",hg_out, "hg[13][64]/I");
    outT->Branch("lg",lg_out, "lg[13][64]/I");
    outT->Branch("tot_fast",tot_fast_out, "tot_fast[64]/I");
    outT->Branch("tot_slow",tot_slow_out, "tot_slow[64]/I");
    outT->Branch("toa_rise",toa_rise_out, "toa_rise[64]/I");
    outT->Branch("toa_fall",toa_fall_out, "toa_fall[64]/I");
}

void ntuplizer::init_rootBranch() {
    
    // initialize root branch
    Chain1->SetBranchAddress("event",&event);
    Chain1->SetBranchAddress("chip",&chip);
    Chain1->SetBranchAddress("roll",&roll);
    Chain1->SetBranchAddress("dacinj",&dacinj);
    Chain1->SetBranchAddress("timesamp",&timesamp);
    Chain1->SetBranchAddress("hg",&hg);
    Chain1->SetBranchAddress("lg",&lg);
    Chain1->SetBranchAddress("tot_fast",&tot_fast);
    Chain1->SetBranchAddress("tot_slow",&tot_slow);
    Chain1->SetBranchAddress("toa_rise",&toa_rise);
    Chain1->SetBranchAddress("toa_fall",&toa_fall);


    // Set object pointer
    skirocID = 0;
    boardID = 0;
    channelID = 0;
    HighGainTS3_CM0 = 0;
    HighGainTS3_CM3 = 0;
    HighGainADC = 0;
    HighGainTmax = 0;
    HighGainChi2 = 0;
    HighGainErrorADC = 0;
    HighGainErrorTmax = 0;
    HighGainStatus = 0;
    HighGainNCalls = 0;
    LowGainTS3_CM0 = 0;
    LowGainTS3_CM3 = 0;
    LowGainADC = 0;
    LowGainTmax = 0;
    LowGainChi2 = 0;
    LowGainErrorADC = 0;
    LowGainErrorTmax = 0;
    LowGainStatus = 0;
    LowGainNCalls = 0;
    TotSlow = 0;
    ToaRise = 0;
    ToaFall = 0;

    Chain2->SetBranchAddress("eventID", &eventID, &b_eventID);
    Chain2->SetBranchAddress("skirocID", &skirocID, &b_skirocID);
    Chain2->SetBranchAddress("boardID", &boardID, &b_boardID);
    Chain2->SetBranchAddress("channelID", &channelID, &b_channelID);
    Chain2->SetBranchAddress("HighGainTS3_CM0", &HighGainTS3_CM0, &b_HighGainTS3_CM0);
    Chain2->SetBranchAddress("HighGainTS3_CM3", &HighGainTS3_CM3, &b_HighGainTS3_CM3);
    Chain2->SetBranchAddress("HighGainADC", &HighGainADC, &b_HighGainADC);
    Chain2->SetBranchAddress("HighGainTmax", &HighGainTmax, &b_HighGainTmax);
    Chain2->SetBranchAddress("HighGainChi2", &HighGainChi2, &b_HighGainChi2);
    Chain2->SetBranchAddress("HighGainErrorADC", &HighGainErrorADC, &b_HighGainErrorADC);
    Chain2->SetBranchAddress("HighGainErrorTmax", &HighGainErrorTmax, &b_HighGainErrorTmax);
    Chain2->SetBranchAddress("HighGainStatus", &HighGainStatus, &b_HighGainStatus);
    Chain2->SetBranchAddress("HighGainNCalls", &HighGainNCalls, &b_HighGainNCalls);
    Chain2->SetBranchAddress("LowGainTS3_CM0", &LowGainTS3_CM0, &b_LowGainTS3_CM0);
    Chain2->SetBranchAddress("LowGainTS3_CM3", &LowGainTS3_CM3, &b_LowGainTS3_CM3);
    Chain2->SetBranchAddress("LowGainADC", &LowGainADC, &b_LowGainADC);
    Chain2->SetBranchAddress("LowGainTmax", &LowGainTmax, &b_LowGainTmax);
    Chain2->SetBranchAddress("LowGainChi2", &LowGainChi2, &b_LowGainChi2);
    Chain2->SetBranchAddress("LowGainErrorADC", &LowGainErrorADC, &b_LowGainErrorADC);
    Chain2->SetBranchAddress("LowGainErrorTmax", &LowGainErrorTmax, &b_LowGainErrorTmax);
    Chain2->SetBranchAddress("LowGainStatus", &LowGainStatus, &b_LowGainStatus);
    Chain2->SetBranchAddress("LowGainNCalls", &LowGainNCalls, &b_LowGainNCalls);
    Chain2->SetBranchAddress("TotSlow", &TotSlow, &b_TotSlow);
    Chain2->SetBranchAddress("ToaRise", &ToaRise, &b_ToaRise);
    Chain2->SetBranchAddress("ToaFall", &ToaFall, &b_ToaFall);

    TotalEntries = Chain1->GetEntries();
    Nevents = TotalEntries/NCHIP;
    cout << "Finisth initializing root branches" << endl;
    cout << "Total Events = " << Nevents << endl;
    
}


void ntuplizer::init_analysisParameter(){

    hg_sig = new double*[NCHANNEL];
    lg_sig = new double*[NCHANNEL];
    tot_sig = new double[NCHANNEL];
    hg_allCh        = new double*[NCHANNEL];
    lg_allCh        = new double*[NCHANNEL];
    tot_allCh       = new double*[NCHANNEL];
    toaf_allCh      = new double*[NCHANNEL];
    toar_allCh      = new double*[NCHANNEL];
    toaf_r_allCh    = new double*[NCHANNEL];
    for ( int ich = 0; ich < NCHANNEL; ich++ ) {
	hg_sig[ich] = new double[NSCA];
	lg_sig[ich] = new double[NSCA];
    }
    for(int i = 0; i < NCHANNEL; i++){
	hg_allCh[i]        = new double[Nevents];
	lg_allCh[i]        = new double[Nevents];
	tot_allCh[i]       = new double[Nevents];
	toaf_allCh[i]      = new double[Nevents];
	toar_allCh[i]      = new double[Nevents];
	toaf_r_allCh[i]    = new double[Nevents];
    }
    goodEventCount = 0;
    
}


///
/// ==================== Pedestal_CM_Subtractor ==================== ///
///
void ntuplizer::Pedestal_CM_Subtractor(){
	
    for (int ich = 0; ich < NCH; ich+=2){
	for (int sca = 0; sca < NSCA; sca++){
	    hg_sig[ich][sca] -= avg_HG[chip][ich][sca];  // Pedestal Subtraction
	    lg_sig[ich][sca] -= avg_LG[chip][ich][sca];
	}
    }

    double *hgCM_sca, *lgCM_sca;
    hgCM_sca = CMCalculator_v2( hg_sig, chip ); // Calculate CM for each sca
    lgCM_sca = CMCalculator_v2( lg_sig, chip );
    
    for (int ich = 0; ich < NCH; ich+=2){
	for (int sca = 0; sca < NSCA; sca++){
	    hg_sig[ich][sca] -= hgCM_sca[sca]; // CM subtraction
	    lg_sig[ich][sca] -= lgCM_sca[sca];
	    //cout << hg_sig[ich][sca] << endl;
	}
    }

}


///
/// ==================== CMCalculator_v2 ==================== ///
///
double* ntuplizer::CMCalculator_v2 ( double **sig_subPed, int chip ) {
    // Calculate CM for each TS
    static double meanChipPedestal[NSCA];

    int scaCount[NSCA];
    for (int sca = 0; sca < NSCA; sca++) {
	meanChipPedestal[sca] = 0;
	scaCount[sca] = 0;
    }

    for (int ich = 0; ich < NCHANNEL; ich+=2) {
	int ichannel = ich + chip*NCH;
	bool noisy_flag = false;
		
	// Determine if the channel is a noisy channel
	for(int i = 0; i < noisyChannel.size(); i++){
	    if ( ichannel == noisyChannel.at(i) ) { noisy_flag = true; }
	}
	if (noisy_flag) continue;

	// Calculate mean pedestal for each sca 
	for (int sca = 0; sca < NSCA; sca++) {
	    meanChipPedestal[sca] += sig_subPed[ich][sca];
	    scaCount[sca]++;
	}
    }

    for (int sca = 0; sca < NSCA; sca++) {
	meanChipPedestal[sca] /= scaCount[sca];
    }
    return meanChipPedestal;
}

///
/// ==================== CMCalculator ==================== ///
///
double ntuplizer::CMCalculator( double **sig_subPed, int *TS ) {
    // Common mode calculation, per event, per chip
    // Assume common mode is identical for all time samples
    // Pick 5-9 TS for CM
    int scaCount = 0;
    double meanChipPedestal = 0;
  
    for(int ich = 0; ich < NCH; ich+=2){
	for( int timesample = 0; timesample <= 9; timesample++){
	    int sca = 0;
	    while(true){
		if ( TS[sca] == timesample ) break;
		else sca++;
	    }
	    meanChipPedestal += sig_subPed[sca][ich];
	    scaCount++;
	}
    }
    meanChipPedestal /= scaCount;
    return meanChipPedestal;
}


///
/// ==================== read_P_and_N ==================== ///
///
void ntuplizer::read_P_and_N(string ped_file){

    int end = input_fileName.find("ana_output");
    string pedPath = input_fileName.substr(0,end-1);
    char pedFileName[100];
    sprintf(pedFileName,"%s/pedestal",pedPath.c_str());

    char HG_name[100],LG_name[100];
    sprintf(HG_name,"%sHG.txt",pedFileName);
    sprintf(LG_name,"%sLG.txt",pedFileName);
    ifstream inHG(HG_name);
    ifstream inLG(LG_name);
    if( !inHG.is_open() || !inLG.is_open()){
	cout << "File not found! Neither " << HG_name << " or " << LG_name
	     << " exist!" << endl;
	return;}
    else{
	cout << "pedFile = " << endl;
	cout << "1. " << HG_name << "\n" << "2. "<< LG_name << endl;
	string line;
	int ichip,ich;
	int testeof;

	
	while(true){
	    inHG >> testeof;
	    if( inHG.eof() ) break;
	    else{
		inHG >> ichip; 
		inHG >> ich;
		for(int sca = 0; sca < NSCA ; ++sca){
		    inHG >> avg_HG[ichip][ich][sca];
		    inHG >> sigma_HG[ichip][ich][sca];
		}
	    }
	}

	while(true){
	    inLG >> testeof;
	    if( inLG.eof() ) break;
	    else{
		inLG >> ichip;
		inLG >> ich;
		for(int sca = 0; sca < NSCA ; ++sca){
		    inLG >> avg_LG[ichip][ich][sca];
		    inLG >> sigma_LG[ichip][ich][sca];
		}
	    }
	}

    
	cout << "Reading pedestal file done!" << endl << endl;
	inHG.close();
	inLG.close();
    }  
}

///
/// ==================== yamlReader ==================== ///
///
void ntuplizer::yamlReader(){

    int start = input_fileName.find_last_of("/");
    int end = input_fileName.find("_pedestal.root");
    string f_substr = input_fileName.substr(start+1,end-start-1);
    string rootFileName(f_substr);
    end = input_fileName.find("ana_output");
    f_substr = input_fileName.substr(0,end-1);
    string yamlPath(f_substr);
    char yamlFileName[200];
    sprintf(yamlFileName,"%s/yaml/%s.yaml",yamlPath.c_str(), rootFileName.c_str());
  
    string searchstr;
    string line;
    ifstream yamlFile(yamlFileName);
    if(!yamlFile.is_open()){
	cout << "Did not find injection file " << yamlFileName
	     << ".\n Take this run as pedestal.(Inj_dac = 0)" << endl;
    }
    if(yamlFile.is_open()){
	cout << "yamlFile = " << yamlFileName << endl;
	while( true ) {
	    if ( yamlFile.eof() ) break;
	    getline (yamlFile, line);
	  
	    if ( line.find("channelIds:") != -1 ){
		string tmp;
		start = line.find("[");
		end = line.find("]");
		searchstr = line.substr(start+1,end-start+1);
		injCh = atoi(searchstr.c_str());

		if (start == -1) {
		    getline(yamlFile, line);
		    start = line.find_last_of("-");
		    searchstr = line.erase(0,start+2);
		    injCh = atoi(searchstr.c_str());
		}
				
		cout << "InjCh = " << injCh << endl;
	    }
	    else if ( line.find("acquisitionType") != -1 ){
		start = line.find(":");
		searchstr = line.erase(0, start+2);
		acquisitionType = searchstr;
		cout << "acquisitionType = " << acquisitionType << endl;
	    }
	    else if ( line.find("moduleNumber") != -1 ){
		start = line.find("moduleNumber:");
		end = line.find(",");
		searchstr = line.substr(start+14,end-start-14);
		if ( searchstr.find("'") != -1 ) {
		    start = searchstr.find("'");
		    end = searchstr.find_last_of("'");
		    cout << start << " " << end << endl;
		    searchstr = searchstr.substr(start+1,end-start-1);
		}
		moduleNumber = searchstr;
		cout << "moduleNumber = " << moduleNumber << endl;
	    }
	    else if ( line.find("chipId:") != -1 ) {
		start = line.find(":");
		searchstr = line.substr(start+1,end-start+1);
		injChip = atoi(searchstr.c_str());
		if ( injChip == -1 ) { 
		    cout << "InjChip = " << injChip << endl;
		    oneChannelInjection_flag = false;
		}
		else {
		    injChip = 3 - injChip;
		    cout << "InjChip = " << injChip << endl;
		    oneChannelInjection_flag = true;
		}
	    }
	}
    }

    cout << endl;
}


bool ntuplizer::totFireCheck(int MaxTS_sca){

    bool totFire_flag = true;
    if ( !oneChannelInjection_flag ) {
	for ( int ichip = 0; ichip < NCHIP; ichip++ ){
	    int inj_channel = ichip*NCH + injCh;
	    if ( lg_sig[inj_channel][MaxTS_sca] > 500 && tot_sig[inj_channel] <= 4 ) {
		totFire_flag = false;
	    }
	}
    }
    else {
	int inj_channel = injChip*NCH + injCh;
	if ( lg_sig[inj_channel][MaxTS_sca] > 500 && tot_sig[inj_channel] <= 4 ) {
	    totFire_flag = false;
	}
    }

    return totFire_flag;
}

bool ntuplizer::totFireCheck_chip(int ichip, double lg, double tot){
    bool totFire_flag;
    if ( lg > 500 && tot <= 4 ){ totFire_flag = false; }
    else { totFire_flag = true; }

    return totFire_flag;
}

void ntuplizer::output_gainFactor() {

    ofstream f;
    sprintf(title,"gainfactor_%s.txt",moduleNumber.c_str());
    f.open(title, ios::out | ios::app );

    for ( int ichip = 0; ichip < NCHIP; ichip++ ) {
	if ( oneChannelInjection_flag && ichip!=injChip ) continue;
	f << ichip << " " << injCh << " " << HG2DAC[ichip][injCh] << " " << HGTP[ichip][injCh] << " "  << LG2DAC[ichip][injCh] << " " << LGTP[ichip][injCh] << " " << TOT2DAC[ichip][injCh] << " " << TOTOffSet[ichip][injCh] << endl;
    }
    
}

void ntuplizer::readNsort_gainFactor() {
     
    int ichip, ich;
    string line;
    ifstream ifile;
    sprintf(title,"gainfactor_%s.txt",moduleNumber.c_str());
    ifile.open(title);

    while ( !ifile.eof() ) {
	ifile >> ichip >> ich >> HG2DAC[ichip][ich] >> HGTP[ichip][ich] >> LG2DAC[ichip][ich] >> LGTP[ichip][ich] >> TOT2DAC[ichip][ich] >> TOTOffSet[ichip][ich];
    }
    ifile.close();

    ofstream ofile;
    ofile.open(title);
    for ( int ichip = 0; ichip < NCHIP; ichip++ ) {
	for ( int ich = 0; ich < NCH; ich+=2 ) {
	    ofile << ichip << '\t' << ich << '\t' << fixed << setprecision(6) << HG2DAC[ichip][ich] << '\t' << fixed << setprecision(0) << HGTP[ichip][ich] << '\t' << fixed << setprecision(6) << LG2DAC[ichip][ich] << '\t' << fixed << setprecision(0) << LGTP[ichip][ich] << '\t' << fixed << setprecision(6) << TOT2DAC[ichip][ich] << '\t' << fixed << setprecision(3) << TOTOffSet[ichip][ich] << endl;
	}
    }
}
    




void ntuplizer::injectionPlots(){

    char pltTit[100];
    string Xtit, Ytit, Opt;
    int MkSty,MkClr, LClr, fitmin, fitmax;
    bool Stat, Wait, SavePlot;

    TGraph* gh[NCHIP];
    TGraph* gl[NCHIP];
    TGraph* gtot[NCHIP];
    TGraph* LG2HG[NCHIP];
    TGraph* TOT2LG[NCHIP];

    for(int ichip = 0; ichip < NCHIP; ichip++){

	if ( oneChannelInjection_flag && ichip!= injChip ) continue;
	
	int inj_channel = (ichip*64) + injCh;

	gh[ichip] = new TGraph( goodEventCount, &dac_ctrl[0], &hg_injCh[ichip][0]);
	gl[ichip] = new TGraph( goodEventCount, &dac_ctrl[0], &lg_injCh[ichip][0]);
	gtot[ichip] = new TGraph( goodEventCount, &dac_ctrl[0], &tot_injCh[ichip][0]);
	gh[ichip]->Fit("pol1","Q ROB","",0, 100);
	gl[ichip]->Fit("pol1","Q ROB","",0, 500);
	gtot[ichip]->Fit("pol1","Q ROB","",500, 3000);
	gh[ichip]->GetFunction("pol1")->SetLineColor(3);
	gl[ichip]->GetFunction("pol1")->SetLineColor(3);
	gtot[ichip]->GetFunction("pol1")->SetLineColor(3);
	Linear_fit_hg = gh[ichip]->GetFunction("pol1");
	Linear_fit_lg = gl[ichip]->GetFunction("pol1");
	Linear_fit_tot = gtot[ichip]->GetFunction("pol1");
	HGTP[ichip][injCh] = findFitEdge(Linear_fit_hg, dac_ctrl, hg_injCh[ichip] );
	LGTP[ichip][injCh] = findFitEdge(Linear_fit_lg, dac_ctrl, lg_injCh[ichip] );

	HG2DAC[ichip][injCh] = 1 / Linear_fit_hg->GetParameter(1); 
	LG2DAC[ichip][injCh] = 1 / Linear_fit_lg->GetParameter(1);
	TOT2DAC[ichip][injCh] = 1 / Linear_fit_tot->GetParameter(1);
	TOTOffSet[ichip][injCh] = Linear_fit_tot->GetParameter(0);
	/*
	Wait = 0;
	sprintf(pltTit,"HG_Chip%d",ichip);
	P->GStd(*gh[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait, SavePlot = 0);

	sprintf(pltTit,"LG_Chip%d",ichip);
	P->GStd(*gl[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait, SavePlot = 0);

	sprintf(pltTit,"TOT_Chip%d",ichip);
	P->GStd(*gtot[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait, SavePlot = 0);
	*/
    }

}

double ntuplizer::findFitEdge(TF1 *Linear_fit, vector<double> x, vector<double> y){
    int i = 5;
    while ( i < x.size() ) {
	int deviation = abs(y.at(i) - Linear_fit->Eval(x.at(i)));
	if ( deviation > 100 ) break;
	
	i++;
    }
    return y.at(i);
}



