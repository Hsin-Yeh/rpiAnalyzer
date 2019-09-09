#include "makePlots.h"
#include "rootFileIntegrator.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "TImage.h"

//#define DEBUG

//Constructor
makePlots::makePlots(TChain* chain1, TChain* chain2):Chain1(chain1),Chain2(chain2)
{
    cout << "Constructor of makePlot ... \n\n" << endl;
}
//Destructor
makePlots::~makePlots()
{
    //deallocate();
    delete c; // delete canvas first before deleting app
    delete app;
    cout << "\n\n";
    cout << "Destructor of makePlot ... " << endl;
}


///
/// ==================== Init ==================== ///
///
void makePlots::Init( string pedfile, string gainfile, string noisyfile ){
    
    cout << "----------Init start----------" << endl;
    // read init files and initialize
    readmap();
    yamlReader();
    GainFactorReader( gainfile );
    noisyChannelReader( noisyfile );
    read_P_and_N( pedfile );
    //Crosstalk(injCh);
    P.root_logon();

    //moduleNumber = "97";
    // init private parameters
    init_outputFile();
    init_rootBranch();
    init_analysisParameter(); // always after init_rootBranch();
    init_rootDir();
    init_histo();

    // make directories
    char command[100];
    sprintf(command,"mkdir plots/%s",moduleNumber.c_str());
    system(command);
    sprintf(plot_dir,"plots/%s",moduleNumber.c_str());
        
    // init Canvas 
    gROOT->SetBatch("kTRUE");
    app = new TApplication("app",0,0);
    c = new TCanvas();
	
    cout << "----------Init complete----------" << endl << endl;
    
}

///
/// ==================== pedestalPlotter ==================== ///
///
void makePlots::pedestalPlotter(){

    /// Define Parameters 
    int MaxTS = 2; //choose this time sample to be the peak


    /// --------------- Start of pedstalPlotter Loop --------------- ///
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
		hg_sig[ich][sca] = hg[sca][ich];
		lg_sig[ich][sca] = lg[sca][ich];
	    }
	}
	if(subPed_flag) { Pedestal_CM_Subtractor( chip ); } // Ped & CM Subtraction


	// TOT & TOA Fire at pedestal Warning
	for(int ich = 0; ich < NCH; ich++){
	    if ( toa_fall[ich] > 4 || toa_rise[ich] > 4 || tot_slow[ich] > 4 ) {
		cout << "!!!!!! TOA Fire at Pedestal Run !!!!!!!" << endl;
		cout << "event: " << event << " chip: " << chip << " ch: " << ich << " "
		     << tot_slow[ich] << " " << toa_rise[ich] << " " << toa_fall[ich] << endl;
	    }
	}
	
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

	/// Fill Pedestal histograms
	for (int ich = 0; ich < NCH; ich++){
	    for (int sca = 0; sca < NSCA; sca++){
		int ichannel = chip*64 + ich;
		h_hgPedestal[sca][ichannel]->Fill( hg_sig[ich][sca] );
		h_lgPedestal[sca][ichannel]->Fill( lg_sig[ich][sca] );
		h_tot[ichannel]->Fill( tot_slow[ich] );
		h_toarise[ichannel]->Fill( toa_rise[ich] );
		h_toafall[ichannel]->Fill( toa_fall[ich] );
		h_toaf_r[ichannel]->Fill( toa_fall[ich] - toa_rise[ich] );
	    }
	}
	
    }

    /// --------------- End of pedestalPlotter Loop --------------- ///


    /// Fit
    fit_pedestalHisto();

    ///
    /// Plots!!!!!
    ///
    Pedestal_poly();


    outfile->Write();
    outfile->Close();
}


///
/// ==================== sweepPlotter ==================== ///
///
void makePlots::sweepPlotter(){
    
    /// Define Parameters 
    int MaxTS = 2; //choose this time sample to be the peak
    int AverageEvents = 0;

  
    /// --------------- Start Loop of sweepPlotter --------------- ///
    for(int entry = 0; entry < TotalEntries ; ++entry) {
    
	if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }
	Chain1 -> GetEntry(entry);
	Chain2 -> GetEntry(entry);
	dac_ctrl[event] = dacinj;

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
		hg_sig[ich][sca] = hg[sca][ich];
		lg_sig[ich][sca] = lg[sca][ich];
	    }
	}
	
	if(subPed_flag) { Pedestal_CM_Subtractor( chip ); } // Ped & CM Subtraction


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

    
	/// Pedestal histograms
	for (int ich = 0; ich < NCH; ich++){
	    for (int sca = 0; sca < NSCA; sca++){
		int ichannel = chip*64 + ich;
		h_hgPedestal[sca][ichannel]->Fill( hg_sig[ich][sca] );
		h_lgPedestal[sca][ichannel]->Fill( lg_sig[ich][sca] );
		h_tot[ichannel]->Fill( tot_slow[ich] );
		h_toarise[ichannel]->Fill( toa_rise[ich] );
		h_toafall[ichannel]->Fill( toa_fall[ich] );
		h_toaf_r[ichannel]->Fill( toaf_r_allCh[ich][event] );
	    }
	}

    
	/// Injection & Cross Talk Analysis
	for(int ich = 0; ich < NCH; ich++){
	    int channel      = ich + chip*64;
	    double hg = hg_sig[ich][MaxTS_sca];
	    double lg = lg_sig[ich][MaxTS_sca];
	    double tot = tot_slow[ich];
			
	    hg_allCh[channel][event]  = hg;
	    lg_allCh[channel][event]  = lg;
	    tot_allCh[channel][event] = tot;
	    toaf_allCh[channel][event] = toa_fall[ich];
	    toar_allCh[channel][event] = toa_rise[ich];
	    toaf_r_allCh[channel][event] = toa_fall[ich] - toa_rise[ich];
	    /// mip conversion
	    double energy_mip = mipConverter( hg, lg, tot, channel);
	    mip_allCh[channel][event] = energy_mip;
	}

	/// Injection XTalk calculation
	if ( chip == 3 ) {
	    for(int ichannel = 0; ichannel < NCHANNEL; ichannel++){

		int ichip = ichannel / NCH;
		int inj_channel;
		if ( oneChannelInjection_flag )
		    inj_channel = ( injChip * NCH ) + injCh;
		else
		    inj_channel = ( ichip * NCH ) + injCh;

		XTalkCoupling[ichannel][event] = mip_allCh[ichannel][event] / mip_allCh[inj_channel][event];

#ifdef DEBUG 
		cout << "event = " << event << " channel = " << ichannel << " energy = " << mip_allCh[ichannel][event] << " Xtalk = " << XTalkCoupling[ichannel][event] << endl;
#endif
	
		if( event>200 && event<=700 ){
		    XTalkCoupling_Average[ichannel] += XTalkCoupling[ichannel][event];
		    AverageEvents++;
		}
		/// Calulate event ring Energy
		int iring;
		iring = ringPositionFinder( inj_channel, ichannel );
		if( iring > -1 ) {
		    if ( oneChannelInjection_flag )
			mip_Ring_1Chip[iring][event] += mip_allCh[ichannel][event];
		    else
			mip_Ring_4Chip[iring][ichip][event] += mip_allCh[ichannel][event];
		}
	    }

	    /// Calculate XTalkCoupling 
	    if ( oneChannelInjection_flag ) {
		for(int iring = 1; iring < NRings; iring++) {
		    XTalkCoupling_Ring_1Chip[iring][event] = mip_Ring_1Chip[iring][event] / mip_Ring_1Chip[0][event];
		}
		if( event>200 && event<=700 ) {
		    XTalkCoupling_Ring_1Chip_average += XTalkCoupling_Ring_1Chip[1][event];
//#ifdef DEBUG
		    cout << "XTalkCoupling_Ring_1Chip_average = "  << XTalkCoupling_Ring_1Chip_average << endl;
//#endif
		}

	    }
	    else {
		for(int ichip = 0; ichip < NCHIP; ichip++){
		    for(int iring = 1; iring < NRings; iring++) {
			XTalkCoupling_Ring_4Chip[iring][ichip][event] = mip_Ring_4Chip[iring][ichip][event] / mip_Ring_4Chip[0][ichip][event];
		    }
		}
	    }
	}
    }
    /// --------------- End Loop of sweepPlotter --------------- ///

    
    /// --------------- Start of Loop --------------- ///
    for(int entry = 0; entry < Chain2->GetEntries(); entry++){
	Chain2->GetEntry(entry);
		
	for (int ihit = 0; ihit < channelID->size(); ihit++){
	    //cout << "" << endl;
	}
    }

    for (int ichannel = 0; ichannel < NCHANNEL; ichannel++) {
    	XTalkCoupling_Average[ichannel] /= (AverageEvents/NCHANNEL);
	XTalkCoupling_Ring_1Chip_average /= (AverageEvents/NCHANNEL);
    }

    /// Fit
    fit_pedestalHisto();
    
    ///
    /// Plots!!!!!
    ///
    if(oneChannelInjection_flag) { oneChannelInjection_injectionPlots(); }    // hglgtot, mip, xtalk_ring VS dac_ctrl
    else { injectionPlots(); }
    injectionPlots_allCh(); // collect hglgtot plots for all channel 
    Pedestal_poly();    // pedestal 2d poly
    Xtalk_1D();
    XTalk_poly();      // XTalk 2d poly
    toa_plot();

    ///
    /// Output 
    ///
    output_xtalkCoupling();
    
    outfile->Write();
    outfile->Close();

}


///
/// ==================== const_injPlotter ==================== ///
///
///void makePlots::const_injPlotter() {

    
//}


///
/// ==================== cosmicAnalyzer ==================== ///
///
void makePlots::cosmicAnalyzer(){
	
    /// Declare Parameters
    int TotalEntries = Chain1->GetEntries();
    int Nevents = TotalEntries/NCHIP;
    cout << "Total Events = " << Nevents << endl;
    int MaxTS = 2;              //choose this time sample to be the peak
    int mipCount = 0;

    /// Define Histograms
    TH1D *h_mipAllCh = new TH1D("h_mipAllCh","",50,0,400);
  
  
    //==================== Loop over the events ====================
   
    for(int entry = 0; entry < TotalEntries ; ++entry){
	if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }
	Chain1 -> GetEntry(entry);

	/// Timesample
	int TS[NSCA];
	int TS0_sca, MaxTS_sca;
	for(int sca = 0 ; sca < NSCA ; sca++) {
	    TS[sca] = timesamp[sca];
	    if (timesamp[sca] == 0) { TS0_sca = sca ; }
	    if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
	}

	/// hg, lg signal passed to self defined parameters 
	for (int ich = 0; ich < NCH; ich++){
	    for (int sca = 0; sca < NSCA; sca++){
		hg_sig[ich][sca] = hg[sca][ich];
		lg_sig[ich][sca] = lg[sca][ich];
	    }
	}
	if(subPed_flag) { Pedestal_CM_Subtractor( chip ); }  // Ped & CM Subtraction


	/// Calculate hit per chip
	int hit = 0;
	for(int ich = 0; ich < NCH; ich++){
	    if ( mipSigCheck(hg_sig[ich], TS ) ) { hit++; }
	}

	/// Fill histogram
	for (int ich = 0; ich < NCH; ich+=2) {
	    if ( ich + chip*NCH == 44 ) continue;
	    //if ( ich != 18 ) continue;
	    if ( mipSigCheck(hg_sig[ich], TS ) && hit == 1 ) {
		h_mipAllCh->Fill( hg_sig[ich][MaxTS_sca] );
		pulsePlotter( hg_sig[ich], TS , event, chip, ich, -1, -1);
		mipCount++;
	    }
	}
    }

    //... ==================== End of Loop ==================== ...

    cout << endl << "totalEvent# = " << Nevents << " signal# = " << mipCount << endl;
    cout << "efficiency = " << (float)mipCount / Nevents << endl;
 
    outfile->Write();
    outfile->Close();

}


///
/// ==================== Pulse_display ==================== ///
///
void makePlots::Pulse_display( int displayChannel, int pulseDisplay_type, int lowerR, int upperR, int startEv ){

    int Nevents = Chain1->GetEntries();
    cout << "Total Events = " << Nevents << endl;

    /// --------------- Loop Over Events --------------- ///
    for(int ev = 0; ev < Nevents ; ++ev){
	Chain1 -> GetEntry(ev);
	if( event < startEv ) continue;
	int TS[NSCA];
	for(int i = 0 ; i < NSCA ; ++i)  TS[i] = timesamp[i];

	/// Passing hg lg signal to self defined array
	for (int ich = 0; ich < NCH; ich++){
	    for (int sca = 0; sca < NSCA; sca++){
		hg_sig[ich][sca] = hg[sca][ich];
		lg_sig[ich][sca] = lg[sca][ich];
	    }
	}
	if(subPed_flag) { Pedestal_CM_Subtractor( chip ); } // Ped & CM Subtraction

	/// --------------- Pulse Plots --------------- ///
		
	/// Loop Over every Channel
	if ( pulseDisplay_type == 0 && displayChannel == -1) { 
	    for( int ich =0; ich < 64; ich+=2){
		pulsePlotter( hg_sig[ich], TS, event, chip, ich, lowerR, upperR);
	    }
	}
	/// injCh display
	else if ( pulseDisplay_type == 1 ) { 
	    pulsePlotter( lg_sig[injCh], TS, event, chip, injCh, lowerR, upperR);
	}
	/// find signal and display 
	else if ( pulseDisplay_type == 2 ) {  
	    for( int ich =0; ich < 64; ich+=2){
		if ( mipSigCheck( hg_sig[ich], TS ) ) {
		    pulsePlotter( hg_sig[ich], TS, event, chip, ich, lowerR, upperR );
		}
	    }
	}
	/// selected channel display
	else {  
	    int ichip = displayChannel / 64;
	    int ich   = displayChannel % 64;
	    if ( chip != ichip ) continue;
	    pulsePlotter( hg_sig[ich], TS, event, ichip, ich, lowerR, upperR);
	}
    }

}
///
/// ==================== mipConverter ==================== ///
///
double makePlots::mipConverter( double hg_SubPed, double lg_SubPed, double tot , int channel){
  
    double mip;
    int ichip = channel / NCH;
    int ich   = channel % NCH;
  
    if( hg_SubPed < HGTP[ichip][ich]){
	mip = hg_SubPed * ADC2MIP;
    }
    else{
	//if( lg_sig < LGTP[channel]){
	if( lg_SubPed < LGTP_default){
	    mip = ( lg_SubPed * LG2HG_Conversion[ichip][ich] * ADC2MIP);
	}
	else{
	    mip = ( (tot - TOTOffSet[ichip][ich]) * TOT2LG_Conversion[ichip][ich] * LG2HG_Conversion[ichip][ich] * ADC2MIP);
	}
    }
    return mip;
    
}

///
/// ==================== ringPositionFinder ==================== ///
///
int makePlots::ringPositionFinder( int inj_channel, int ichannel){

    int formatCh = ichannel / 2;
    int formatInjCh = inj_channel / 2;
    double X     = CHmap[formatCh].first;
    double Y     = CHmap[formatCh].second;
    double inj_X = CHmap[formatInjCh].first;
    double inj_Y = CHmap[formatInjCh].second;

    int ring = -1;
  
    if ( ichannel % 2 == 0 ){ 
	double dx = X - inj_X;
	double dy = Y - inj_Y;
	double dR = sqrt(dx*dx + dy*dy);
	
	if ( formatCh == formatInjCh ){ ring = 0; }
	else {
	    ring = 0;
	    while(true) {
		ring++;
		if( ring == 5 ) {
		    ring = -1;
		    break; 
		}
		if( dR < 1.12455*ring*1.2 ) break;
	    }
	}
    }

    return ring;
}

///
/// ==================== yamlReader ==================== ///
///
void makePlots::yamlReader(){

    int start = input_fileName.find_last_of("/");
    int end = input_fileName.find("_pedestal.root");
    string f_substr = input_fileName.substr(start+1,end-start-1);
    string rootFileName(f_substr);
    end = input_fileName.find("ana_output");
    f_substr = input_fileName.substr(0,end-1);
    string yamlPath(f_substr);
    char yamlFileName[100];
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

///
/// ==================== GainFactorReader ==================== ///
///
void makePlots::GainFactorReader( string gainfile ){

    string TB_GainFactors("src_txtfile/TPro_fittingoutput.txt");
  
    ifstream GainFile(gainfile);
    string line;
    char tmp[50];
    int ichip, ich;
    if(!GainFile.is_open()){
	cout << "Did not find GainFactor file " << gainfile
	     << ".\n Take this run's GainFactor as default & calculate the Gainfactors " << endl;
	for(ichip = 0; ichip < NCHIP; ichip++){
	    for(ich = 0; ich < NCH; ich++){
		HGTP[ichip][ich] = 1500;
		LG2HG_Conversion[ichip][ich] = 8.5;
		LGTP[ichip][ich] = 900;
		TOT2LG_Conversion[ichip][ich] = 3.8;
		TOTOffSet[ichip][ich] = 180;
	    }
	}
	//Gain_factor_producer();
    }
  
    if(GainFile.is_open()){
	cout << "gainFile = " << gainfile << endl;
	getline(GainFile,line);
	while(!GainFile.eof()){
	    GainFile >> tmp >> tmp >> ichip >> ich >> tmp;
	    GainFile >> HGTP[ichip][ich] >> LG2HG_Conversion[ichip][ich] >> LGTP[ichip][ich] >> TOT2LG_Conversion[ichip][ich] >> TOTOffSet[ichip][ich];
	    getline(GainFile,line);
	}
    }
    cout << endl;
}
///
/// ==================== noisyChannelReader ==================== ///
///
void makePlots::noisyChannelReader( string noisyFileName ) {
    ifstream noisyFile( noisyFileName );
    string line;
    int ichip, ich;
    if(!noisyFile.is_open()){
	cout << "Did not find noisyChannelFile " << noisyFileName
	     << ".\n not setting noisyChannels " << endl;
    }
    else{
	cout << "noisyFile = " << noisyFileName << endl;
	while(!noisyFile.eof()){
	    noisyFile >> ichip >> ich;
	    int ichannel = ichip*NCH + ich;
	    noisyChannel.push_back(ichannel);
	}
	noisyChannel.pop_back();
    }
    cout << "noisy channels = ";
    for(int i = 0; i < noisyChannel.size(); i++){
	cout << noisyChannel.at(i) << " " ;
    }
    cout << endl << endl;
}


///
/// ==================== Pedestal_CM_Subtractor ==================== ///
///
void makePlots::Pedestal_CM_Subtractor( int chip ){
	
    for (int ich = 0; ich < NCH; ich++){
	for (int sca = 0; sca < NSCA; sca++){
	    hg_sig[ich][sca] -= avg_HG[chip][ich][sca];  // Pedestal Subtraction
	    lg_sig[ich][sca] -= avg_LG[chip][ich][sca];
	}
    }

    double *hgCM_sca, *lgCM_sca;
    hgCM_sca = CMCalculator_v2( hg_sig, chip ); // Calculate CM for each sca
    lgCM_sca = CMCalculator_v2( lg_sig, chip );
	

    for (int ich = 0; ich < NCH; ich++){
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
double* makePlots::CMCalculator_v2 ( double **sig_subPed, int chip ) {
    // Calculate CM for each TS
    static double meanChipPedestal[NSCA];

    int scaCount[NSCA];
    for (int sca = 0; sca < NSCA; sca++) {
	meanChipPedestal[sca] = 0;
	scaCount[sca] = 0;
    }

    for (int ich = 0; ich < NCH; ich+=2) {
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
double makePlots::CMCalculator( double **sig_subPed, int *TS ) {
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
/// ==================== mipSegCheck ==================== ///
///
bool makePlots::mipSigCheck( double *sig, int *TS ) {
    int p_noisy_cut  = 2000;
    int n_noisy_cut  = -1000;
    double noSignal_cut = 30;
    bool sig_flag = false;
    double sig_ts[NSCA];
  
    for( int sca = 0; sca < NSCA; sca++) {
	sig_ts[ TS[sca] ] = sig[sca];
	if ( TS[sca] > 1 && TS[sca] < 5 && sig[sca] < p_noisy_cut && sig[sca] > noSignal_cut ) sig_flag = true;
    }
    for( int sca = 0; sca < NSCA; sca++) {
	if ( sig[sca] < n_noisy_cut ) sig_flag = false;
    }
    if ( sig_ts[0] > sig_ts[2] || sig_ts[0] > sig_ts[3] ) sig_flag = false;
    if ( sig_ts[0] < -200 || sig_ts[0] > 150 ) sig_flag = false;
    if ( sig_ts[11] < -200 || sig_ts[11] > 150 ) sig_flag = false;
    if ( sig_ts[2] < 15 ) sig_flag = false;
  
    return sig_flag;
}

///
/// ==================== pulsePlotter ==================== ///
///
void makePlots::pulsePlotter( double *sig, int *TS, int ev, int ichip, int ich, int lowerR, int upperR ) {
	
    // This function plots the input 13 timesamples and show it on the screen 
    double double_TS[NSCA];
    for(int i = 0; i < NSCA; i++) double_TS[i] = (double)TS[i];
    TGraph *gr = new TGraph(NSCA, double_TS, sig );
    char plot_title[50];
    gr->SetMarkerColor(ichip+1);
    gr->SetMarkerStyle(22);
    gr->SetMarkerSize(1.2);
    if ( upperR != -1 ) gr->GetYaxis()->SetRangeUser(lowerR,upperR);
    gr->Draw("AP");
    sprintf(plot_title,"evt %d chip %d channel%d", ev, ichip, ich);
    gr->SetTitle(plot_title);
    gr->GetXaxis()->SetTitle("TS");
    gr->GetYaxis()->SetTitle("ADC");
    c->Update();
    gPad->WaitPrimitive();
  
}

///
/// ==================== read_P_and_N ==================== ///
///
void makePlots::read_P_and_N(string ped_file){

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
/// ==================== readmap ==================== ///
///
void makePlots::readmap(){
    ifstream file("./src_txtfile/CH_map.txt");
    string line;
    int ichip,ich,itype,iformatCH;
    double iposx, iposy;
    while(true){
	getline(file,line);
	if( file.eof() ) break;
	file >> ichip >> ich >> iposx >> iposy >> itype;
	iformatCH = ichip*32 + ich/2;
	CHmap[iformatCH] = make_pair(iposx,iposy);}
    file.close();
    //Since there is no such pad, assign a unreasonable value
    CHmap[2*32+60/2] = make_pair(1000.,1000.);

}

///
/// ==================== InitTH2Poly ==================== ///
///
void makePlots::InitTH2Poly(TH2Poly& poly)
{
    int MAXVERTICES = 6;
    double HexX[MAXVERTICES];
    double HexY[MAXVERTICES];
    int iu,iv,CellXYsize;
    ifstream file("src_txtfile/poly_frame.txt");
    string line;
  
    for(int header = 0; header < 4; ++header )     getline(file,line);
  
    while(true){
	getline(file,line);
	if( file.eof() ) break;
	file >> iu >> iv >> CellXYsize;    
	for(int i = 0; i < CellXYsize ; ++i){
	    getline(file,line);
	    file >> HexX[i] >> HexY[i];
	}
	poly.AddBin(CellXYsize, HexX, HexY);
    }
    file.close();
}



///
/// ==================== Gain_factor_producer ==================== ///
///
void makePlots::Gain_factor_producer(){
  
    //-------------------- Define Parameters --------------------

    int TotalEntries = Chain1->GetEntries();
    int Nevents = TotalEntries/NCHIP;
    int MaxTS = 2;

    int hg_injCh[NCHIP][Nevents], lg_injCh[NCHIP][Nevents], tot_injCh[NCHIP][Nevents];
    int dac_ctrl[Nevents];
    int HGLGfitmax[NCHIP], TOTLGfitmax[NCHIP];
    bool HGTP_flag[NCHIP], LGTP_flag[NCHIP];

    //==================== Set Output File ====================

  


    //==================== Initialize ====================

    for(int ichip = 0; ichip < NCHIP; ichip++){
	HGTP_flag[ichip] = false;
	LGTP_flag[ichip] = false;
    }


    //==================== Loop Over Events ====================

    for(int ev = 0; ev < TotalEntries; ev++){
	if(ev%1000 == 0){ cout << "Now Processing = " << ev << endl;}
	Chain1->GetEntry(ev);
	dac_ctrl[event] = dacinj;
    
	int TS0_sca, MaxTS_sca;
	for(int sca = 0 ; sca < NSCA ; sca++) {
	    if (timesamp[sca] == 0) { TS0_sca = sca ; }
	    if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
	}
    
    
	hg_injCh[chip][event] = ( hg[MaxTS_sca][injCh] - hg[TS0_sca][injCh] );
	lg_injCh[chip][event] = ( lg[MaxTS_sca][injCh] - lg[TS0_sca][injCh] );
	tot_injCh[chip][event] = tot_slow[injCh];


	if(hg_injCh[chip][event] > HGTP[chip][injCh] && HGTP_flag[chip] == false){
	    HGTP_flag[chip] = true;
	    HGLGfitmax[chip] = lg_injCh[chip][event];
	}
	if(lg_injCh[chip][event] > LGTP[chip][injCh] && LGTP_flag[chip] == false){
	    LGTP_flag[chip] = true;
	    TOTLGfitmax[chip] = tot_injCh[chip][event];
	}
    }

  
    //==================== End Loop ====================

    //...

    //==================== Plots ====================

    char pltTit[100];
    string Xtit, Ytit, Opt;
    int MkSty, MkClr, LClr, fitmin, fitmax;
    float MkSize, LWid;
    bool Stat, Wait, SavePlot;
    TCanvas *c1 = new TCanvas();

    TGraph** gh = new TGraph*[NCHIP];
    TGraph** gl = new TGraph*[NCHIP];
    TGraph** gtot = new TGraph*[NCHIP];
    TGraph** LG2HG = new TGraph*[NCHIP];
    TGraph** TOT2LG = new TGraph*[NCHIP];

    for(int ichip = 0; ichip < NCHIP; ichip++){
	gh[ichip] = new TGraph(Nevents,dac_ctrl,hg_injCh[ichip]);
	gl[ichip] = new TGraph(Nevents,dac_ctrl,lg_injCh[ichip]);
	gtot[ichip] = new TGraph(Nevents,dac_ctrl,tot_injCh[ichip]);
	LG2HG[ichip] = new TGraph(Nevents,lg_injCh[ichip],hg_injCh[ichip]);
	TOT2LG[ichip] = new TGraph(Nevents,tot_injCh[ichip],lg_injCh[ichip]);
    
	LG2HG[ichip]->Fit("pol1","","",fitmin = 0,fitmax = HGLGfitmax[ichip]);
	//TOT2LG[ichip]->Fit("pol1","","",fitmin = TOTOffSet,fitmax = TOTLGfitmax[ichip]);
	TOT2LG[ichip]->Fit("pol1","","",fitmin = 200,fitmax = 300);
    
	TF1* Linear_fit_LG2HG = LG2HG[ichip]->GetFunction("pol1");
	TF1* Linear_fit_TOT2LG = TOT2LG[ichip]->GetFunction("pol1");
	LG2HG_Conversion[ichip][injCh] = Linear_fit_LG2HG->GetParameter(1);
	TOT2LG_Conversion[ichip][injCh] = Linear_fit_TOT2LG->GetParameter(1);
	TOTOffSet[ichip][injCh] = -Linear_fit_TOT2LG->GetParameter(0)/Linear_fit_TOT2LG->GetParameter(1);
	cout << LG2HG_Conversion[ichip][injCh] << " " <<  TOT2LG_Conversion[ichip][injCh] << " " << TOTOffSet[ichip][injCh] << endl;
	/*
	  TOT2LG[ichip]->Draw("AP");
	  c1->Update();
	  gPad->WaitPrimitive();
	*/
	sprintf(pltTit,"HG_Chip%d",ichip);
	P.GStd(*gh[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);

	sprintf(pltTit,"LG_Chip%d",ichip);
	P.GStd(*gl[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);

	sprintf(pltTit,"TOT_Chip%d",ichip);
	P.GStd(*gtot[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);
    
	sprintf(pltTit,"LG2HG_Chip%d",ichip);
	P.GStd(*LG2HG[ichip], pltTit, Xtit = "LG", Ytit = "HG", Opt = "AP", Wait = 0, SavePlot = 1);
    
	sprintf(pltTit,"TOT2LG_Chip%d",ichip);
	P.GStd(*TOT2LG[ichip], pltTit, Xtit = "TOT", Ytit = "LG", Opt = "AP", Wait = 0, SavePlot = 1);
    }
}

///
/// ==================== Crosstalk ==================== ///
///
void makePlots::Crosstalk(Int_t CH){

    TCanvas* c1 = new TCanvas();
    for(int ichip = 0; ichip < NCHIP; ichip++){
    
	int cross_num = 6;
	int formatInj_Ch = CH/2+ichip*32;
	float Xdist = 0.974452;
	float Ydist = 0.5626;
	float cross_posx[cross_num], cross_posy[cross_num];
	float X = CHmap[formatInj_Ch].first;
	float Y = CHmap[formatInj_Ch].second;
  
	cross_posx[0] = X - Xdist;
	cross_posy[0] = Y + Ydist;
	cross_posx[1] = X;
	cross_posy[1] = Y + 2*Ydist;
	cross_posx[2] = X + Xdist;
	cross_posy[2] = Y + Ydist;
	cross_posx[3] = X + Xdist;
	cross_posy[3] = Y - Ydist;
	cross_posx[4] = X;
	cross_posy[4] = Y - 2*Ydist;
	cross_posx[5] = X - Xdist;
	cross_posy[5] = Y - Ydist;
		
	cout << "Chip" << ichip << " FirstRing Channels = " ;
	for(int icross = 0; icross < 6; icross++){
	    //cout << cross_posx[i] << " " << cross_posy[i] << endl;
	    int ch = 0;
	    bool good_channel = true;
	    while(abs(CHmap[ch].first-cross_posx[icross]) > 1e-4 || abs(CHmap[ch].second-cross_posy[icross]) > 1e-4){
		if( ch > 127) {
		    good_channel = false;
		    break;
		}
		ch++;	    
		//cout << ch << " "<<abs(CHmap[ch].first - cross_posx[i])<< " " << abs(CHmap[ch].second - cross_posy[i])<< endl;
		/*
		  forCH = chip*32 + cross_ch[i]/2;
		  X = CHmap[forCH].first;
		  Y = CHmap[forCH].second;	      
		  poly->Fill(X,Y,cross_ch[i]);
		*/
      
	    }
	    cross_ch_FirstRing[ichip][icross] = ch * 2;
	    cout << cross_ch_FirstRing[ichip][icross] << " " ;
	}
	cout << endl;
    }
    cout << endl;
    TH2Poly *poly = new TH2Poly;
    InitTH2Poly(*poly);
    /*  
	for(int i=0; i<6; i++){
	poly->Draw("colztext0");
	c1->Update();
	//gPad->WaitPrimitive();
	}
    */
    delete c1;
  
}


void makePlots::init_outputFile() {
    
    /// Set Output Root File
    int start = input_fileName.find_last_of("/");
    int end   = input_fileName.find(".root");
    string outf = input_fileName.substr(start+1,end-start-1);

    integrate_flag = false;
    if ( integrate_flag ) {
	sprintf(title,"root_plot/%s.root",moduleNumber.c_str());
	outfile = new TFile(title,"update");
    }
    else {
	sprintf(title,"root_plot/plot_%s.root",outf.c_str());
	outfile = new TFile(title,"recreate");
    }
    char output_filename[200];
    sprintf(output_filename,"%s",title);
    cout << "output file = " <<  title << endl << endl;

}

void makePlots::init_analysisParameter() {

    XTalkCoupling_Average    = new double[NCHANNEL];
    dac_ctrl                 = new double[Nevents];
    hg_NoisyChannel          = new double[Nevents];
    hg_allCh        = new double*[NCHANNEL];
    lg_allCh        = new double*[NCHANNEL];
    tot_allCh       = new double*[NCHANNEL];
    toaf_allCh      = new double*[NCHANNEL];
    toar_allCh      = new double*[NCHANNEL];
    toaf_r_allCh    = new double*[NCHANNEL];
    mip_allCh       = new double*[NCHANNEL];
    XTalkCoupling   = new double*[NCHANNEL];
    hgFitMean       = new double*[NCHANNEL];
    lgFitMean       = new double*[NCHANNEL];
    hgFitSigma      = new double*[NCHANNEL];
    lgFitSigma      = new double*[NCHANNEL];
    hgFitChisquare  = new double*[NCHANNEL];
    lgFitChisquare  = new double*[NCHANNEL];
    hgMean          = new float*[NCHANNEL];
    lgMean          = new float*[NCHANNEL];
    hgSigma         = new float*[NCHANNEL];
    lgSigma         = new float*[NCHANNEL];
    for(int i = 0; i < NCHANNEL; i++){
	hg_allCh[i]        = new double[Nevents];
	lg_allCh[i]        = new double[Nevents];
	tot_allCh[i]       = new double[Nevents];
	toaf_allCh[i]      = new double[Nevents];
	toar_allCh[i]      = new double[Nevents];
	toaf_r_allCh[i]    = new double[Nevents];
	mip_allCh[i]       = new double[Nevents];
	XTalkCoupling[i]   = new double[Nevents];
	hgFitMean[i]       = new double[NSCA];
	lgFitMean[i]       = new double[NSCA];
	hgFitSigma[i]      = new double[NSCA];
	lgFitSigma[i]      = new double[NSCA];
	hgFitChisquare[i]  = new double[NSCA];
	lgFitChisquare[i]  = new double[NSCA];
	hgMean[i]          = new float[NSCA];
	lgMean[i]          = new float[NSCA];
	hgSigma[i]         = new float[NSCA];
	lgSigma[i]         = new float[NSCA];
    }
    hg_SubPed = new double*[NSCA];
    lg_SubPed = new double*[NSCA];
    for(int i = 0; i < NSCA; i++){
	hg_SubPed[i] = new double[NCH];
	lg_SubPed[i] = new double[NCH];
    }
    hg_SubPedCM = new double*[NCH];
    lg_SubPedCM = new double*[NCH];
    for(int i = 0; i < NCH; i++){
	hg_SubPedCM[i] = new double[NSCA];
	lg_SubPedCM[i] = new double[NSCA];
    }
    mip_Ring_4Chip           = new double**[NRings];
    XTalkCoupling_Ring_4Chip = new double**[NRings];
    for(int i = 0; i < NRings; i++){
	mip_Ring_4Chip[i]           = new double*[NCHIP];
	XTalkCoupling_Ring_4Chip[i] = new double*[NCHIP];
	for(int j = 0; j < NCHIP; j++){
	    mip_Ring_4Chip[i][j]           = new double[Nevents];
	    XTalkCoupling_Ring_4Chip[i][j] = new double[Nevents];
	}
    }
    mip_Ring_1Chip           = new double*[NRings];
    XTalkCoupling_Ring_1Chip = new double*[NRings];
    for(int i = 0; i < NRings; i++){
	mip_Ring_1Chip[i]           = new double[Nevents];
	XTalkCoupling_Ring_1Chip[i] = new double[Nevents];
    }

    hg_sig = new double*[NCH];
    lg_sig = new double*[NCH];
    for(int i = 0; i < NCH; i++){
	hg_sig[i] = new double[NSCA];
	lg_sig[i] = new double[NSCA];
    }

    
    for( int ichannel = 0; ichannel < NCHANNEL; ichannel++) {
	XTalkCoupling_Average[ichannel] = 0;
    }
    XTalkCoupling_Ring_1Chip_average = 0;
    
}

void makePlots::init_rootBranch() {
    
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
    cout << "Total Events = " << Nevents << endl;
    
}

void makePlots::init_rootDir() {
    
    /// Declare TDirectories
    if ( acquisitionType == "standard" ) {
	sprintf(title,"Pedestal_Run");
	cdinjCh = outfile->mkdir(title);
    }
    else if ( acquisitionType == "sweep" ) {
	sprintf(title,"InjCh%d",injCh);
	cdinjCh = outfile->mkdir(title);
	cdallCh = cdinjCh->mkdir("allCh_hglgtot");
	cdinj = cdinjCh->mkdir("injection_analysis_plots");
    }

    cdPedestal_histo = cdinjCh->mkdir("Pedestal_histo");
    cdtot = cdPedestal_histo->mkdir("tot_toa_histo");
    cdhglgPedestal = cdPedestal_histo->mkdir("hglg_histo");
    cdPedestal_poly = cdinjCh->mkdir("Pedestal_poly");
    cdhgPedestal = cdPedestal_poly->mkdir("hgPedestal");
    cdlgPedestal = cdPedestal_poly->mkdir("lgPedestal");
    cdhgNoise = cdPedestal_poly->mkdir("hgNoise");
    cdlgNoise = cdPedestal_poly->mkdir("lgNoise");
    cdhgChi = cdPedestal_poly->mkdir("hgChisquare");
    cdlgChi = cdPedestal_poly->mkdir("lgChisquare");
    
}

void makePlots::init_histo() {
    
    /// Initialize Histograms
    for(int ichannel = 0; ichannel < NCHANNEL; ichannel++) {
	char h_title[50];
	for(int sca = 0; sca < NSCA; ++sca) {
	    cdhglgPedestal->cd();
	    sprintf(h_title,"h_hgPedestal_Ch%d_SCA%d", ichannel, sca);
	    h_hgPedestal[sca][ichannel] = new TH1D(h_title,h_title,500,-500,500); 
	    sprintf(h_title,"h_lgPedestal_Ch%d_SCA%d", ichannel, sca);
	    h_lgPedestal[sca][ichannel] = new TH1D(h_title,h_title,500,-500,500);
	}
	cdtot->cd();
	sprintf(h_title,"h_tot_Ch%d", ichannel);
	h_tot[ichannel] = new TH1D(h_title,h_title,150,0,3000);
	sprintf(h_title,"h_toarise_Ch%d", ichannel);
	h_toarise[ichannel] = new TH1D(h_title,h_title,150,0,3000);
	sprintf(h_title,"h_toafall_Ch%d", ichannel);
	h_toafall[ichannel] = new TH1D(h_title,h_title,150,0,3000);
	sprintf(h_title,"h_toaf_r_Ch%d", ichannel);
	h_toaf_r[ichannel] = new TH1D(h_title,h_title,200,700,1200);
    }
}


void makePlots::fit_pedestalHisto() {
    
        /// Fit 
    for(int ichannel = 0; ichannel < NCHANNEL; ichannel++){
	
	for (int sca = 0; sca < NSCA; sca++){
	    if (ichannel%2 == 1) continue;
	    h_hgPedestal[sca][ichannel]->Fit("gaus","Q");
	    hgFitMean [ichannel][sca] = h_hgPedestal[sca][ichannel]->GetFunction("gaus")->GetParameter(1);
	    hgFitSigma[ichannel][sca] = h_hgPedestal[sca][ichannel]->GetFunction("gaus")->GetParameter(2);
	    hgFitChisquare[ichannel][sca] = h_hgPedestal[sca][ichannel]->GetFunction("gaus")->GetChisquare();
	    
	    h_lgPedestal[sca][ichannel]->Fit("gaus","Q");
	    lgFitMean [ichannel][sca] = h_lgPedestal[sca][ichannel]->GetFunction("gaus")->GetParameter(1);
	    lgFitSigma[ichannel][sca] = h_lgPedestal[sca][ichannel]->GetFunction("gaus")->GetParameter(2);
	    lgFitChisquare[ichannel][sca] = h_lgPedestal[sca][ichannel]->GetFunction("gaus")->GetChisquare();

#ifdef DEBUG
	    //if ( hgFitMean[ichannel][sca] > 200 ) 
		cout << " ichannel " << ichannel << " sca " << sca
		     << " Mean " << hgFitMean[ichannel][sca] << " Sigma " << hgFitSigma[ichannel][sca] << " Chi " << hgFitChisquare[ichannel][sca]
		     << " Mean " << lgFitMean[ichannel][sca] << " Sigma " << lgFitSigma[ichannel][sca] << " Chi " << lgFitChisquare[ichannel][sca]
		     << endl;
		    //<< " Chi " << hgFitChisquare[ichannel][sca]<< endl;
#endif
      

	}
    }
    
}


void makePlots::Pedestal_poly() {

    /// 2D Pedestal after CM
    TH2Poly *polyhgPed[NSCA];
    TH2Poly *polylgPed[NSCA];
    TH2Poly *polyhgErr[NSCA];
    TH2Poly *polylgErr[NSCA];
    TH2Poly *polyhgChi[NSCA];
    TH2Poly *polylgChi[NSCA];
    for (int sca = 0; sca < NSCA; sca++){
	cdhgPedestal->cd();
	polyhgPed[sca] = new TH2Poly();
	InitTH2Poly(*polyhgPed[sca]);
	sprintf(title,"hgPedestal_SCA%d",sca);
	polyhgPed[sca]->SetTitle(title);
	polyhgPed[sca]->SetName(title);
	polyhgPed[sca]->SetMarkerSize(1);

	cdhgNoise->cd();
	polyhgErr[sca] = new TH2Poly();
	InitTH2Poly(*polyhgErr[sca]);
	sprintf(title,"hgSigma_SCA%d",sca);
	polyhgErr[sca]->SetTitle(title);
	polyhgErr[sca]->SetName(title);
	polyhgErr[sca]->SetMarkerSize(1);

	cdhgChi->cd();
	polyhgChi[sca] = new TH2Poly();
	InitTH2Poly(*polyhgChi[sca]);
	sprintf(title,"hgChisquare_SCA%d",sca);
	polyhgChi[sca]->SetTitle(title);
	polyhgChi[sca]->SetName(title);
	polyhgChi[sca]->SetMarkerSize(1);
	
	cdlgPedestal->cd();
	polylgPed[sca] = new TH2Poly();
	InitTH2Poly(*polylgPed[sca]);
	sprintf(title,"lgPedestal_SCA%d",sca);
	polylgPed[sca]->SetTitle(title);
	polylgPed[sca]->SetName(title);
	polylgPed[sca]->SetMarkerSize(1);

	cdlgNoise->cd();
	polylgErr[sca] = new TH2Poly();
	InitTH2Poly(*polylgErr[sca]);
	sprintf(title,"lgSigma_SCA%d",sca);
	polylgErr[sca]->SetTitle(title);
	polylgErr[sca]->SetName(title);
	polylgErr[sca]->SetMarkerSize(1);

	cdlgChi->cd();
	polylgChi[sca] = new TH2Poly();
	InitTH2Poly(*polylgChi[sca]);
	sprintf(title,"lgChisquare_SCA%d",sca);
	polylgChi[sca]->SetTitle(title);
	polylgChi[sca]->SetName(title);
	polylgChi[sca]->SetMarkerSize(1);

    }

    for(int sca = 0; sca < NSCA; sca++){
	for(int ichannel = 0; ichannel < NCHANNEL; ichannel+=2){
	    float X, Y;
	    int forCH = ichannel / 2;
	    X = CHmap[forCH].first;
	    Y = CHmap[forCH].second;
	    polyhgPed[sca]->Fill( X, Y, (int)hgFitMean [ichannel][sca]);
	    polyhgErr[sca]->Fill( X, Y, (int)hgFitSigma[ichannel][sca]);
	    polyhgChi[sca]->Fill( X, Y, (int)lgFitChisquare[ichannel][sca]);
	    polylgPed[sca]->Fill( X, Y, (int)lgFitMean [ichannel][sca]);
	    polylgErr[sca]->Fill( X, Y, (int)lgFitSigma[ichannel][sca]);
	    polylgChi[sca]->Fill( X, Y, (int)lgFitChisquare[ichannel][sca]);
	}
	cdhgPedestal->cd();
	polyhgPed[sca]->Write();
	cdlgPedestal->cd();
	polylgPed[sca]->Write();
	cdhgNoise->cd();
	polyhgErr[sca]->Write();
	cdlgNoise->cd();
	polylgErr[sca]->Write();
	cdhgChi->cd();
	polyhgChi[sca]->Write();
	cdlgChi->cd();
	polylgChi[sca]->Write();
    }

}


void makePlots::injectionPlots(){

    cdinj->cd();
    // hglgtot vs dac
    // mip vs dac
    // XTalk_ring vs daq
    TMultiGraph *multig_XTalkCoupling_ring = new TMultiGraph();

    for(int ichip = 0; ichip < NCHIP; ichip++){
	int inj_channel = (ichip*64) + injCh;
	
	TGraph* ginjCh_hg  = new TGraph( Nevents, dac_ctrl, hg_allCh[inj_channel] );
	sprintf(title,"hg");
	ginjCh_hg->SetTitle(title);
	ginjCh_hg->SetName(title);
	ginjCh_hg->SetMarkerColor(P.Color(0));
	TGraph* ginjCh_lg  = new TGraph( Nevents, dac_ctrl, lg_allCh[inj_channel] );
	sprintf(title,"lg");
	ginjCh_lg->SetTitle(title);
	ginjCh_lg->SetName(title);
	ginjCh_lg->SetMarkerColor(P.Color(1));
	TGraph* ginjCh_tot = new TGraph( Nevents, dac_ctrl, tot_allCh[inj_channel] );
	sprintf(title,"tot");
	ginjCh_tot->SetTitle(title);
	ginjCh_tot->SetName(title);
	ginjCh_tot->SetMarkerColor(P.Color(2));
	TGraph* ginjCh_mip = new TGraph( Nevents, dac_ctrl, mip_allCh[inj_channel] );
	sprintf(title,"mip_InjCh%d_chip%d", injCh, ichip);
	ginjCh_mip->SetTitle(title);
	ginjCh_mip->SetName(title);
	ginjCh_mip->Write();
	
	TMultiGraph *multig_InjCh_hltot = new TMultiGraph();
	multig_InjCh_hltot->Add(ginjCh_hg);
	multig_InjCh_hltot->Add(ginjCh_lg);
	multig_InjCh_hltot->Add(ginjCh_tot);
	sprintf(title,"hglgtot_InjCh%d_chip%d", injCh, ichip);
	multig_InjCh_hltot->SetTitle(title);
	multig_InjCh_hltot->SetName(title);
	multig_InjCh_hltot->Write();


	/// Xtalk vs dac_ctrl
	for(int iring = 1; iring < 2; iring++){
	    TGraph* gXTalkCoupling = new TGraph(Nevents, mip_allCh[inj_channel], XTalkCoupling_Ring_4Chip[iring][ichip] );
	    sprintf(title,"chip%d", ichip);
	    gXTalkCoupling->SetTitle(title);
	    gXTalkCoupling->SetName(title);
	    gXTalkCoupling->SetMarkerColor(P.Color(ichip));
	    gXTalkCoupling->SetLineWidth(0);
	    gXTalkCoupling->SetFillColor(0);
	    multig_XTalkCoupling_ring->Add(gXTalkCoupling);
	}
    }
    sprintf(title,"EfirstRing/EInj InjCh%d", injCh);
    multig_XTalkCoupling_ring->SetTitle(title);
    multig_XTalkCoupling_ring->SetName(title);
    multig_XTalkCoupling_ring->Draw("AP");
    c->Update();
    c->BuildLegend(0.75,0.2,0.85,0.35);
    multig_XTalkCoupling_ring->GetXaxis()->SetTitle("Injected Charge [DAC]");
    multig_XTalkCoupling_ring->GetYaxis()->SetTitle("EfirstRing / EInj");
    multig_XTalkCoupling_ring->GetYaxis()->SetTitleOffset(1.2);
    multig_XTalkCoupling_ring->GetYaxis()->SetRangeUser(-0.1,0.1);
    multig_XTalkCoupling_ring->GetXaxis()->SetRangeUser(200,1000);
    multig_XTalkCoupling_ring->Write();
    sprintf(title,"%s/XtalkCoupling_InjCh%d.png", plot_dir, injCh);
    c->SaveAs(title);
    sprintf(title,"%s/XtalkCoupling_InjCh%d.pdf", plot_dir, injCh);
    c->SaveAs(title);

}


void makePlots::oneChannelInjection_injectionPlots(){

    cdinj->cd();
    // hglgtot vs dac
    // mip vs dac
    // XTalk_ring vs daq 
    
    int inj_channel = injCh + (injChip * NCH);
    TGraph* ginjCh_hg  = new TGraph( Nevents, dac_ctrl, hg_allCh[inj_channel] );
    sprintf(title,"hg");
    ginjCh_hg->SetTitle(title);
    ginjCh_hg->SetName(title);
    ginjCh_hg->SetMarkerColor(P.Color(0));
    TGraph* ginjCh_lg  = new TGraph( Nevents, dac_ctrl, lg_allCh[inj_channel] );
    sprintf(title,"lg");
    ginjCh_lg->SetTitle(title);
    ginjCh_lg->SetName(title);
    ginjCh_lg->SetMarkerColor(P.Color(1));
    TGraph* ginjCh_tot = new TGraph( Nevents, dac_ctrl, tot_allCh[inj_channel] );
    sprintf(title,"tot");
    ginjCh_tot->SetTitle(title);
    ginjCh_tot->SetName(title);
    ginjCh_tot->SetMarkerColor(P.Color(2));
    TMultiGraph *multig_InjCh_hltot = new TMultiGraph();
    multig_InjCh_hltot->Add(ginjCh_hg);
    multig_InjCh_hltot->Add(ginjCh_lg);
    multig_InjCh_hltot->Add(ginjCh_tot);
    sprintf(title,"hglgtot_InjCh%d_chip%d", injCh, injChip);
    multig_InjCh_hltot->SetTitle(title);
    multig_InjCh_hltot->SetName(title);
    multig_InjCh_hltot->Write();

    TGraph* ginjCh_mip = new TGraph( Nevents, dac_ctrl, mip_allCh[inj_channel] );
    sprintf(title,"mip_InjCh%d_chip%d", injCh, injChip);
    ginjCh_mip->SetTitle(title);
    ginjCh_mip->SetName(title);
    ginjCh_mip->Write();

    TMultiGraph *multig_XTalkCoupling_ring = new TMultiGraph();
    for(int iring = 1; iring < 2; iring++){
	TGraph* gXTalkCoupling = new TGraph(Nevents, mip_allCh[inj_channel], XTalkCoupling_Ring_1Chip[iring]);
	sprintf(title,"ring %d", iring);
	gXTalkCoupling->SetTitle(title);
	gXTalkCoupling->SetName(title);
	gXTalkCoupling->SetMarkerColor(P.Color(iring-1));
	gXTalkCoupling->SetLineWidth(0);
	gXTalkCoupling->SetFillColor(0);
	multig_XTalkCoupling_ring->Add(gXTalkCoupling);
    }
    sprintf(title,"XtalkCoupling_InjCh%d_chip%d", injCh, injChip);
    multig_XTalkCoupling_ring->SetTitle(title);
    multig_XTalkCoupling_ring->SetName(title);
    multig_XTalkCoupling_ring->Draw("AP");
    multig_XTalkCoupling_ring->GetYaxis()->SetRangeUser(-0.01,0.3);
    multig_XTalkCoupling_ring->Write();
    
}



void makePlots::injectionPlots_allCh() {

    // Collect hglgtot plots for all channels ( even they are not the injected channel )
    cdallCh->cd();
    for(int ichannel = 0; ichannel < NCHANNEL; ichannel++){
	int ichip = ichannel / 64;
	int inj_channel = (ichip*64) + injCh;
	
	TGraph* ginjCh_hg  = new TGraph( Nevents, dac_ctrl, hg_allCh[ichannel] );
	sprintf(title,"hg");
	ginjCh_hg->SetTitle(title);
	ginjCh_hg->SetName(title);
	ginjCh_hg->SetMarkerColor(P.Color(0));
	TGraph* ginjCh_lg  = new TGraph( Nevents, dac_ctrl, lg_allCh[ichannel] );
	sprintf(title,"lg");
	ginjCh_lg->SetTitle(title);
	ginjCh_lg->SetName(title);
	ginjCh_lg->SetMarkerColor(P.Color(1));
	TGraph* ginjCh_tot = new TGraph( Nevents, dac_ctrl, tot_allCh[ichannel] );
	sprintf(title,"tot");
	ginjCh_tot->SetTitle(title);
	ginjCh_tot->SetName(title);
	ginjCh_tot->SetMarkerColor(P.Color(2));
	TGraph* ginjCh_mip = new TGraph( Nevents, dac_ctrl, mip_allCh[ichannel] );
	sprintf(title,"mip_Ch%d", ichannel);
	ginjCh_mip->SetTitle(title);
	ginjCh_mip->SetName(title);
	ginjCh_mip->Write();
	
	TMultiGraph *multig_InjCh_hltot = new TMultiGraph();
	multig_InjCh_hltot->Add(ginjCh_hg);
	multig_InjCh_hltot->Add(ginjCh_lg);
	multig_InjCh_hltot->Add(ginjCh_tot);
	sprintf(title,"hglgtot_Ch%d", ichannel);
	multig_InjCh_hltot->SetTitle(title);
	multig_InjCh_hltot->SetName(title);
	multig_InjCh_hltot->Write();

	TGraph* gXTalkCoupling = new TGraph(Nevents, mip_allCh[inj_channel], XTalkCoupling[ichannel] );
	sprintf(title,"xtalk_Ch%d", ichannel);
	gXTalkCoupling->SetTitle(title);
	gXTalkCoupling->SetName(title);
	gXTalkCoupling->Write();
    }

}

void makePlots::XTalk_poly() {

    cdinj->cd();
    /// 2D Average Xtalk 
    int NNoisy = 8;
    int NoisyChannel[8] = {248,186,214,120,126,42,254,190};

    TH2Poly *poly = new TH2Poly;
    InitTH2Poly(*poly);
    poly->SetMinimum(-0.1);
    for(int ichannel = 0; ichannel < NCHANNEL; ichannel+=2){
	int ichip = ichannel / NCH;
	float X, Y;
	int forCH = ichannel / 2;
	bool NoisyBool = false;
	X = CHmap[forCH].first;
	Y = CHmap[forCH].second;
	if(ichannel%64 == injCh){
	    if(oneChannelInjection_flag==true && ichip==injChip) { poly->Fill(X,Y,-2); }
	    else if (oneChannelInjection_flag==true && ichip!=injChip) { poly->Fill(X,Y,XTalkCoupling_Average[ichannel]); }
	    else { poly->Fill(X,Y,-2); }
	}
	else {
	    if(!NoisyBool){
		poly->Fill(X,Y,XTalkCoupling_Average[ichannel]);
		//poly->Fill(X,Y,forCH);
	    }
	}
    }
    sprintf(title,"<E / EInj> InjCh%d", injCh);
    gStyle->SetPaintTextFormat("2.3f");
    poly->SetTitle(title);
    poly->SetName(title);
    poly->SetMarkerSize(1);
    poly->SetMinimum(-0.02);
    poly->Draw("colztext");
    latex.SetTextSize(0.05);
    latex.SetTextAlign(13);  //align at top
    latex.SetTextSize(0.03);
    latex.DrawLatex(-7,-7,"*White pads are the charge injection pads");
    //latex.DrawLatex(-7,-7.5,"*Each pad is filled with its energy divided by the enrgy of the injection channel energy on the same chip");
    c->Update();
    sprintf(title,"%s/XtalkCoupling_Poly_InjCh%d.png", plot_dir, injCh);
    c->SaveAs(title);
    c->Write();
}

void makePlots::Xtalk_1D(){

    cdinj->cd();
    double labelCnct[NCHANNEL/2], labelUnCnct[NCHANNEL/2];
    double XTalkCoupling_Cnct[NCHANNEL/2], XTalkCoupling_UnCnct[NCHANNEL/2];
    for( int ichannel = 0; ichannel < NCHANNEL; ichannel++ ) {
	if ( ichannel%2 == 1 ) {
	    labelCnct[ichannel/2] = ichannel;
	    XTalkCoupling_UnCnct[ichannel/2] = XTalkCoupling_Average[ichannel];
	}
	else {
	    labelUnCnct[ichannel/2] = ichannel;
	    XTalkCoupling_Cnct[ichannel/2] = XTalkCoupling_Average[ichannel];
	}
    }
    TGraph* gXTalkCoupling_Cnct   = new TGraph(NCHANNEL/2, labelCnct , XTalkCoupling_Cnct);
    sprintf(title,"InjCh%d_XtalkCoupling_Connected_Channel", injCh);
    gXTalkCoupling_Cnct->SetTitle(title);
    gXTalkCoupling_Cnct->SetName(title);
    gXTalkCoupling_Cnct->GetXaxis()->SetTitle("ChannelId");
    gXTalkCoupling_Cnct->GetYaxis()->SetTitle("E / EInj");
    gXTalkCoupling_Cnct->Write();
    sprintf(title,"%s/InjCh%d_XtalkCoupling_ConnectedChannels.png",plot_dir,injCh);
    gXTalkCoupling_Cnct->Draw("AP");
    c->Update();
    c->SaveAs(title);
    TGraph* gXTalkCoupling_UnCnct = new TGraph(NCHANNEL/2, labelUnCnct , XTalkCoupling_UnCnct);
    sprintf(title,"InjCh%d_XtalkCoupling_UnConnected_Channel", injCh);
    gXTalkCoupling_UnCnct->SetTitle(title);
    gXTalkCoupling_UnCnct->SetName(title);
    gXTalkCoupling_UnCnct->Write();

}

void makePlots::toa_plot(){

    // maximum bin
    cdinj->cd();
    int toaf_r_max[NCHANNEL];
    int ChannelID[NCHANNEL];
    for (int ichannel = 0; ichannel < NCHANNEL; ichannel++){
	ChannelID[ichannel] = ichannel;
	toaf_r_max[ichannel] = h_toaf_r[ichannel]->GetMaximumBin();
    }
    TGraph* gtoaf_r = new TGraph(NCHANNEL, ChannelID, toaf_r_max);
    sprintf(title,"InjCh%d_toaFminR", injCh);
    gtoaf_r->SetTitle(title);
    gtoaf_r->SetName(title);
    gtoaf_r->SetLineWidth(0);
    gtoaf_r->SetLineColor(0);
    gtoaf_r->Write();

    // vs dac_ctrl
    TMultiGraph *multig_InjCh_toaFminR = new TMultiGraph();
    TMultiGraph *multig_InjCh_toaF = new TMultiGraph();
    TMultiGraph *multig_InjCh_toaR = new TMultiGraph();
    for(int ichip = 0; ichip < NCHIP; ichip++){
	int inj_channel = (ichip*64) + injCh;
	
	TGraph* ginjCh_toaF  = new TGraph( Nevents, dac_ctrl, toaf_allCh[inj_channel] );
	sprintf(title,"chip%d",ichip);
	ginjCh_toaF->SetTitle(title);
	ginjCh_toaF->SetName(title);
	ginjCh_toaF->SetMarkerColor(P.Color(ichip));
	ginjCh_toaF->SetLineWidth(0);
	ginjCh_toaF->SetFillColor(0);

	TGraph* ginjCh_toaR  = new TGraph( Nevents, dac_ctrl, toar_allCh[inj_channel] );
	sprintf(title,"chip%d",ichip);
	ginjCh_toaR->SetTitle(title);
	ginjCh_toaR->SetName(title);
	ginjCh_toaR->SetMarkerColor(P.Color(ichip));
	ginjCh_toaR->SetLineWidth(0);
	ginjCh_toaR->SetFillColor(0);

	TGraph* ginjCh_toaFminR  = new TGraph( Nevents, dac_ctrl, toaf_r_allCh[inj_channel] );
	sprintf(title,"chip%d",ichip);
	ginjCh_toaFminR->SetTitle(title);
	ginjCh_toaFminR->SetName(title);
	ginjCh_toaFminR->SetMarkerColor(P.Color(ichip));
	ginjCh_toaFminR->SetLineWidth(0);
	ginjCh_toaFminR->SetFillColor(0);

	multig_InjCh_toaF->Add(ginjCh_toaF);
	multig_InjCh_toaR->Add(ginjCh_toaR);
	multig_InjCh_toaFminR->Add(ginjCh_toaFminR);
    }
    sprintf(title,"toaFminR_InjCh%d", injCh);
    multig_InjCh_toaFminR->SetTitle(title);
    multig_InjCh_toaFminR->SetName(title);
    multig_InjCh_toaFminR->Write();
    sprintf(title,"toaF_InjCh%d", injCh);
    multig_InjCh_toaF->SetTitle(title);
    multig_InjCh_toaF->SetName(title);
    multig_InjCh_toaF->Write();
    sprintf(title,"toaR_InjCh%d", injCh);
    multig_InjCh_toaR->SetTitle(title);
    multig_InjCh_toaR->SetName(title);
    multig_InjCh_toaR->Write();
}


void makePlots::output_xtalkCoupling() {
    ofstream f;
    sprintf(title,"xtalkCoupling_%s.txt",moduleNumber.c_str());
    f.open(title, ios::out | ios::app );
    cout << "write " << title << endl;
    cout << XTalkCoupling_Ring_1Chip_average << endl;
    f << injChip << " " << injCh << " " << XTalkCoupling_Ring_1Chip_average << endl;
}


void makePlots::deallocate() {
    
    /// deallocate
    /*
    for (int i = 0; i < NRings; i++){
	for (int j = 0; j < NCHIP; j++) {
	    delete[] mip_Ring_4Chip[i][j];
	    delete[] XTalkCoupling_Ring_4Chip[i][j];
	}
	delete[] mip_Ring_4Chip[i];
	delete[] XTalkCoupling_Ring_4Chip[i];
    }
    delete[] mip_Ring_4Chip;
    delete[] XTalkCoupling_Ring_4Chip;
    */
    

    delete XTalkCoupling_Average;   
    delete dac_ctrl;                
    delete hg_NoisyChannel;

}

