map<int,pair < double,double > > CHmap;

void readmap() {
    
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

void InitTH2Poly(TH2Poly& poly)
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



void xtalk_labdata_compare() {

    readmap();
    ifstream infile;
    infile.open("xtalkCoupling_124oneChannelInjection_all.txt");


    int chip, ch, pad;
    double xtalk_average[128], xtalk_intersect[128], xtalk_slope[128];
    vector<double> xtalk_intersect_vector;
    double tmp, intersect;
    for ( int ichannel = 0; ichannel < 128; ichannel++ ) {
        xtalk_average[ichannel] = 0;
        xtalk_intersect[ichannel] = 0;
        xtalk_slope[ichannel] = 0;
    }
    TH1D *hist = new TH1D("hist","",30,0,0.03);
    hist->GetXaxis()->SetTitle("E_oneFirstRingChannel / E_injectionChannel");
   
    while ( !infile.eof() ) {
        // infile >> chip >> ch >> pad >> xtalk_average[(chip*32) + (ch/2)] >> xtalk_intersect[(chip*32) + (ch/2)] >> xtalk_slope[(chip*32) + (ch/2)];
         infile >> chip >> ch >> pad >> tmp >> intersect >> tmp;
         xtalk_intersect_vector.push_back(intersect);
         hist->Fill(intersect);
         // std::cout << chip << " " << ch << " " << " " << pad << " " << " " << intersect << " " << std::endl;
    }

    std::cout << "hi" << std::endl;
    TH2Poly *poly = new TH2Poly;
    InitTH2Poly(*poly);

    for(int ichannel = 0; ichannel < 128; ichannel++){
        if(xtalk_average[ichannel]==0) {continue;}
        float X, Y;
        int forCH = ichannel;
        X = CHmap[forCH].first;
        Y = CHmap[forCH].second;
        // poly->Fill( X, Y, xtalk_intersect[ichannel] );
        // hist->Fill(xtalk_intersect[ichannel]);
    }

    hist->Fit("gaus");
    // Plot
    TCanvas* c1 = new TCanvas();
    poly->Draw("colztext");
    c1->Update();
    gPad->WaitPrimitive();
    hist->Draw();
    c1->Update();
    gPad->WaitPrimitive();
    
}

