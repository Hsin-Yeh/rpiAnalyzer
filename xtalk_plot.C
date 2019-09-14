map<int,pair < double,double > > CHmap;

void readmap(){
    
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

void xtalk_plot(){
    ifstream f;
    ifstream f_datainput;
    char filename[100];
    f_datainput.open("data_input.txt");
    f_datainput >> filename;
    f.open(filename);

    double xtalk_avg[128];
    double xtalk_intersect[128];
    double xtalk_slope[128];
    double channelID[128];
    while(!f.eof()){
	int chip, ch, icell;
	double avg, posx, posy, intersect, slope;
	f >> chip >> ch >> icell >> avg >> intersect >> slope;
	cout << chip << " " << ch << " " << icell << " " << avg << " " << intersect << " " << slope << endl;
	xtalk_avg[chip*32 + ch/2] = avg;
	xtalk_intersect[chip*32 + ch/2]  = intersect;
	xtalk_slope[chip*32 + ch/2]  = slope;
    }
    for(int ich = 0; ich < 128; ich++){
	channelID[ich] = ich*2;
    }

    TCanvas *c = new TCanvas();
    TGraph *g_avg = new TGraph(128, channelID, xtalk_avg);
    g_avg->SetMarkerStyle(22);
    g_avg->GetXaxis()->SetTitle("channelID");
    g_avg->GetYaxis()->SetTitle(" < EFirstRing / EInj > ");
    g_avg->SetTitle("charge injection xtalk");
    g_avg->Draw("AP");
    c->Update();
    gPad->WaitPrimitive();
    c->SaveAs("output.pdf");

    TH2Poly *poly = new TH2Poly;
    InitTH2Poly(*poly);
    poly->SetMinimum(-0.1);
    for(int ichannel = 0; ichannel < 256; ichannel+=2){
	int ichip = ichannel / 64;
	float X, Y;
	int forCH = ichannel / 2;
	bool NoisyBool = false;
	X = CHmap[forCH].first;
	Y = CHmap[forCH].second;
	poly->Fill(X,Y,xtalk_intersect[ichannel/2]);
    }
    poly->Draw("colztext");
    c->Update();
    gPad->WaitPrimitive();
    
    TGraph *g_intersect = new TGraph(128, channelID, xtalk_intersect);
    g_intersect->SetMarkerStyle(22);
    g_intersect->GetXaxis()->SetTitle("channelID");
    g_intersect->GetYaxis()->SetTitle(" < EFirstRing / EInj > ");
    g_intersect->SetTitle("charge injection xtalk");
    g_intersect->Draw("AP");
    c->Update();
    gPad->WaitPrimitive();
    c->SaveAs("output.pdf");

    TGraph *g_slope = new TGraph(128, channelID, xtalk_slope);
    g_slope->SetMarkerStyle(22);
    g_slope->GetXaxis()->SetTitle("channelID");
    g_slope->GetYaxis()->SetTitle(" slope ");
    g_slope->SetTitle("charge injection xtalk");
    g_slope->Draw("AP");
    c->Update();
    gPad->WaitPrimitive();
    c->SaveAs("output.pdf");

}

