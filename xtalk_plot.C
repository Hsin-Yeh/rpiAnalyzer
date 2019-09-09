void xtalk_plot(){
    ifstream f;
    ifstream f_datainput;
    char filename[100];
    f_datainput.open("data_input.txt");
    f_datainput >> filename;
    f.open(filename);

    double xtalk[128];
    double channelID[128];
    while(!f.eof()){
	int chip, ch;
	double x;
	f >> chip >> ch >> x;
	xtalk[chip*32 + ch/2] = x;
    }
    for(int ich = 0; ich < 128; ich++){
	channelID[ich] = ich*2;
    }

    TCanvas *c = new TCanvas();
    TGraph *g = new TGraph(128, channelID, xtalk);
    g->SetMarkerStyle(22);
    g->Draw("AP");
    c->Update();
    c->SaveAs("output.pdf");
}
