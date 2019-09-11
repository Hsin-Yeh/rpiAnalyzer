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
	f >> chip >> ch >> posx >> posy >> avg >> intersect >> slope;
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

