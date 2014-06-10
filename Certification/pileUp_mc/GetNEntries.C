{

	vector<string> samples;
	samples.push_back("ZJetsToNuNu_50_HT_100_7TeV-madgraph_Summer11.root");
	samples.push_back("ZJetsToNuNu_100_HT_100_7TeV-madgraph_Summer11.root");
	samples.push_back("ZJetsToNuNu_200_HT_inf_7TeV-madgraph_Summer11.root");

	for(unsigned int i=0; i<samples.size(); ++i){
		TFile *files;
		files      = TFile::Open(samples[i].c_str());
		TH1D *histo= (TH1D*) files->Get("pileup");
		TString nentries = TString::Format("%.0f",histo ->GetEntries());
		cout << samples[i] << " nentries: " << nentries<< endl;
	}	


}
