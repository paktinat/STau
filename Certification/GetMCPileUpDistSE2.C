{
	string base       = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/susy/ntuples/mc/V02-05-00/";
	const int  nfiles = 451;
	//string basestring = "/shome/haweber/MT2Analysis/Code/Macros/";
	//char *textfile; = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/LM6_UpgradeStdGeom2.dat";
	//char *textfile = "/shome/haweber/MT2Analysis/Code/Macros/test.dat";
	//char *textfile;
	//vector<char*> text;
	vector<string> samples;

	string basestring = "/shome/haweber/MT2Analysis/Code/Certification/rootfileslocation/";
	string endstring  = "/shome/haweber/MT2Analysis/Code/Certification/pileUp_mc/";
	string filebase   = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/"
	//samples.push_back("DoubleMu-Run2011A-PromptReco-v4-AOD_RETRANSFER");
	samples.push_back("LM6_UpgradeStdGeom2"); //V02-05-00

	for(unsigned int i=0; i<samples.size(); ++i){
		TChain *t=new TChain("analyze/Analysis");
		//string samplepath=base+samples[i];
		//t->Add(samplepath.c_str());
		//TChain *u = new TChain("analyze");


		//TFile *file;
		//TH1I  *histo;
		//TH1I  *pileupI = new TH1I("pileupI", "pileupI", 40, 0, 40);
		TH1D  *pileup = new TH1D("pileup", "pileup", 51, -0.5, 50.5);
		//TH1D  *pileup2 = new TH1D("pileup2", "pileup2", 51, -0.5, 50.5);

		//TTree * ntree;

		
		cout << " processing sample " << samples[i] << " -----------------" << endl;
			TString rootFile;
			//ifstream is(textfile);
			ifstream is((basestring + samples[i] + string(".dat")).c_str());
			while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
				if(rootFile[0] == '#') continue;
				//theChain->Add(rootFile);
				printf(" Adding file: %s\n", rootFile.Data());
				string currfile = rootFile.Data();
				currfile = 
				//file = TFile::Open(currfile.c_str());
				//histo= (TH1I*) file->Get("analyze/PileUpStats");
				//pileupI->Add(histo);
				//if(t->GetNtrees()==3) continue;
				t->AddFile(currfile.c_str());
				//cout << "entries " << t->GetEntries() << endl;
				//cout << "NTrees  " << t->GetNtrees() << endl;
				//cout << "NBranch " << t->GetNbranches() << endl;
				//ntree = t->GetTree();
				//ntree->Draw("PUnumInteractions>>pileup2", "");
			}
				t->Draw("PUnumInteractions>>pileup", "");
	//LoadTree

		string outname = endstring + samples[i]+".root";
		TFile *outfile = new TFile(outname.c_str(), "RECREATE");
		pileup->Write();
		//pileup2->Write();
		//pileupI->Write();
		outfile->Close();


		// cleanup
		delete pileup;
		delete outfile;
		delete t;

		cout << endl << "finished " << samples[i] << endl << endl;


		/*TFile *files[nfiles];
		TH1D  *histos[nfiles];
		TH1D  *pileup = new TH1D("pileup", "pileup", 51, -0.5, 50.5);
		for(int n=0; n<nfiles; ++n){

			
			std::stringstream numstream;
			numstream << n;
			string currfile = base+samples[i]+"/output_"+numstream.str()+".root";
			cout << "adding file " <<  currfile << endl;
			files[n] = TFile::Open(currfile.c_str());
			histos[n]= (TH1D*) files[n]->Get("pileup");
			pileup->Add(histos[n]);
		}	

		string outname = samples[i]+".root";
		TFile *outfile = new TFile(outname.c_str(), "RECREATE");
		pileup->Write();
		outfile->Close();


		// cleanup
		delete pileup;
		delete outfile;*/
	}
}
