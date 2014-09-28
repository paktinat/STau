#include "TList.h"
#include "MassPlotter.hh"
#include "Corrector.h"
#include "BaseMassPlotter.hh"
#include <vector>
#include "TChain.h"

BaseMassPlotter::BaseMassPlotter(TString outputdir){
  fOutputDir = Util::MakeOutputDir(outputdir);
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);
}

void BaseMassPlotter::init(TString filename){
  if(fVerbose > 0) cout << "------------------------------------" << endl;
  if(fVerbose > 0) cout << "Initializing MassPlotter ... " << endl;
  loadSamples(filename);
}

void BaseMassPlotter::loadSamples(const char* filename){
  fSamples.clear();
  char buffer[200];
  ifstream IN(filename);

  char ParName[100], StringValue[1000];
  float ParValue;

  if(fVerbose > 0) cout << "------------------------------------" << endl;
  if(fVerbose > 0) cout << "Sample File  " << filename << endl;
  int counter(0);
	
  while( IN.getline(buffer, 200, '\n') ){
    // ok = false;
    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'
    }
    if( !strcmp(buffer, "GENERAL") ){
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Path\t%s", StringValue);
      fPath = StringValue;	
      cout <<"My path: " <<fPath << endl;
			
      if(fVerbose >0){
	cout << " ----  " << endl;
	cout << "  Path " << fPath << endl;
      }

    }


    if( !strcmp(buffer, "SAMPLE")){

      sample s;
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Name\t%s", StringValue);
      s.name = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "SName\t%s", StringValue);
      s.sname = TString(StringValue);
			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "ShapeName\t%s", StringValue);
      s.shapename = TString(StringValue);

      IN.getline(buffer, 400, '\n');
      sscanf(buffer, "File\t%s", StringValue);
      TString file =fPath+StringValue;
      if(fVerbose > 3)
	cout<<"my file: "<<file<<endl;

		
      s.tree = new TChain("MassTree"); //(TTree*)f->Get("MassTree");
      ((TChain*)(s.tree))->Add( file , 0 );
      ((TChain*)(s.tree))->LoadTree(0);

      cout << s.tree->GetEntries() << endl;
      s.file = ((TChain*)(s.tree))->GetFile() ;

			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Xsection\t%f", &ParValue);
      s.xsection = ParValue;
			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Kfact\t%f", &ParValue);
      s.kfact = ParValue;
			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Lumi\t%f", &ParValue);
      s.lumi = ParValue;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Type\t%s", StringValue);
      s.type = StringValue;
			
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Color\t%f", &ParValue);
      s.color = ParValue;

      if(s.type!="data"){
	TH1F *h_PUWeights = (TH1F*) s.file->Get("h_PUWeights");
	TH1F *h_Events    = (TH1F*) s.file->Get("h_Events");
	if(h_PUWeights==0 || h_Events==0){
	  cout << "ERROR: sample " << (s.file)->GetName() << " does not have PU and NEvents histos! " << endl;
	  exit(1);
	}
	s.type!="data" ? s.PU_avg_weight = h_PUWeights->GetMean()    : s.PU_avg_weight =1;
	s.type!="data" ? s.nevents       = h_Events   ->GetEntries() : s.nevents       =1;
	delete h_PUWeights;
	delete h_Events;
      } else{
	s.PU_avg_weight =1;
	s.nevents       =1;
      }

      // DON'T DO THIS AT HOME !!!!
      if ( s.name == "T1tttt_mGlu-1000_mLSP-400" ) s.nevents = 20000;
      if ( s.name.Contains("T2bb") ) s.nevents = 10000;
      if ( s.name.Contains("T2tt") ) s.nevents = 50000;
      // DON'T DO THAT AT HOME !!!!

      if ( s.type == "susy" && s.name.Contains("Tau")){
	TH2F * h_SMSEvents = (TH2F*) s.file->Get("h_SMSEvents");
	int binNumber = h_SMSEvents->FindBin(150.0, 150.0);//350,50//saeid
			  
	s.nevents = h_SMSEvents->GetBinContent(binNumber); //saeid
	s.nevents = 1000;//350,50//saeid
	s.PU_avg_weight =1;		//saeid	  
      }//saeid
      else{

	if ( s.type == "susy" && !s.name.Contains("TStau")){//saeid
	  TH2F * h_SMSEvents = (TH2F*) s.file->Get("h_SMSEvents");
	  int binNumber = h_SMSEvents->FindBin(350.0, 50.0);//350,50//saeid
			  
	  s.nevents = h_SMSEvents->GetBinContent(binNumber); //saeid
	  //			  s.nevents = 128333;//350,50//saeid
	  s.PU_avg_weight =1;		//saeid	  
	}//saeid
	else
	  if ( s.type == "susy" && s.name.Contains("TStau")){
	    s.PU_avg_weight =1;	
	    s.nevents = 10000;
	  }
      }
      //s.ReconstructSName();
      if(fVerbose > 0){
	cout << " ---- " << endl;
	cout << "  New sample added: " << s.name << endl;
	cout << "   Sample no.      " << counter << endl;
	cout << "   Short name:     " << s.sname << endl;
	cout << "   File:           " << (s.file)->GetName() << endl;
	cout << "   Events:         " << s.nevents  << endl;
	cout << "   Events in tree: " << s.tree->GetEntries() << endl; 
	cout << "   Xsection:       " << s.xsection << endl;
	cout << "   Lumi:           " << s.lumi << endl;
	cout << "   kfactor:        " << s.kfact << endl;
	cout << "   avg PU weight:  " << s.PU_avg_weight << endl;
	cout << "   type:           " << s.type << endl;
	cout << "   Color:          " << s.color << endl;
      }

      fSamples.push_back(s);
      counter++;
    }
  }
  if(fVerbose > 0) cout << "------------------------------------" << endl;
}




