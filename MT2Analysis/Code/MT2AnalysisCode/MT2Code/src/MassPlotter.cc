/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#include "TList.h"
#include "MassPlotter.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TMath.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TEventList.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TClonesArray.h"

#include "TGraph2D.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TPaveStats.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <time.h> // access to date/time

using namespace std;

 Objectproperties::Objectproperties(TString type_):type(type_), nBins(10),lowVal(-1000.),highVal(0.){
		if(type == "pt"){	
			nBins = 100;
			lowVal = -2.;
			highVal = 500.;
		}
                else if(type == "phi"){
                        nBins = 70;
                        lowVal = -3.5;
                        highVal = 3.5;
                }
                else if(type == "mult"){
                        nBins = 20;
                        lowVal = -0.5;
                        highVal = 19.5;
                }
		else std::cout<<"Bad Type!!!"<<std::endl;


	}
 Objectproperties::Objectproperties():type(""), nBins(10),lowVal(-1000.),highVal(0.){}
//____________________________________________________________________________
MassPlotter::MassPlotter(){
	fSave=true;
	fisPhoton=false;
	fPUReweight=true;
	fStitching=true;
	fbSFReWeight=true;
	myChannel = "TBD";
	
	//Put all channel specific weights true, You can 
	//muTau
	fMuIdSF=true;
	fMuIsoSF=true;
	fMuTrgSF=true;
	fTauTrgSF=true;
	fTauWjetsSF=true;
	fTauEnergySF=true;

        fE0IdIsoSF=true;
        fE1IdIsoSF=true;
        fETrgSF =true;

        //eleMu channel
       fMuTrgSFeleMu=true;
       fMuIdIsoSF=true;
       fEleTrgSF=true;
       fEleIdIsoSF=true;



// Default constructor, no samples are set
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);

}

//____________________________________________________________________________
MassPlotter::MassPlotter(TString outputdir){
// Explicit constructor with output directory
       	fSave=true;
	fisPhoton=false;
	fPUReweight=true;
	fStitching=true;
	fbSFReWeight=true;
	myChannel = "TBD";
	
	//Put all channel specific weights true, You can 
	//muTau
	fMuIdSF=true;
	fMuIsoSF=true;
	fMuTrgSF=true;
	fTauTrgSF=true;
	fTauWjetsSF=true;
	fTauEnergySF=true;

        fE0IdIsoSF=true;
        fE1IdIsoSF=true;
        fETrgSF =true;

       //eleMu channel
       fMuTrgSFeleMu=true;
       fMuIdIsoSF=true;
       fEleTrgSF=true;
       fEleIdIsoSF=true;




// Default constructor, no samples are set
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);
	setOutputDir(outputdir);
}

//____________________________________________________________________________
MassPlotter::MassPlotter(TString outputdir, TString outputfile){
// Explicit constructor with output directory and output file
//	MassPlotter();
	fSave=true;
	fisPhoton=false;
	fPUReweight=true;
	fStitching=true;
	fbSFReWeight=true;
	myChannel = "TBD";
	
	//Put all channel specific weights true, You can 
	//muTau
	fMuIdSF=true;
	fMuIsoSF=true;
	fMuTrgSF=true;
	fTauTrgSF=true;
	fTauWjetsSF=true;
	fTauEnergySF=true;

        fE0IdIsoSF=true;
        fE1IdIsoSF=true;
        fETrgSF =true;

        //eleMu channel
       fMuTrgSFeleMu=true;
       fMuIdIsoSF=true;
       fEleTrgSF=true;
       fEleIdIsoSF=true;



// Default constructor, no samples are set
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);


	setOutputDir(outputdir);
	setOutputFile(outputfile);
}

//____________________________________________________________________________
MassPlotter::~MassPlotter(){
	fOutputFile->Close();
	delete fOutputFile;
}

//____________________________________________________________________________
void MassPlotter::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing MassPlotter ... " << endl;
	//Util::SetStyle();
	loadSamples(filename);

	// map for Z->nunu plots
	// to be used for Z->nunu           to be used for Z->ll
	RemoveLeptMap["misc.MET"]       = "Znunu.METplusLeptsPtReco";
}

//___________________________________________________________________________
void MassPlotter::makeSmallCopy(unsigned int nevents, unsigned int mysample, TString cuts, TString trigger){
  // 	if(fSamples.size()<mysample) return;


  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;
  /*
   TFile *pileup_data = new TFile(GETDATALOCALPATH(Certification/pileUp_data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.root),"READ");
  
   TH1F* pileup_data_histo = (TH1F*) pileup_data->Get("pileup");

   TFile *pileup_mc = new TFile(GETDATALOCALPATH(Certification/pileUp_mc/MC2012PU_S10_60bins.root),"READ");
  
   TH1F* pileup_mc_histo = (TH1F*) pileup_mc->Get("pileup");
  
   pileup_data_histo->Scale(1.0/pileup_data_histo->Integral());

   pileup_mc_histo->Scale(1.0/pileup_mc_histo->Integral());
  
// //   pileup_data_histo->Divide(pileup_mc_histo);

*/

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    if(fSamples.size() > mysample && ii != mysample)
      continue;
    //  if(fSamples[ii].type != "data" || fSamples[ii].name != "Run2012D1F")
    // continue;


    fMT2tree = new MT2tree();
    fSamples[ii].tree->SetBranchAddress("MT2tree", &fMT2tree);


    cout << "making small copy of tree " << fSamples[ii].name << endl;
    //TFile *newfile = new TFile(fOutputDir+"/"+fSamples[ii].name+".root","recreate");
    TString fileName = fSamples[ii].file->GetName();

    //    fileName =  fileName.ReplaceAll(".root", "_NBJetsCSVM0_MET30.root");
    fileName =  fileName.ReplaceAll(".root", "_DoubleTauNew.root");

    TFile *newfile = new TFile(fOutputDir+"/"+fileName,"recreate");
    TTree *newtree = fSamples[ii].tree->CloneTree(0);
    TH1F *h_PUWeights = (TH1F*) fSamples[ii].file->Get("h_PUWeights");
    TH1F *h_Events    = (TH1F*) fSamples[ii].file->Get("h_Events");
    
    //    h_PUWeights->Reset();

    if(fSamples[ii].file->Get("h_SMSEvents")){
      h_SMSEvents    = (TH2F*) fSamples[ii].file->Get("h_SMSEvents");
      h_mSugraEvents = (TH2F*) fSamples[ii].file->Get("h_mSugraEvents");
    }

    TString myCuts = cuts;
    if( fSamples[ii].type=="data") myCuts += " && " + trigger;
 

    fSamples[ii].tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();

    for (unsigned int i=0;i<nentries; i++) {
      if(i > nevents) continue;

      if ( i % 100000 == 0 ){ 	
	fprintf(stdout, "\rProcessed events: %6d of %6d ", i + 1, nentries);
	fflush(stdout);
      }
 
      fSamples[ii].tree->GetEntry(myEvtList->GetEntry(i));
      /*   
      float oldWeight = fMT2tree->pileUp.Weight;
      `
//    if(fSamples[ii].sname == "Wtolnu" || (fSamples[ii].shapename == "ZJetsToLL" && fSamples[ii].name != "DYToLL_M10To50")){
// 	int VBosonID = 24;
// 	if(fSamples[ii].sname == "Wtolnu")
// 	  VBosonID = 24;
// 	else
// 	  VBosonID = 23;

// 	fMT2tree->pileUp.Weight = newWeight * fMT2tree->NewWeightFromStitching(VBosonID);
//       }else
//  	fMT2tree->pileUp.Weight = newWeight;

       // if(fMT2tree->eleMu[0].ele0Ind >= 0)
	 //	 fMT2tree->eleMu[0].EleIdIsoSF  = fMT2tree->ele[fMT2tree->eleMu[0].ele0Ind].GetEleIDISOSFelemu();
     
       fMT2tree->pileUp.Weight = oldWeight * fMT2tree->weightTauTau();
*/
       newtree->Fill();

       //    h_PUWeights->Fill(newWeight);
   
//       cout<<"old weight "<<oldWeight<<" new weight "<<newWeight<<endl;
    }
    cout << endl;
	
    h_PUWeights->Write();
    h_Events->Write();
    if(fSamples[ii].file->Get("h_SMSEvents")){
      h_SMSEvents->Write();
      h_mSugraEvents->Write();
    }

    newfile->Write();
	
    delete newfile;
  }
}

//___________________________________________________________________________
void MassPlotter::makePlots(){
	double dPhisplit[4]={0., 1.0, 1.5, 3.142};
	MakeMT2PredictionAndPlots(false, dPhisplit);  //, 2.5);  // not cleaned, fudgefactor 2
}


//__________________________________________________________________________
void MassPlotter::MakeMT2PredictionAndPlots(bool cleaned , double dPhisplit[]){//, double fudgefactor){
	TString cleanflag;
	if(cleaned) cleanflag="cleaned";
	else        cleanflag="notcleaned";
	std::ostringstream o0, o1, o2, o3;
	o0 << dPhisplit[0];
	o1 << dPhisplit[1];
	o2 << dPhisplit[2];
	o3 << dPhisplit[3];
	TString dPhirange1 = (TString) o0.str() + "to" +(TString) o1.str();
	TString dPhirange2 = (TString) o2.str() + "to" +(TString) o3.str();
	TString dPhirange  = (TString) o0.str()+(TString) o1.str() + (TString) o2.str() +(TString) o3.str();
	
	std::vector<sample>  QCDSamples;
	std::vector<sample>  Samples_NOTData;
	std::vector<sample>  QCDandDataSamples;
	for(unsigned int i=0; i< fSamples.size(); ++i){
		if(fSamples[i].sname=="QCD"){
			QCDSamples.push_back(fSamples[i]);
			QCDandDataSamples.push_back(fSamples[i]);
		}
		if(fSamples[i].type!="data"){
			Samples_NOTData.push_back(fSamples[i]);
		}
		if(fSamples[i].type=="data"){
			QCDandDataSamples.push_back(fSamples[i]);
		}
	}

	//	double ZerotoPi[4]    ={0., 3.142, 3.142, 3.142};
	//	double controlsplit[4]={dPhisplit[2], dPhisplit[3], dPhisplit[3], dPhisplit[3]};

	std::ostringstream cutStream;
	cutStream  
		  << "misc.MET > 30"                                << "&&"
		  << "misc.HT > 300"                                << "&&"
		  << "misc.Jet0Pass == 1"                           << "&&"
		  << "misc.Jet1Pass == 1"                           << "&&"
		  << "misc.PassJetID == 1"                          << "&&"
		  << "misc.Vectorsumpt<70"                          << "&&"
		  << "misc.MinMetJetDPhi>0.3"                       << "&&" 
		  << "misc.EcalDeadCellBEFlag==1"                   << "&&"
		  << "misc.HBHENoiseFlag == 1"                   ;

	TString cuts = cutStream.str().c_str();
	
	//                 variable                            cuts njets  nlepts HLT title     bins               flip_order  log  composite    ratio  stacked overlay 
	//	MakePlot(fSamples,"misc.MT2" ,                         cuts, -3,   0,    "",  "MT2" , gNMT2bins, gMT2bins , false,  true ,  true,      true,  true,  false);
//	MakePlot(fSamples,"hemi[0].lv1.M()" ,                  cuts, -3,   0,    "",  "h1.M", 30,       0, 1000  , false,  true ,  true,      true,  true,  false);
//	MakePlot(fSamples,"hemi[0].lv2.M()" ,                  cuts, -3,   0,    "",  "h1.M", 30,       0, 1000  , false,  true ,  true,      true,  true,  false);

}

//________________________________________________________________________

void MassPlotter::makeplot(TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,  TString xtitle,
			   const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type){

  MakePlot(fSamples, var, maincuts, basecuts, njets, nbjets, nleps, HLT, xtitle, nbins, min, max, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro, type);
  
}

void MassPlotter::makePlot(TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT,  TString xtitle,
			   const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type){

  MakePlot(fSamples, var, cuts, njets, nbjets, nleps, HLT, xtitle, nbins, min, max, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro, type);
  
}

//________________________________________________________________________
void MassPlotter::Makeplot(TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT, TString xtitle, 
			   const int nbins,  const double *bins,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type){

  MakePlot(fSamples, var, maincuts, basecuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro, type);

}

void MassPlotter::MakePlot(TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT, TString xtitle, 
			   const int nbins,  const double *bins,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type){

  MakePlot(fSamples, var, cuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro, type);

}

// ________________________________________________________________________
void MassPlotter::PrintWEfficiency(int sample_index ,TString process,  std::string lept, unsigned int nevents, bool includeTaus){
	sample Sample = fSamples[sample_index];

	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	cout << "printing W efficiency: \n"
	     << Sample.name << endl;	
	std::cout << setfill('-') << std::setw(70) << "" << std::endl;

	enum counters_t { count_begin, all=count_begin, presel, NJetsIDLoose , PassJetID, MinMetJetDPhi, VectorSumPt, MT2Cut,  count_end };
	Monitor counters[count_end];
	TString lablx[count_end] = {"all events", "presel", "NJetsIDLoose" , "PassJetID", "MinMetJetDPhi", "VectorSumPt", "MT2Cut"};

	string to_measure;
	if     (lept == "ele" && process=="W" )       {to_measure = "W->enu  recoed";} 
	else if(lept == "muo" && process=="W" )       {to_measure = "W->munu recoed";} 
	else if(lept == "ele" && process=="Top")      {to_measure = "Top W->enu  recoed";}
	else if(lept == "muo" && process=="Top")      {to_measure = "Top W->munu  recoed";}

	fMT2tree = new MT2tree();
	Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
	unsigned int nentries =  Sample.tree->GetEntries();
	unsigned int nbytes = 0, nb = 0;
	for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
		nb = Sample.tree->GetEntry(jentry);   nbytes += nb;
		Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
		if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;

		bool leptfound(false);
		bool eventgood(false);
		bool acceptance(false);
		if(process=="W"){
			if(lept=="ele" && fMT2tree->GenLeptFromW(11, 0 , 1000,includeTaus)==true)                                              eventgood =true;
			if(lept=="muo" && fMT2tree->GenLeptFromW(13, 0 , 1000,includeTaus)==true)                                              eventgood =true;
			if(lept=="ele" && fMT2tree->GenLeptFromW(11,10,  2.4 ,includeTaus)==true)                                              acceptance=true;
			if(lept=="muo" && fMT2tree->GenLeptFromW(13,10,  2.4 ,includeTaus)==true)                                              acceptance=true;
			if(lept=="ele" && fMT2tree->NEles==1 && fMT2tree->NMuons==0 && fMT2tree->GenLeptFromW(11, 0, 1000,includeTaus)==true && fMT2tree->ele[0].lv.Pt()>10)  leptfound =true;
			if(lept=="muo" && fMT2tree->NMuons==1&& fMT2tree->NEles== 0 && fMT2tree->GenLeptFromW(13, 0, 1000,includeTaus)==true && fMT2tree->muo[0].lv.Pt()>10)  leptfound =true;
		}else if(process=="Top"){ 
			if(lept=="ele" && fMT2tree->TopDecayModeResult(11)==true)                                                 eventgood =true;
			if(lept=="muo" && fMT2tree->TopDecayModeResult(13)==true)                                                 eventgood =true;
			if(lept=="ele" && fMT2tree->TopDecayModeResult(11)==true && fMT2tree->SLTopAccept(10, 2.4)==true )        acceptance =true;
			if(lept=="muo" && fMT2tree->TopDecayModeResult(13)==true && fMT2tree->SLTopAccept(10, 2.4)==true )        acceptance =true;
			if(lept=="ele" && fMT2tree->NEles==1 && fMT2tree->NMuons==0 && eventgood &&fMT2tree->ele[0].lv.Pt()>10)             leptfound =true;
			if(lept=="muo" && fMT2tree->NMuons==1&& fMT2tree->NEles== 0 && eventgood &&fMT2tree->muo[0].lv.Pt()>10)             leptfound =true;
		}

		Double_t weight = fMT2tree->pileUp.Weight;

		// all events
		if(eventgood)   counters[all].fill("all events", weight);
		if(acceptance)  counters[all].fill("acceptance", weight);
		if(leptfound)   counters[all].fill(to_measure, weight);
		
		// presel
		if(fMT2tree->misc.MET                     < 30    )  continue;
		if(fMT2tree->misc.HT                      < 300   )  continue;
		if(fMT2tree->misc.Jet0Pass                ==0     )  continue;
		if(fMT2tree->misc.Jet1Pass                ==0     )  continue;
//		if(fMT2tree->misc.LeadingJPt              < 150   )  continue;
		if(fMT2tree->misc.SecondJPt               < 100   )  continue;
//		if(fMT2tree->NBJets                       < 1     )  continue;
//		if(fMT2tree->misc.HBHENoiseFlag           ==0     )  continue;
		if(fMT2tree->misc.CrazyHCAL               ==1     )  continue;
		if(eventgood)   counters[presel].fill("presel", weight);
		if(acceptance)  counters[presel].fill("acceptance", weight);
		if(leptfound)   counters[presel].fill(to_measure, weight);
		
		
		// NJetsIDLoose
		if(fMT2tree->NJetsIDLoose                 <3     ) continue;  
		if(eventgood)   counters[NJetsIDLoose].fill("NJetsIDLoose", weight);
		if(acceptance)  counters[NJetsIDLoose].fill("acceptance", weight);
		if(leptfound)   counters[NJetsIDLoose].fill(to_measure, weight);

		// PassJetID
		if(fMT2tree->misc.PassJetID               ==0     ) continue;
		if(eventgood)   counters[PassJetID].fill("PassJetID", weight);
		if(acceptance)  counters[PassJetID].fill("acceptance", weight);
		if(leptfound)   counters[PassJetID].fill(to_measure, weight);
		
		// MinMetJetDPhi
		if(fMT2tree->misc.MinMetJetDPhi <0.3 )       continue;
		if(eventgood)   counters[MinMetJetDPhi].fill("MinMetJetDPhi", weight);
		if(acceptance)  counters[MinMetJetDPhi].fill("acceptance", weight);
		if(leptfound)   counters[MinMetJetDPhi].fill(to_measure, weight);
		
		// VectorSumPt
		if(fMT2tree->misc.Vectorsumpt             >70     )  continue;
		if(eventgood)   counters[VectorSumPt].fill("VectorSumPt", weight);
		if(acceptance)  counters[VectorSumPt].fill("acceptance", weight);
		if(leptfound)   counters[VectorSumPt].fill(to_measure, weight);
		
		// MT2
//		if(fMT2tree->misc.MT2             <400     )  continue;
//		if(fMT2tree->misc.MT2             <150     )  continue;
		if(fMT2tree->misc.MT2 <200 || fMT2tree->misc.MT2>400 )  continue;
//		if(fMT2tree->misc.MT2 <150 || fMT2tree->misc.MT2>300 )  continue;
//		if(fMT2tree->misc.MT2 <100 || fMT2tree->misc.MT2>150  )  continue;
		if(eventgood)   counters[MT2Cut].fill("MT2Cut", weight);
		if(acceptance)  counters[MT2Cut].fill("acceptance", weight);
		if(leptfound)   counters[MT2Cut].fill(to_measure, weight);
		
	}

	// print stats     
	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	std::cout << "Statistics" << std::endl;
	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
	for ( counters_t iCount=count_begin; iCount<count_end; iCount = counters_t(iCount+1) ) {
		counters[iCount].print();
		std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;  
	}
	std::cout << setfill('=') << std::setw(70) << "" << std::endl;

	// fill histo
	TH1D* h_all = new TH1D("h_all", "", count_end, 0., (double) count_end );
	TH1D* h_acc = new TH1D("h_acc", "", count_end, 0., (double) count_end );
	TH1D* h_rec = new TH1D("h_rec", "", count_end, 0., (double) count_end );
	TH1D* h_1   = new TH1D("h_1"  , "", count_end, 0., (double) count_end );
	TH1D* h_2   = new TH1D("h_2"  , "", count_end, 0., (double) count_end );
	TH1D* h_3   = new TH1D("h_3"  , "", count_end, 0., (double) count_end );
	for(int i=0; i<count_end; ++i){
		h_1    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_2    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_3    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_all->SetBinContent(i+1,counters[i].counts((string) lablx[i]));
		h_acc->SetBinContent(i+1,counters[i].counts("acceptance"));
		h_rec->SetBinContent(i+1,counters[i].counts(to_measure));
	}
	h_1    ->Sumw2();
	h_2    ->Sumw2();
	h_3    ->Sumw2();
	h_1    ->Divide(h_acc,h_all);	
	h_2    ->Divide(h_rec,h_acc);	
	h_3    ->Divide(h_rec,h_all);	
	TCanvas *col  = new TCanvas("Wevents", "", 0, 0, 900, 700);
	col -> cd();
	h_1    ->SetLineColor(kBlue);
	h_1    ->SetMarkerStyle(20);
	h_1    ->SetMinimum(0);
	h_1    ->SetMaximum(1);
	h_1    ->SetDrawOption("E");
	h_1    ->Draw();
	string name  ="Wevents_"+lept+"_acceptance_"+(string) Sample.name; 
	h_2    ->SetLineColor(kBlue);
	h_2    ->SetMinimum(0);
	h_2    ->SetMaximum(1);
	h_2    ->SetMarkerStyle(20);
	h_2    ->SetDrawOption("E");
	h_2    ->Draw();
	name  ="Wevents_"+lept+"_reco_"+(string) Sample.name; 
	h_3    ->SetLineColor(kBlue);
	h_3    ->SetMinimum(0);
	h_3    ->SetMaximum(1);
	h_3    ->SetMarkerStyle(20);
	h_3    ->SetDrawOption("E");
	h_3    ->Draw();
	name  ="Wevents_"+lept+"_prob_"+(string) Sample.name; 

	cout << "lept " << lept << " process " << process << endl;	
	//__________________________________
	if        (lept=="ele" && process=="W") {
		fWpred.Wenu_acc        =h_1->GetBinContent(count_end);
		fWpred.Wenu_acc_err    =h_1->GetBinError(count_end); 
		fWpred.Wenu_rec        =h_2->GetBinContent(count_end);
		fWpred.Wenu_rec_err    =h_2->GetBinError(count_end);  
		fWpred.Wenu_prob       =h_3->GetBinContent(count_end);
		fWpred.Wenu_prob_err   =h_3->GetBinError(count_end);  
	}else if  (lept=="ele" && process=="Top"){
		fWpred.TopWenu_acc     =h_1->GetBinContent(count_end);
		fWpred.TopWenu_acc_err =h_1->GetBinError(count_end); 
		fWpred.TopWenu_rec     =h_2->GetBinContent(count_end);
		fWpred.TopWenu_rec_err =h_2->GetBinError(count_end);  
		fWpred.TopWenu_prob    =h_3->GetBinContent(count_end);
		fWpred.TopWenu_prob_err=h_3->GetBinError(count_end);  
	}else if  (lept=="muo" && process=="W") {
		fWpred.Wmunu_acc       =h_1->GetBinContent(count_end);
		fWpred.Wmunu_acc_err   =h_1->GetBinError(count_end); 
		fWpred.Wmunu_rec       =h_2->GetBinContent(count_end);
		fWpred.Wmunu_rec_err   =h_2->GetBinError(count_end); 
		fWpred.Wmunu_prob      =h_3->GetBinContent(count_end);
		fWpred.Wmunu_prob_err  =h_3->GetBinError(count_end);  
	}else if  (lept=="muo" && process=="Top") {
		fWpred.TopWmunu_acc    =h_1->GetBinContent(count_end);
		fWpred.TopWmunu_acc_err=h_1->GetBinError(count_end); 
		fWpred.TopWmunu_rec    =h_2->GetBinContent(count_end);
		fWpred.TopWmunu_rec_err=h_2->GetBinError(count_end); 
		fWpred.TopWmunu_prob    =h_3->GetBinContent(count_end);
		fWpred.TopWmunu_prob_err=h_3->GetBinError(count_end);  
	}

	//_________________________________
	delete h_1;
	delete h_2;
	delete h_3;
	delete h_all;
	delete h_acc;
	delete h_rec;
	delete col;
}

//________________________________________________________________________


//________________________________________________________________________

void MassPlotter::plotSig(TString var, TString theCuts, TString xtitle, int nbins, double min, double max, bool flip_order, int type , int LowerCut){

        cout<<" flip_order "<<flip_order<<" LowerCut "<<LowerCut<<endl;

        TString varname = Util::removeFunnyChar(var.Data());

	TH1D*    h_mc_sum    = new TH1D   (varname+"mc_sum", "", nbins, min, max );
	TH1D*    h_susy      = new TH1D   (varname+"susy"  , "", nbins, min, max );	
	h_mc_sum  ->Sumw2();
	h_susy    ->Sumw2();
	// vector of all histos
	vector<TH1D*> h_samples;

	for(size_t i = 0; i < fSamples.size(); ++i){
	        if(!(fSamples[i].type=="susy") && !(fSamples[i].type=="mc")) continue;

		h_samples.push_back(new TH1D(varname+"_"+fSamples[i].name, "", nbins, min, max));
		h_samples[i] -> Sumw2();
                
             
               
                   Float_t weight =0;
  
                if  (fSamples[i].type=="mc") weight= fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
			 
		if(fSamples[i].type=="susy")
	        weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);

	        if(fSamples[i].type=="data")
		weight = 1.0;

	        TString btagweight = "SFWeight.BTagCSV40eq0";//"1.0";

	        TString selection;	
		if(     fSamples[i].type=="mc" && fPUReweight && fbSFReWeight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), theCuts.Data());
		else if(fSamples[i].type=="mc" && fPUReweight                ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    theCuts.Data());
		else if(fSamples[i].type=="mc" &&                fbSFReWeight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), theCuts.Data());
		else                                                            selection = TString::Format("(%.15f) * (%s)",                 weight,                    theCuts.Data()); 

              
                TString variable  = TString::Format("%s>>%s",var.Data(),h_samples[i]->GetName());
	  
		if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;


		int nev = fSamples[i].tree->Draw( variable.Data(), selection.Data(),"goff");
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				   << "\t->Integral() : "  <<  h_samples[i]->Integral(0,nbins+1) << endl;
		
		if( fSamples[i].type=="mc"        )   h_mc_sum->Add(h_samples[i]);
                  
		else if( fSamples[i].type=="susy" )   h_susy  ->Add(h_samples[i]);

	}
        if (fVerbose>2)cout<<"----------------------------------bg----"<<h_mc_sum->Integral(0,nbins+1) <<endl;
   TCanvas *MyC = new TCanvas("MyC", "MyC");
  MyC->Divide(2,1);
  MyC->cd(1);
  TGraphAsymmErrors* sig1 = plotSig(h_susy, h_mc_sum, xtitle, "Lower Cut", type, /*sys = */ 0.1);
  sig1->Draw("ACP");

  MyC->cd(2);
  TGraphAsymmErrors* sig2 = plotSig(h_susy, h_mc_sum, xtitle, "Upper Cut", type, /*sys = */ 0.1);
  sig2->Draw("ACP");

}
//__________________________________________________________________________
void MassPlotter::CompSamples(TString var, TString cuts, TString optcut, bool RemoveLepts, 
		              TString xtitle, const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	CompSamples(fSamples, var, cuts, optcut, RemoveLepts, xtitle, nbins, bins, add_underflow, logflag, scale_factor, normalize);
}
//___________________________________________________________________________________________________________________
void MassPlotter::compSamples(TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	CompSamples(fSamples, var, cuts, optcut, RemoveLepts, xtitle, nbins, min, max, add_underflow, logflag, scale_factor, normalize);

}
// _______________________________________________________________________
void MassPlotter::CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	double bins[nbins];
  	bins[0] = min;
  	for(int i=1; i<=nbins; i++)
    	bins[i] = min+i*(max-min)/nbins;
  	CompSamples(Samples, var, cuts, optcut, RemoveLepts, xtitle, nbins, bins, add_underflow, logflag, scale_factor, normalize);

}
// __________________________________________________________________________
void MassPlotter::CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize){

        TString varname = Util::removeFunnyChar(var.Data());
	TH1D*    histos[2];
	cout << Samples.size() << endl;
	histos[0] = new TH1D(varname+Samples[0].sname  , "", nbins, bins);
	histos[1] = new TH1D(varname+Samples[1].sname  , "", nbins, bins);

	histos[0] -> Sumw2(); 
	histos[1] -> Sumw2();
	histos[0] -> SetLineColor(kBlue); 
	histos[1] -> SetLineColor(kRed); 
	histos[0] -> SetLineWidth(4);
	histos[1] -> SetLineWidth(4);

	// legend
	TLegend* Legend1 = new TLegend(.71,.54,.91,.92);

	for(int i=0; i<2; ++i){
		Double_t weight = scale_factor * Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);
		if(fVerbose>2) cout << "GetHisto: looping over " << Samples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
		// exchange variable to plot with the corresponding one with removed leptons;	
		TString variable; 
		TString selection; 
		if(RemoveLepts && (Samples[i].sname == "DYToLL" || Samples[i].sname=="QCD_2")){
			MapType::iterator iter = RemoveLeptMap.begin();
			iter = RemoveLeptMap.find(var);
			if (iter != RemoveLeptMap.end() ){
				variable = iter -> second;
			} else {cout << "found RemoveLepts==true, but no corresponding map" << endl; variable = var;}
		}else {
			variable = var;
		}
		if(optcut!="_" && Samples[i].sname == "DYToLL"){
			TString newcuts = cuts + "&&" + optcut;  
			if(Samples[i].type!="data" && fPUReweight==true ) selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,newcuts.Data());
			else                                              selection      = TString::Format("(%.15f) * (%s)"              ,weight,newcuts.Data()); 
		} else {
			if(Samples[i].type!="data" && fPUReweight==true ) selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,cuts.Data());
			else                                              selection      = TString::Format("(%.15f) * (%s)"              ,weight,cuts.Data()); 
		}

		variable          = TString::Format("%s>>%s",variable.Data(),histos[i]->GetName());
		
			  
		if(fVerbose>2) cout << "+++++ Drawing " << variable  << endl
				    << "\twith cuts: "  << selection << endl;

		int nev = Samples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case
		if(add_underflow) {
			histos[i]->SetBinContent(1,
				    histos[i]->GetBinContent(0) + histos[i]->GetBinContent(1));
			histos[i]->SetBinError(1,
				  sqrt(histos[i]->GetBinError(0)*histos[i]->GetBinError(0)+
				       histos[i]->GetBinError(1)*histos[i]->GetBinError(1) ));
		}
		histos[i]->SetBinContent(histos[i]->GetNbinsX(),
					    histos[i]->GetBinContent(histos[i]->GetNbinsX()  )+ 
					    histos[i]->GetBinContent(histos[i]->GetNbinsX()+1) );
		histos[i]->SetBinError(histos[i]->GetNbinsX(),
					  sqrt(histos[i]->GetBinError(histos[i]->GetNbinsX()  )*
					       histos[i]->GetBinError(histos[i]->GetNbinsX()  )+
					       histos[i]->GetBinError(histos[i]->GetNbinsX()+1)*
					       histos[i]->GetBinError(histos[i]->GetNbinsX()+1)  ));
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				    << "\t->Integral() : "  <<  histos[i]->Integral() << endl;
	     	Legend1 ->AddEntry(histos[i], Samples[i].sname, "l"); 
	}
	
	plotRatio(histos[0], histos[1], logflag, normalize, varname+"_"+Samples[0].sname+"_"+Samples[1].sname, Legend1 , xtitle, "");
	
	
}

//________________________________________________________________________

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  MakePlot(Samples, var, maincuts, basecuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro, type);

}

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  MakePlot(Samples, var, cuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, 	saveMacro, type);

}

//________________________________________________________________________

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double *bins, 
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type){

  TString basecuts = "";
  TString maincuts = cuts;
  MakePlot(Samples, var, maincuts, basecuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro, type);

}

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double *bins, 
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro, int type){

        TH2D *h_PN_MLSP_MChi = new TH2D("h_PN_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
	h_PN_MLSP_MChi->Sumw2();

        TH2D *h_N_MLSP_MChi = new TH2D("h_N_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
	h_N_MLSP_MChi->Sumw2();

	int SusySampleInd = -1;

        TString varname = Util::removeFunnyChar(var.Data());

	TString nMinCut = var;
	nMinCut += TString::Format(">=%g", bins[0]);

	TString nJets, nJetsVar = "NJetsIDLoose40";
	if (njets>=10) {
	  nJets =  "(" + nJetsVar + TString::Format(">=%d",njets/10);
	  nJets += "&&"+ nJetsVar + TString::Format("<=%d",njets%10)+")";
	}
	else{
	  nJets = nJetsVar + (njets < 0 ? ">=" : "==");
	  nJets = nJets + TString::Format("%d",abs(njets));
	}

	if(njets == -10)
	  nJets = nJetsVar + ">= 0";


	//saeid
// 	TString nBJets = "NBJets40CSVM";    // nbjets = -10  --> >=0 b-tags
	TString nBJets = "NBJetsCSVT";    // nbjets = -10  --> >=0 b-tags
	//saeid
	nBJets += nbjets < 0 ? ">=" : "==";
	nBJets += nbjets==-10 ? "0" : TString::Format("%d",abs(nbjets));

	TString  nLeps;
	if     (nleps < 0 )  nLeps = " && (NEles + NMuons) >=";
	else if(nleps >=0  ) nLeps = " && (NEles + NMuons) ==";
	nLeps += TString::Format("%d",abs(nleps));
	if     (nleps ==-10) nLeps = " "; 
	if     (nleps ==-11) nLeps = " && NEles ==1 && NMuons ==0"; 
	if     (nleps ==-13) nLeps = " && NEles ==0 && NMuons ==1"; 

	TString varname2 = "h_";
	
	THStack* h_stack     = new THStack(varname2, "");
  	TH1D*    h_data      = new TH1D   (varname2+"data"  , "", nbins, bins );
	TH1D*    h_mc_sum    = new TH1D   (varname2+"mc_sum", "", nbins, bins );
	TH1D*    h_susy      = new TH1D   (varname2+"susy"  , "", nbins, bins );	

	// h_data
	h_data -> Sumw2();
	h_data -> SetMarkerStyle(20);
	h_data -> SetMarkerColor(kBlack);
	h_data -> SetLineColor(kBlack);
	h_data -> SetStats(false);

	h_mc_sum -> SetFillStyle(3004);
	h_mc_sum -> SetFillColor(kBlack);
	h_mc_sum -> SetStats(0);
	
	h_mc_sum    ->Sumw2();
	h_susy      ->Sumw2();


	// vector of all histos
	vector<TH1D*> h_samples;

// 	// arrays of composited stuff
// 	TH1D    *h_composited[7];
// 	//TString  cnames[5] = {"QCD", "W/Z/#gamma production", "Top production", "susy", "data"};
// 	//int      ccolor[5] = {401, 418, 602, 0, 632};
// 	TString  cnames[7] = {"QCD", "W+jets", "Z+jets", "Top",   "WW+jets", "susy", "data"};
// 	TString  cnames2[7] = {"QCD", "Wjets", "Zjets",  "Top",   "WWjets",  "susy", "data"};
// 	int      ccolor[7] = { 401,   417,       419,      855,    603,     0,     632};


	int SignalInd = NumberOfSamples-2;//index of signal in the canmes[]

	// arrays of composited stuff
	TH1D    *h_composited[NumberOfSamples];
	//TString  cnames[5] = {"QCD", "W/Z/#gamma production", "Top production", "susy", "data"};
	//int      ccolor[5] = {401, 418, 602, 0, 632};
	TString  cnames[NumberOfSamples]  = {"QCD",     "W",    "ZX",  "Top",       "WW",  "Higgs", "susy", "data"};
	TString  cnames2[NumberOfSamples] = {"QCD", "Wjets", "Zjets",  "Top",   "WWjets",  "Higgs", "susy", "data"};
	int      ccolor[NumberOfSamples]  = {  401,     417,     419,    855,        603,   kRed,        0,    632};



	vector<TH1D*> h_signals;
	int n_sig=0;
	TString  snames[5] = {"signal1","signal2","signal3","signal4","signal5"};
	int      scolor[5] = {     kRed+3,   kOrange-3,   kViolet-6,   kBlack,   kRed-4};
	for (int i=0; i<NumberOfSamples; i++){
	  h_composited[i] = new TH1D(varname2+"_"+cnames2[i], "", nbins, bins);
	  h_composited[i] -> Sumw2();
	  h_composited[i] -> SetFillColor  (stacked ? ccolor[i] : 0);
	  h_composited[i] -> SetLineColor  (ccolor[i]);
	  if (!stacked) {
	    h_composited[i] -> SetLineWidth(4);
	  }
	  h_composited[i] -> SetMarkerColor(ccolor[i]);
	  h_composited[i] -> SetStats(false);
	}
	if(fisPhoton){
		h_composited[4] -> SetFillColor(kViolet-3); 
		h_composited[4] -> SetLineColor(kViolet-3);
		cnames[4]="Photons";
	}

	// legend
	TLegend* Legend1 = new TLegend(.71,.68,.91,.92);
	Legend1->SetName(varname2+"_legend");
	//TLegend* Legend1 = new TLegend(.3,.5,.6,.88);
	bool SetfPUReweight = fPUReweight; // To be able to control the pile weight from the run_MassPlot.C
	for(size_t i = 0; i < Samples.size(); ++i){

		h_samples.push_back(new TH1D(varname2+"_"+Samples[i].name, "", nbins, bins));
		h_samples[i] -> Sumw2();
		h_samples[i] -> SetFillColor(stacked ? Samples[i].color : 0);
		h_samples[i] -> SetLineColor(Samples[i].color);
		if (!stacked) {
		  h_samples[i] -> SetLineWidth(4);
		}
		h_samples[i] -> SetMarkerColor(Samples[i].color);
		h_samples[i] -> SetStats(false);
		if(Samples[i].type == "susy" ){
			h_samples[i] -> SetLineColor(kBlack);
			h_samples[i] -> SetLineStyle(kDotted);
			h_composited[SignalInd] -> SetLineColor(kBlack);
			h_composited[SignalInd] -> SetLineStyle(kDotted);
		}
	
		Double_t weight=0;

		

		//------added be esmaeel
		if (Samples[i].name == "QCD-Pt-15-20-MuEnriched"){
		  fPUReweight = false;}
		else fPUReweight = SetfPUReweight;
		//------added be esmaeel
		if(fPUReweight) weight = Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents*Samples[i].PU_avg_weight);
		else            weight = Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);

		if(fVerbose>2) cout << "===========================================================" <<  endl;
		if(fVerbose>2) cout << "MakePlot: looping over " << Samples[i].sname << "-----------------------------------" <<  endl;
		if(fVerbose>2) cout << "  +++++++ xsection:    "    << Samples[i].xsection << " k-fact " << Samples[i].kfact << endl;
		if(fVerbose>2) cout << "  +++++++ tot events:  "  << Samples[i].nevents  << " avg pu weight " << Samples[i].PU_avg_weight << endl;
		if(fVerbose>2) cout << "  +++++++ PU reweight: "  << fPUReweight << endl;
		if(fVerbose>2) cout << "  +++++++ weight:      " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s>>%s",var.Data(),h_samples[i]->GetName());
		TString theCuts = nJets + "&&" + nBJets + nLeps;
		if(!add_underflow) theCuts = theCuts + "&&" + nMinCut;
		theCuts = theCuts + "&&" + maincuts;
		if(basecuts!="") theCuts = theCuts + "&&" + basecuts;
		
		TString selection;

		if(Samples[i].type=="data"){
		  if(HLT!="")
		    theCuts += " &&("+HLT+")"; // triggers for data
		  selection = TString::Format("(%.15f) * (%s)",                 weight,   theCuts.Data());
		}else{
		  
		  TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
		  if(nbjets>=0 && nbjets<=3) btagweight = TString::Format("SFWeight.BTagCSV40eq%d",abs(nbjets));
		  else if(nbjets>=-3)        btagweight = TString::Format("SFWeight.BTagCSV40ge%d",abs(nbjets));

		  TString ChannelSpecificSF = "1.00";
		  
		  if(myChannel == "muTau"){
		    if(fMuIdSF)
		      ChannelSpecificSF += "*muTau[0].muIdSF";
		    if(fMuIsoSF)
		      ChannelSpecificSF += "*muTau[0].muIsoSF";
		    if(fMuTrgSF)
		      ChannelSpecificSF += "*muTau[0].muTrgSF";
		    if(fTauTrgSF)
		      ChannelSpecificSF += "*muTau[0].tauTrgSF";
		    if(fTauWjetsSF && (Samples[i].sname == "Wtolnu"))
		      ChannelSpecificSF += "*muTau[0].tauWjetsSF";	  
		    if(fTauEnergySF )
		      ChannelSpecificSF += "*muTau[0].tauEnergySF";	  
		  }
		  else
		    if(myChannel == "ee"){
		      if(fE0IdIsoSF)
			ChannelSpecificSF += "*doubleEle[0].Ele0IdIsoSF";
		      if(fE1IdIsoSF)
			ChannelSpecificSF += "*doubleEle[0].Ele1IdIsoSF";
		      if(fETrgSF)
			ChannelSpecificSF += "*doubleEle[0].DiEleTrgSF";
		    }
                 else if (myChannel == "eleMu"){
                   
                  if (fMuTrgSFeleMu)
                  ChannelSpecificSF += "*eleMu[0].MuTrgSF";
                  if (fMuIdIsoSF)
                  ChannelSpecificSF += "*eleMu[0].MuIdIsoSF";
                  if (fEleTrgSF)
                  ChannelSpecificSF += "*eleMu[0].EleTrgSF";
                   if ( fEleIdIsoSF )
                  ChannelSpecificSF += "*eleMu[0].EleIdIsoSF";

		 }
            

		  if(fStitching && (Samples[i].sname == "Wtolnu" || (Samples[i].shapename == "ZJetsToLL" && Samples[i].name != "DYToLL_M10To50"))){
		    weight = Samples[i].lumi;
		    if(fPUReweight) weight /= Samples[i].PU_avg_weight;
		  }
		  
		  if(fPUReweight && fbSFReWeight) selection = TString::Format("(%.15f*pileUp.Weight*%s*%s) * (%s)",weight, btagweight.Data(), ChannelSpecificSF.Data(), theCuts.Data());
		  else if(fPUReweight)  selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",   weight,                    ChannelSpecificSF.Data(), theCuts.Data());
		  else if(fbSFReWeight) selection = TString::Format("(%.15f*%s*%s) * (%s)",              weight, btagweight.Data(), ChannelSpecificSF.Data(), theCuts.Data());
		  else                  selection = TString::Format("(%.15f*%s) * (%s)",                 weight,                    ChannelSpecificSF.Data(), theCuts.Data()); 
		  
		  if(     Samples[i].type=="susy"){
		    if(fPUReweight) weight = Samples[i].lumi / (Samples[i].PU_avg_weight);
		    else            weight = Samples[i].lumi;
		    TString SUSYselection;
		    if(fPUReweight && fbSFReWeight) SUSYselection = TString::Format("(%.15f*pileUp.Weight*%s*%s) * (%s)",weight, btagweight.Data(), ChannelSpecificSF.Data(), theCuts.Data());
		    else if(fPUReweight)  SUSYselection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",   weight,                    ChannelSpecificSF.Data(), theCuts.Data());
		    else if(fbSFReWeight) SUSYselection = TString::Format("(%.15f*%s*%s) * (%s)",              weight, btagweight.Data(), ChannelSpecificSF.Data(), theCuts.Data());
		    else                  SUSYselection = TString::Format("(%.15f*%s) * (%s)",                 weight,                    ChannelSpecificSF.Data(), theCuts.Data()); 

		    TString SUSYvariable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_PN_MLSP_MChi->GetName());

		    if(fVerbose>2) cout << "  +++++++ Drawing SUSY " << SUSYvariable  << endl
					<< "  +++++++ with cuts:   " << setw(40)  << SUSYselection << endl;
		
		    fSamples[i].tree->Draw(SUSYvariable, SUSYselection, "goff");

		    SusySampleInd = i; 	  
		    
		    SUSYselection = TString::Format(" (%s)",                 theCuts.Data()); 

		    SUSYvariable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_N_MLSP_MChi->GetName());

		    if(fVerbose>2) cout << "  +++++++ Drawing SUSY " << SUSYvariable  << endl
					<< "  +++++++ with cuts:   " << setw(40)  << SUSYselection << endl;
		
		    fSamples[i].tree->Draw(SUSYvariable, SUSYselection, "goff");
		  
		  }
		 
		}
		if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;

		int nev = Samples[i].tree->Draw(variable.Data(),selection.Data(),"goff");


		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case


		AddOverAndUnderFlow(h_samples[i], true, add_underflow);

		if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;

		/// event count with errors
		TH1F * clone = (TH1F*)h_samples[i]->Clone();
		clone->Rebin(clone->GetNbinsX());
		if(fVerbose>2) cout << "  +++++++ Events:       " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
		
		if (Samples[i].sname.Contains("QCD")) {
		  h_composited[0]->Add(h_samples[i]);
		}
		else if(Samples[i].sname.Contains("Top") ){
		  h_composited[3]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="Wtolnu") {
		  h_composited[1]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="DY") {
		  h_composited[2]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="VV" || Samples[i].sname=="VGamma" || Samples[i].sname=="Photons") {
		  h_composited[4]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="Higgs") {
		  h_composited[5]->Add(h_samples[i]);  
		}
		else if(Samples[i].type.Contains("susy")){
		  h_composited[SignalInd]->Add(h_samples[i]);
		  cnames[SignalInd] = Samples[i].sname;
		  snames[n_sig] = Samples[i].sname;

		  char buf[10];
		  sprintf(buf,"no%d", n_sig);
		  //itoa(n_sig, buf, 10);
		  TString numstring = TString(buf);
		  h_signals.push_back(new TH1D(varname2+"_"+snames[n_sig]+"_"+numstring, "", nbins, bins));
		  h_signals[n_sig] -> Sumw2();
		  h_signals[n_sig] -> SetFillColor  (stacked ? scolor[n_sig] : 0);
		  h_signals[n_sig] -> SetLineColor  (scolor[n_sig]);
		  h_signals[n_sig] -> SetLineWidth  (3);
		  if (!stacked) {
		    h_signals[n_sig] -> SetLineWidth(4);
		  }
		  h_signals[n_sig] -> SetMarkerColor(scolor[n_sig]);
		  h_signals[n_sig] -> SetStats(false);

		  h_signals[n_sig]->Add(h_samples[i]);
		  n_sig++;
		}
		else if(Samples[i].type.Contains("data")){
		  h_composited[NumberOfSamples - 1]->Add(h_samples[i]);
		}

	}

	if (composited){
	  if(flip_order){
	    for (int i=(NumberOfSamples-2); i>=0; --i){
	      if (!stacked)  {
	        h_composited[i]->Scale(1./h_composited[i]->Integral());
	        h_composited[i] -> SetMinimum(0.00005);
	        //if (fVerbose>2) cout << " h_composited[" << i<< "]->Integral = " << h_composited[i]->Integral()  << endl;
	      }
	      if( i!=SignalInd || !overlaySUSY) h_stack  -> Add(h_composited[i]);
	      if( i==SignalInd)                 h_susy   -> Add(h_composited[i]);
	      else        	        h_mc_sum -> Add(h_composited[i]);
	      if( i==SignalInd && overlaySUSY)
	        //Legend1 ->AddEntry(h_composited[i], cnames[i]  , "f");  
		for (int iii=0; iii<n_sig; iii++) Legend1 ->AddEntry(h_signals[iii], snames[iii]  , "f");  
	      else
	        Legend1 ->AddEntry(h_composited[i], cnames[i], stacked ? "f" : "l");
	    }
	  }
	  else{
	    for (int i=0; i<(NumberOfSamples - 1); ++i){
	      if (!stacked)  {
	        h_composited[i]->Scale(1./h_composited[i]->Integral());
	        h_composited[i] -> SetMinimum(0.00005);
	        //if (fVerbose>2) cout << " h_composited[" << i<< "]->Integral = " << h_composited[i]->Integral()  << endl;
	      }
	      if( i!=SignalInd || !overlaySUSY) h_stack  -> Add(h_composited[i]);
	      if( i==SignalInd)                 h_susy   -> Add(h_composited[i]);
	      else        	      h_mc_sum -> Add(h_composited[i]);
	      if( i==SignalInd && overlaySUSY){
	        //if(h_composited[i]->GetEntries()>0) Legend1 ->AddEntry(h_composited[i], cnames[i] + (overlayScale ? overlayScale==1? "":

		//  					      TString::Format(" x %.0f",overlayScale) :
		//  					      " scaled to data")   , "f");  
		for (int iii=0; iii<n_sig; iii++) Legend1 ->AddEntry(h_signals[iii], snames[iii]+  (overlayScale ? overlayScale==1? "":
		  					      TString::Format(" x %.0f",overlayScale) :
		  					      " scaled to data")   , "f");  
	      }else{
	        if(h_composited[i]->GetEntries()>0) Legend1 ->AddEntry(h_composited[i], cnames[i], stacked ? "f" : "l");
	      }
	    }
	  }
	  h_data->Add(h_composited[NumberOfSamples - 1]);
	  if(h_data->Integral()>0 && stacked){
	    Legend1     ->AddEntry(h_data, "data", "p");
	  }
	  if(fVerbose > 2 && composited) {
		           cout << "------------------------------------"                << endl
	                        << "QCD Integral:      " << h_composited[0]->Integral()  << endl
			        << "W+jets Integral:   " << h_composited[1]->Integral()  << endl
			        << "Z+jets Integral:   " << h_composited[2]->Integral()  << endl
			        << "Top Integral:      " << h_composited[3]->Integral()  << endl
			        << "WW +jets Integral: " << h_composited[4]->Integral()  << endl
			        << "Higgs    Integral: " << h_composited[5]->Integral()  << endl
				<< "TOTAL BG:          " << h_composited[0]->Integral()+h_composited[1]->Integral()+h_composited[2]->Integral()+h_composited[3]->Integral()+h_composited[4]->Integral()+h_composited[5]->Integral() <<endl
			        << "SUSY:              " << h_composited[NumberOfSamples - 2]->Integral()  << endl
			        << "Data:              " << h_composited[NumberOfSamples - 1]->Integral()  << endl;
	  
			   //Same, but with errors
			   string OverallSamples[NumberOfSamples + 1] = {"QCD","W+jets","Z+Jets","Top","WW+jets","Higgs","SUSY","Data","TOTAL BG"};
			   for(int os=0; os<=NumberOfSamples; os++){
			     string tSample = OverallSamples[os];
			     if(tSample!="TOTAL BG"){
			       TH1F* clone = (TH1F*) h_composited[os]->Clone();
			       clone->Rebin(clone->GetNbinsX());
			       if(fVerbose>2) cout << tSample << ": " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
			     }
			     else{
			       TH1F* clone = (TH1F*) h_composited[0]->Clone();
			       clone->Add(h_composited[1]);
			       clone->Add(h_composited[2]);
                               clone->Add(h_composited[3]);
			       clone->Add(h_composited[4]);
			       clone->Add(h_composited[5]);
			       clone->Rebin(clone->GetNbinsX());
                               if(fVerbose>2) cout << tSample << ": " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
			     }
			   }
			  // save nevents for W background prediction
			  if(nleps==-11){
			  	fWpred.QCD_bg_e   = h_composited[0]->Integral();	
				fWpred.W_bg_e     = h_composited[1]->Integral();
				fWpred.Z_bg_e     = h_composited[2]->Integral();
				fWpred.Top_bg_e   = h_composited[3]->Integral();
				fWpred.Other_bg_e = h_composited[4]->Integral();
				fWpred.Other_bg_e = h_composited[5]->Integral();
			  }else if(nleps==-13){
			  	fWpred.QCD_bg_mu   = h_composited[0]->Integral();	
				fWpred.W_bg_mu     = h_composited[1]->Integral();
				fWpred.Z_bg_mu     = h_composited[2]->Integral();
				fWpred.Top_bg_mu   = h_composited[3]->Integral();
				fWpred.Other_bg_mu = h_composited[4]->Integral();
				fWpred.Other_bg_mu = h_composited[5]->Integral();
			  } 
	  }
	}
	else{
	  for(unsigned int i=0; i<Samples.size(); ++i){
	    if (!stacked) {
	      h_samples[i]->Scale(1./h_samples[i]->Integral());
	      h_samples[i] -> SetMinimum(0.0001);
	    }
	    if((Samples[i].sname=="QCD" )){
	      h_samples[i] -> SetMinimum(0.01);
	      h_stack->Add(h_samples[i]);
	      h_mc_sum->Add(h_samples[i]);
	      Legend1 ->AddEntry(h_samples[i], Samples[i].sname, "f");
	    }
	    if(Samples[i].sname!="QCD" && Samples[i].type!="data"){
	      if(Samples[i].type=="susy") h_susy->Add(h_samples[i]);
	      if(Samples[i].type!="susy" || !overlaySUSY) h_stack->Add(h_samples[i]);
	      if(Samples[i].type!="susy") h_mc_sum->Add(h_samples[i]);
	      if(Samples[i].type=="susy" && overlaySUSY){
	         if(h_samples[i]->Integral()>0) Legend1 ->AddEntry(h_samples[i], Samples[i].sname + (overlayScale ? 
								       TString::Format(" x %.0f",overlayScale) :
								       " scaled to data")   , "f") ;
	      }else{
		 if(h_samples[i]->Integral()>0) Legend1 ->AddEntry(h_samples[i], Samples[i].sname, stacked ? "f" : "l");
	      }	
	    }
	    if(Samples[i].type == "data"){
	      h_data -> Add(h_samples[i]);
	    }
	    
	  }
	  if(h_data->Integral()>0){
	    Legend1     ->AddEntry(h_data, "data", "l");
	  }
	}

	TH1F *h_PN_Bkg  = (TH1F*)h_mc_sum->Clone();
	h_PN_Bkg->SetName("h_PN_Bkg");
	h_PN_Bkg->Rebin(h_PN_Bkg->GetNbinsX());
 
	TH1F *h_PN_Data  = (TH1F*)h_data->Clone();
	h_PN_Data->SetName("h_PN_Data");
	h_PN_Data->Rebin(h_PN_Data->GetNbinsX());

 	h_SMSEvents->Rebin2D(4, 4);
 	h_PN_MLSP_MChi->Divide(h_SMSEvents);


	TH2* hXsec = (TH2*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");//CharginoChargino
	//TH2* hXsec = (TH2*) TFile::Open("referenceXSecs.root")->Get("StSt_8TeV_NLONLL_LSP"); hXsec->Rebin2D(2, 2);hXsec->Scale(0.25);//StauStau
	h_PN_MLSP_MChi->Multiply(hXsec);

	TString fileName = fOutputDir;
	if(!fileName.EndsWith("/")) fileName += "/";
	Util::MakeOutputDir(fileName);
	fileName = fileName + "countingForExclusion_" + myChannel +"_" +xtitle + "_Histos.root";
	TFile *savefile = new TFile(fileName.Data(), "RECREATE");
	savefile ->cd();

	h_stack->SetName("h_stack");

	h_stack->Write();
	h_PN_Bkg->Write();
	h_PN_Data->Write();
	h_PN_MLSP_MChi->Write();
	h_N_MLSP_MChi->Write();

	savefile->Close();
  
	std::cout << "Saved histograms in " << savefile->GetName() << std::endl;

	TString ytitle = "Events";

	TString outname = myChannel;
	if(saveHistos){
		TString fileName = fOutputDir;
		if(!fileName.EndsWith("/")) fileName += "/";
		Util::MakeOutputDir(fileName);
		fileName = fileName + outname + ".root";
		TFile *savefile = new TFile(fileName.Data(), "RECREATE");
		savefile ->cd();
		for(unsigned int nh = 0; nh<7; ++nh)                h_composited[nh]->Write();
		for(unsigned int nh = 0; nh<h_signals.size(); ++nh) h_signals[nh]->Write();
		h_stack->Write();
		h_mc_sum->Write();
		h_data->Write();
		h_susy->Write();
		Legend1->Write();
		savefile->Close();
		std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
	}

	if(!stacked) {
	  outname = outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + "_shape";
	  printHisto(h_stack, h_data, h_mc_sum, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps, false);

	}
	else if (!overlaySUSY){
	  outname = outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "");
	  printHisto(h_stack, h_data, h_mc_sum, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps);
	  if (ratio) 
	    plotRatioStack(h_stack,  h_mc_sum, h_data,  logflag, false, outname, Legend1, xtitle, ytitle, njets, nleps);
	}
	else {
	  outname =  outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "") + "_overlay";	  
	  //printHisto(h_stack, h_data, h_mc_sum, h_susy, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nbjets, nleps, overlayScale);
	  printHisto(h_stack, h_data, h_mc_sum, h_signals, n_sig, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nbjets, nleps, overlayScale);
	  if (ratio) 
	    plotRatioStack(h_stack,  h_mc_sum, h_data, h_susy,  logflag, false, outname, Legend1, xtitle, ytitle, njets, nbjets, nleps, overlayScale, saveMacro);

	}

  TCanvas *MyC = new TCanvas("MyC", "MyC");
  MyC->Divide(2,1);
  MyC->cd(1);
  TGraphAsymmErrors* sig1 = plotSig(h_susy, h_mc_sum, xtitle, "Lower Cut", type, /*sys = */ 0.1);
  sig1->Draw("ACP");

  MyC->cd(2);
  TGraphAsymmErrors* sig2 = plotSig(h_susy, h_mc_sum, xtitle, "Upper Cut", type, /*sys = */ 0.1);
  sig2->Draw("ACP");

  TString fileName1 = fOutputDir;
  if(!fileName1.EndsWith("/")) fileName1 += "/";
  Util::MakeOutputDir(fileName1);
  fileName1 = fileName1 + "plotsig_" + xtitle +".root";
  TFile *savefile1 = new TFile(fileName1.Data(), "RECREATE");
  savefile1->cd();
  sig1->Write();
  sig2->Write();
  savefile1->Close();



}
 

//_________________________________________________________________________________
void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets, int nleps, TString saveMacro){
	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	
	TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);
	c1 -> cd();
	
	float border = 0.3;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(name+"_plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(name+"_ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);

	hstack->SetMinimum(0.02);
	hstack ->Draw("hist");
	//h2     ->Draw("sameEX0");
        h2     ->Draw("sameE1"); //LEO MOD

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.1);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d Jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	

	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	h_ratio->SetLineWidth(2);
	//h_ratio->SetFillColor(kBlue);//leo
	//h_ratio->SetLineColor(kBlue);
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
	//h_ratio ->SetMarkerSize(0.1);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	h_ratio ->DrawCopy("E"/*2"*/);
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.15);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "data / MC");
	
	// x axis title
	lat.SetTextSize(0.2);
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	//	TString save=name+"_ratio";
	if(fSave) Util::Print(c1, xtitle, fOutputDir, fOutputFile);

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + xtitle +".root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  c1->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;

  TString fileName1 = fOutputDir;
  if(!fileName1.EndsWith("/")) fileName1 += "/";
  Util::MakeOutputDir(fileName1);
  fileName1 = fileName1 + xtitle +".png";
  TFile *savefile1 = new TFile(fileName1.Data(), "RECREATE");
  savefile1 ->cd();
  c1->Write();
  savefile1->Close();
  std::cout << "Saved histograms in " << savefile1->GetName() << std::endl;

}
//_________________________________________________________________________________
void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale, TString saveMacro){
  //LEO TRUE USE THIS

	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);

	
	//TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);
	TCanvas* c1 = new TCanvas(name+"c_ratio","",0,0,600,600 /*37, 60,636,670*/);
	c1->SetFrameLineWidth(1);
	c1 -> cd();
	
// 	float border = 0.2;
	// 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(name+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1 /*0.00, border, 1.00, 1.00, 0, 0*/);
 	//p_plot->SetBottomMargin(0.05);
	//p_plot->SetTopMargin(0.09);
	//p_plot->SetLeftMargin(0.1669107);
      	//p_plot->SetRightMargin(0.02);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(name+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441/*     0.00, 0.05, 1.00, border, 0, 0*/);
 	//p_ratio->SetTopMargin(0.03);
 	//p_ratio->SetBottomMargin(0.05/*5*/);
	//p_ratio->SetRightMargin(0.02);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);

	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2.5*max;
	else max = 1.5*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h1_orig->GetBinWidth(1) -h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1))<0.01){
		double binwidth = h1_orig->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h1_orig->GetBinWidth(1) << " " << h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	hstack->GetYaxis()->SetTitle(yTitle.str().c_str());
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.3);

	
	//MT2_bSel[0]->SetTitleSize(0.03);
	///MT2_bSel[0]->SetTitleOffset(1.);
	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	h2    ->Draw("sameE");
	h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
	h3->SetFillColor(0);
	h3->SetLineStyle(kDotted);
	h3->SetLineWidth(4);
	h3->Draw("samehist");

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);

	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";


// 	TString text ="";
// 	text = fMT2Analysis?  "M_{T2} Analysis                                          ":"";
// 	text +=fMT2bAnalysis? "M_{T2}b Analysis                                         ":"";
// 	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
// 	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
	TitleBox.DrawLatex(0.13,0.943,text.Data());
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0305);
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	LumiBox.DrawLatex(0.68,0.943,"#sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//standard
	//LumiBox.DrawLatex(0.49,0.943,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS Preliminary
	//LumiBox.DrawLatex(0.62,0.943,"CMS, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)h1_orig->Clone("h1_copy");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->GetYaxis()->SetRangeUser(0,2);
	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
	h_ratio_mc->GetXaxis()->SetTitle(xtitle);
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc	->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
	h_ratio_mc->GetYaxis()->SetNdivisions(509);

	
	h_ratio_mc->SetFillStyle(3001);
	h_ratio_mc->Draw("E2");
	h_ratio ->DrawCopy("Esame");//LEO MOD
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

// 	TString save=name+"_ratio";
// 	if(fSave)Util::Print(c1, save, fOutputDir, fOutputFile);
// 	if(saveMacro != ""){
// 	  TString newXtitle = Util::removeFunnyChar(xtitle.Data());
// 	  c1->SaveAs(newXtitle + "_" + save + "." + saveMacro);
// 	}


	if(fSave)Util::Print(c1, xtitle, fOutputDir, fOutputFile);
	//ssss


	//        c1->SaveAs(xtitle+".root");
	//        c1->SaveAs(xtitle+".png");


  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + xtitle +".root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  c1->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;

  TString fileName1 = fOutputDir;
  if(!fileName1.EndsWith("/")) fileName1 += "/";
  Util::MakeOutputDir(fileName1);
  fileName1 = fileName1 + xtitle +".png";
  TFile *savefile1 = new TFile(fileName1.Data(), "RECREATE");
  savefile1 ->cd();
  c1->Write();
  savefile1->Close();
  std::cout << "Saved histograms in " << savefile1->GetName()<< std::endl;

	//	if(saveMacro != "")
	//	  	  c1->SaveAs(xtitle + "." + saveMacro);
	//	  	  c1->SaveAs(xtitle + ".png");
}

//_________________________________________________________________________________
void MassPlotter::plotRatio(TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle){
	// define canvas and pads
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	
	TCanvas* c1 = new TCanvas("c1","", 20,100,1000,700);
	c1 -> cd();
	
	float border = 0.3;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	gPad->SetFillStyle(0);
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);
	h2     ->Draw();
	h1     ->Draw("same");
	
	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	

	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	h_ratio->SetLineWidth(2);
	h_ratio->SetFillColor(kBlue);
	h_ratio->SetLineColor(kBlue);
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
	h_ratio ->SetMarkerSize(0.1);

 	h_ratio ->Divide(h1, h2);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	gPad->SetFillStyle(0);
	h_ratio ->DrawCopy("E2");
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.07);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "ratio");
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	//gPad->SetLogy(1);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	if(fSave)Util::PrintEPS(c1, save, fOutputDir);	

}
//___________________________________________________________________________
void MassPlotter::printHisto(TH1* h_orig, TString canvname,  Option_t *drawopt, bool logflag){
	TH1D* h = (TH1D*)h_orig->Clone("h1_copy");
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	if(logflag) gPad->SetLogy(1);
	h->SetMinimum(0.01);
	gPad->SetRightMargin(0.15);
	gPad->SetLogz();
	h->SetStats(false);
	h->DrawCopy(drawopt);
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	h->Write();
	delete col;
	delete h;

}

//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TString canvname, Option_t *drawopt, bool logflag){
	
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) gPad->SetLogy(1);
	h->Draw(drawopt);
	gPad->RedrawAxis();
	col ->Update();
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets, int nleps, bool stacked){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(stacked ? 0.05 : 0.0001);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	if(stacked){
	  h_data  ->SetMaximum(max);
	  h_mc_sum->SetMaximum(max);
	  h       ->SetMaximum(max);
	}

	h->Draw(stacked ? drawopt : "histnostack");
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0 && stacked) {
		h_data       ->Draw("same");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d Jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;


	  //	TString newXtitle = Util::removeFunnyChar(xtitle.Data());
	  //	col->SaveAs(newXtitle + "_" + "NoRatio.C");
	
}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, TH1* h_susy, TLegend* legend,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets, int nbjets, int nleps, float overlayScale){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	TLegend *leg = (TLegend*) legend->Clone("leg");

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(0.05);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		h_susy   -> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0) {
		h_data       ->Draw("sameE");
	}
	h_susy->Scale(overlayScale ? overlayScale : h_data->Integral() / h_susy->Integral());
	h_susy->SetLineStyle(kDotted);
	h_susy->SetFillColor(0);
	h_susy->Draw("samehist");
	if(leg != NULL ){
		leg -> SetY1NDC(0.68);
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.0305);
	//TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
	//text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

// 	TLatex TitleBox;
// 	TitleBox.SetNDC();
// 	TitleBox.SetTextSize(0.0305);
// 	TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
// 	text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
// 	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
// 	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
// 	TitleBox.DrawLatex(0.18,0.943,text.Data());


	h->GetXaxis()->SetTitle(xtitle);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.1);
	
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h_data->GetBinWidth(1) -h_data->GetBinWidth(h_data->GetNbinsX()-1))<0.01){
		double binwidth = h_data->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h_data->GetBinWidth(1) << " " << h_data->GetBinWidth(h_data->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	h->GetYaxis()->SetTitle(yTitle.str().c_str());
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;


	//
	//	TString newXtitle = Util::removeFunnyChar(xtitle.Data());
	//	col->SaveAs(newXtitle + "_" + "NoRatio.C");

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, vector<TH1D*> h_sig, int nsig, TLegend* legend,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	TLegend *leg = (TLegend*) legend->Clone("leg");

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(0.05);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		for (int i=0; i<nsig; i++)  h_sig[i]-> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0) {
		h_data       ->Draw("sameE");
	}
	for (int i=0; i<nsig; i++) {
	  h_sig[i]->Scale(overlayScale ? overlayScale : h_data->Integral() / h_sig[i]->Integral());
	  h_sig[i]->SetLineStyle(2);
	  //h_sig[i]->SetLineWidth(2);
	  h_sig[i]->SetFillColor(0);
	  h_sig[i]->Draw("samehist");
	}
	if(leg != NULL ){
		leg -> SetY1NDC(0.68);
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.0305);
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());//standard
//	TitleBox.DrawLatex(0.247,0.895,text.Data());//for CMS Preliminary
//	TitleBox.DrawLatex(0.257,0.895,text.Data());//for CMS, and 1 lepton
//	TitleBox.DrawLatex(0.18,0.943,text.Data());//for CMS
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0305);
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	LumiBox.DrawLatex(0.615,0.943,"#sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//standard
//	LumiBox.DrawLatex(0.18,0.943,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";//for CMS Preliminary
//	LumiBox.DrawLatex(0.18,0.943,"CMS, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS, and 1 lepton
//	LumiBox.DrawLatex(0.55,0.943,"CMS, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS


// 	TLatex TitleBox;
// 	TitleBox.SetNDC();
// 	TitleBox.SetTextSize(0.0305);
// 	TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
// 	text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
// 	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
// 	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
// 	TitleBox.DrawLatex(0.18,0.943,text.Data());

	h->GetXaxis()->SetTitle(xtitle);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.1);
	
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h_data->GetBinWidth(1) -h_data->GetBinWidth(h_data->GetNbinsX()-1))<0.01){
		double binwidth = h_data->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h_data->GetBinWidth(1) << " " << h_data->GetBinWidth(h_data->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	h->GetYaxis()->SetTitle(yTitle.str().c_str());
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

	//	TString newXtitle = Util::removeFunnyChar(xtitle.Data());
	//	col->SaveAs(newXtitle + "_" + "NoRatio.C");

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag , TString xtitle, TString ytitle){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h->SetMinimum(0.01);
		h_data->SetMinimum(0.01);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h     ->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	if(h_data->Integral()>0) {
		h_data  ->Draw("same, E1");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	col ->Update();
	if(fSave) Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	//	TString newXtitle = Util::removeFunnyChar(xtitle.Data());
	//	col->SaveAs(newXtitle + "_" + "NoRatio.C");
	//	delete col;




}
//____________________________________________________________________________
void MassPlotter::loadSamples(const char* filename){
	fSamples.clear();
	char buffer[200];
	ifstream IN(filename);

	char StringValue[1000];
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
			//			if(fVerbose > 3)
			  cout<<"my file: "<<file<<endl;
		
			s.tree = new TChain("MassTree"); //(TTree*)f->Get("MassTree");
			((TChain*)(s.tree))->Add( file , 0 );
			((TChain*)(s.tree))->LoadTree(0);

			if(fVerbose > 3)
			  cout << s.tree->GetEntries() << endl;

			s.file = ((TChain*)(s.tree))->GetFile() ;
			if(s.file == NULL)
			  s.file = new TFile(file);

// 			s.file->Print("ALL");
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
			  TH1F *h_PUWeights = NULL;
			  TH1F *h_Events    = NULL;

			  TObjArray *fileElements=((TChain*)(s.tree))->GetListOfFiles();
			  TIter next(fileElements);
			  TChainElement *chEl=0;
			  while (( chEl=(TChainElement*)next() )) {
			   TFile* f = TFile::Open(chEl->GetTitle());
			  //   cout<<chEl->GetTitle()<<endl;
			    //    f.Print("ALL");
			    if( h_PUWeights == NULL ){
			      gROOT->cd();
			      h_PUWeights = (TH1F*)(f->Get("h_PUWeights")->Clone("ClonedhPUWeights"));
			    }
			    else
			      h_PUWeights->Add( (TH1F*)(f->Get("h_PUWeights") ) );

		// 	    h_PUWeights->Print("");

			    if( h_Events == NULL )
			      h_Events = (TH1F*) (f->Get("h_Events")->Clone("ClonedhEvents") ) ;
			    else
			      h_Events->Add( (TH1F*) (f->Get("h_Events")) );
			    delete f;
			  }
			  if(fileElements->GetEntries() == 0){
			    h_PUWeights =(TH1F*) s.file->Get("h_PUWeights") ;
			    h_Events = (TH1F*) s.file->Get("h_Events") ;
			  }
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

			if ( s.type == "susy" && s.name.Contains("Stau")){
			  h_SMSEvents = NULL ; 
			  
			  TObjArray *fileElements=((TChain*)(s.tree))->GetListOfFiles();
			  TIter next(fileElements);
			  TChainElement *chEl=0;
			  while (( chEl=(TChainElement*)next() )) {
			    TFile* f = TFile::Open(chEl->GetTitle());
			    
			    if( h_SMSEvents == NULL ){
			      gROOT->cd();
			      h_SMSEvents = (TH2F*) (f->Get("h_SMSEvents")->Clone("ClonedhSMSEvents") );}
			    else
			      h_SMSEvents->Add( (TH2F*)(f->Get("h_SMSEvents") ) );
			  }

			  if(fileElements->GetEntries() == 0){
			    h_SMSEvents=(TH2F*)(s.file->Get("h_SMSEvents"));
			  }

			  int binNumber = h_SMSEvents->FindBin(150.0, 150.0);//350,50//saeid
			  
			  s.nevents = h_SMSEvents->GetBinContent(binNumber); //saeid
			  s.nevents = 10000;//350,50//saeid
			  s.PU_avg_weight =1;		//saeid	  
			}//saeid
			else{

			if ( s.type == "susy" && !s.name.Contains("TStau")){//saeid
			  h_SMSEvents = (TH2F*) s.file->Get("h_SMSEvents");
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


void MassPlotter::TauContamination(int sample_index, unsigned int nevents, int flag, TString cuts, TString trigger, int signalContamination){
  float sys = 0.05; //systematic uncertanty on the fake rate and efficiency
  float IsoCut = 1.5;//cut on the isolation 0.5:VLoose, 1.5:Loose...
//   float IsoCutLoose = IsoCut - 1.0;
  float maxDPhiTauMET = 11.5;//To cut on the DeltaPhi between the tau and MET

  setFlags(flag); 
  TH1::SetDefaultSumw2();

  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    MT2[i] = new TH1F(varname+"_"+cnames[i], "", 150, 0, 750);
    MT2[i] -> SetFillColor  (ccolor[i]);
    MT2[i] -> SetLineColor  (ccolor[i]);
    MT2[i] -> SetLineWidth  (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[5] -> SetMarkerStyle(20);
  MT2[5] -> SetMarkerColor(kBlack);
  MT2[5] -> SetLineColor(kBlack);
  
  MT2[4] -> SetFillStyle(3004);
  MT2[4] -> SetFillColor(kBlack);

  TH2F* hEtaPtJet           = new TH2F("hEtaPtJet",           "hEtaPtJet",          1,0,1000.0,1,0,2.5);//tauable jets in MC in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFake       = new TH2F("hEtaPtTauFake",       "hEtaPtTauFake",      1,0,1000.0,1,0,2.5);//taus in MC in QCD dominated region  (MT2 < 70)
  TH2F* hEtaPtJetData       = new TH2F("hEtaPtJetData",       "hEtaPtJetData",      1,0,1000.0,1,0,2.5);//tauable jets in data in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFakeData   = new TH2F("hEtaPtTauFakeData",   "hEtaPtTauFakeData",  1,0,1000.0,1,0,2.5);//taus in data in QCD dominated region  (MT2 < 70)

  TH1F* hMT2TauFound        = new TH1F("hMT2TauFound",   "hMT2TauFound",   NumberOfMT2Bins - 1, xMT2bin);//Rec/Id taus
  TH1F* hMT2JetSig          = new TH1F("hMT2JetSig",     "hMT2JetSig",     NumberOfMT2Bins - 1, xMT2bin);//tauable jets in signal region (MC)
  TH1F* hMT2TauSig          = new TH1F("hMT2TauSig",     "hMT2TauSig",     NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (MC)
  TH1F* hMT2JetSigData      = new TH1F("hMT2JetSigData", "hMT2JetSigData", NumberOfMT2Bins - 1, xMT2bin);//tauable jets in signal region (data)
  TH1F* hMT2TauSigData      = new TH1F("hMT2TauSigData", "hMT2TauSigData", NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (data)
  TH1F* hRatioTight         = new TH1F("hRatioTight",       "hRatioTight",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut
  TH1F* hRatioRelaxed       = new TH1F("hRatioRelaxed",       "hRatioRelaxed",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut

  TH1F* hNAllTaus           = new TH1F("hNAllTaus",      "hNAllTaus",      NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus
  TH1F* hNAllTauJets        = new TH1F("hNAllTauJets",   "hNAllTauJets",   NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus matched with a rec jet

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(int ii = 0; ii < (int) fSamples.size(); ii++){
    float nAllTausBeforeThisSample = hNAllTaus->Integral();
    if(sample_index != -1 && ii != sample_index)
      continue;
    
    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data"){
      data = 1;
      myCuts += " && " + trigger;
    }

    if(Sample.type == "susy" && signalContamination == 0)// || Sample.sname == "QCD") 
      continue;
     
    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   

    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry); 
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
 
      float weight = Weight;
      if(data == 1)
 	weight = 1.0;
      else
	if(Sample.type != "susy")
	  weight *= (fMT2tree->pileUp.Weight * fMT2tree->SFWeight.BTagCSV40ge1 * fMT2tree->SFWeight.TauTagge1/Sample.PU_avg_weight);//

      int Filled = 0;

      float METPhi = fMT2tree->misc.METPhi; //To be used if a cone close to met is used to look for taus
      if(Sample.type != "susy"){
      for(int t = 0; t < fMT2tree->NTaus; t++){

 	if(fMT2tree->tau[t].MT > 100 || fMT2tree->tau[t].Isolation3Hits < IsoCut || Filled == 1)
 	  continue;

	Filled = 1;

 	float DeltaPhi_ = Util::DeltaPhi(fMT2tree->tau[t].lv.Phi(), fMT2tree->misc.METPhi);

	if(DeltaPhi_ > maxDPhiTauMET)
	  continue;

	if(data == 1){
	    
	  MT2[5]->Fill(fMT2tree->misc.MT2, weight);//data
	    
	}else{
	  MT2[4]->Fill(fMT2tree->misc.MT2, weight);
	    
	  if(Sample.sname == "Top")  
	    MT2[3]->Fill(fMT2tree->misc.MT2, weight);
	  else
	    if(Sample.sname == "DY")	
	      MT2[2]->Fill(fMT2tree->misc.MT2, weight);
	    else
	      if(Sample.sname == "Wtolnu")
		MT2[1]->Fill(fMT2tree->misc.MT2, weight);
	      else
		if(Sample.sname == "QCD")
		  MT2[0]->Fill(fMT2tree->misc.MT2, weight);
	  //	  if(fMT2tree->misc.MT2 >= MT2Min)
	  //cout<<"Tau:: "<<fMT2tree->misc.MT2<<endl;
	}}}

      Filled = 0;
      int GenTau = 0;
      if(Sample.type != "susy"){
       if(fMT2tree->misc.TopDecayMode != 48){//skip di_tau_had  events

      for(int j = 0; j < fMT2tree->NGenLepts; j++){
	if(abs(fMT2tree->genlept[j].ID) != 15 || abs(fMT2tree->genlept[j].MID) != 24)
	  continue;
		

	int HadTau = 1;
	    
	for(int i = 0; i < fMT2tree->NGenLepts; i++){
	  if((abs(fMT2tree->genlept[i].ID) == 11 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID && abs(fMT2tree->genlept[i].GMID) == 24) || 
	     (abs(fMT2tree->genlept[i].ID) == 13 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID && abs(fMT2tree->genlept[i].GMID) == 24)){
	    HadTau = 0;
	  }
	}

 	if(HadTau == 0 || Filled != 0){
// 	  if(Sample.sname == "Top") cout<<"ele "<<fMT2tree->TopDecayModeResult(11)<<endl;
// 	  if(Sample.sname == "Top") cout<<"mu "<<fMT2tree->TopDecayModeResult(13)<<endl;
 	  continue;
	}

	GenTau++;

	Filled = 1;
	if((mT2 && fMT2tree->misc.MinMetJetDPhi   >  0.3) || (mT2b && fMT2tree->misc.MinMetJetDPhi4   >  0.3))
	  hRatioTight->Fill(fMT2tree->misc.MT2, weight);
	else
	  hRatioRelaxed->Fill(fMT2tree->misc.MT2, weight);	  
	  
	
	hNAllTaus->Fill(fMT2tree->misc.MT2, weight);
	
	//To find the jet matched to this gen had tau
	int jetIndex = -1;
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();

	float genleptPhi = fMT2tree->genlept[j].lv.Phi();

	for(int k = 0; k < fMT2tree->NJets; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->jet[k].lv.Eta(), genleptPhi, fMT2tree->jet[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}	

	if(minDR < 0.1){
	  hNAllTauJets->Fill(fMT2tree->misc.MT2, weight);
	}

	jetIndex = -1;
	minDR = 100.0;

	//To find the rec/id tau matched to this gen had tau
	for(int k = 0; k < fMT2tree->NTaus; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->tau[k].lv.Eta(), genleptPhi, fMT2tree->tau[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}

	float DeltaPhi_ = Util::DeltaPhi(fMT2tree->tau[jetIndex].lv.Phi(), METPhi);

	if(DeltaPhi_ > maxDPhiTauMET)
	  continue;
	
	if(minDR < 0.1){
	  if(fMT2tree->tau[jetIndex].Isolation3Hits > IsoCut && fMT2tree->tau[jetIndex].MT < 100)
	    hMT2TauFound->Fill(fMT2tree->misc.MT2, weight);
	}
      }}//TopDeayMode
      Filled = 0;
      if(GenTau == 0 && data != 1){
	for(int j = 0; j < fMT2tree->NTaus; j++){
	  if(fMT2tree->tau[j].Isolation3Hits > IsoCut && fMT2tree->tau[j].MT < 100 && Filled == 0){
	    Filled = 1;
	    hMT2JetSig->Fill(fMT2tree->misc.MT2, weight);
	    hMT2JetSigData->Fill(fMT2tree->misc.MT2, weight);
	  }}}
      }//!=susy


      /*
      //QCD dominated region, use these events to find the fake rate
      if(fMT2tree->misc.MT2 < 70){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    
	    if(fMT2tree->jet[j].isTauMatch < 0  || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;
	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCutLoose)
	      hEtaPtJetData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);

	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){  
	      hEtaPtTauFakeData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    }
	  }
	}else{
	  for(int j = 0; j < fMT2tree->NJets; j++){

	    if(fMT2tree->jet[j].isTauMatch < 0  || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;	    

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;

	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCutLoose)
	      hEtaPtJet->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    
	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){ 
	      
	      hEtaPtTauFake->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    }
	  }
	}
      }*/
      
      Filled = 0;

      //Fill The number of Rec/Id taus and tauable jets in the signal region
      if(fMT2tree->misc.MT2 >= MT2Min){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    if(fMT2tree->jet[j].isTauMatch < 0 || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;
 
	    
// 	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCutLoose){
// 	      hMT2JetSigData->Fill(fMT2tree->misc.MT2, weight);
// 	      cout<<"JetSigData:: "<<fMT2tree->misc.MT2<<endl;
// 	    }
	  }
	
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;
	    
	    if(fMT2tree->tau[j].MT < 100 && fMT2tree->tau[j].Isolation3Hits > IsoCut &&  Filled == 0){
	      Filled = 1;
	      hMT2TauSigData->Fill(fMT2tree->misc.MT2, weight);
	      //	      cout<<"TauSigData:: "<<fMT2tree->misc.MT2<<endl;
	    }
	  }
	}else{
	  
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    if(fMT2tree->jet[j].isTauMatch < 0 || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100) 
	      continue;

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;
// 	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCutLoose){
// 	      hMT2JetSig->Fill(fMT2tree->misc.MT2, weight);
// 	      cout<<"JetSig:: "<<fMT2tree->misc.MT2<<endl;
// 	    }
	  }
	
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;
	
	    if(fMT2tree->tau[j].MT < 100 && fMT2tree->tau[j].Isolation3Hits > IsoCut &&  Filled == 0){
	      Filled = 1;
	      hMT2TauSig->Fill(fMT2tree->misc.MT2, weight);
	      //	      cout<<"TauSig:: "<<fMT2tree->misc.MT2<<endl;
	     
	    }
	  }
	}
      }
    }//for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) 
    cout<<" hNAllTaus->Integral() "<<hNAllTaus->Integral()<<endl;
    cout<<" nAllTaus In ThisSample "<<(hNAllTaus->Integral() - nAllTausBeforeThisSample)<<endl;
  }// for(int i = 0; i < fSamples.size(); i++){

 for(int j = 1; j < NumberOfMT2Bins; j++){
   float id = hMT2TauFound  ->GetBinContent(j);
   float jet= hNAllTauJets->GetBinContent(j);
   float gen= hNAllTaus   ->GetBinContent(j);
   float EffTotal = -1.0;
   if(gen != 0) EffTotal = id/gen;
   float Eff = -1.0;
   if(jet != 0) Eff = id/jet;
   float Acc = -1.0;
   if(gen != 0) Acc = jet/gen;

   cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<" Eff = "<<EffTotal<<", Eff * Acc = "<<Eff<<" * "<<Acc<<endl;
 }

 int p = 0;
 hMT2TauFound->Divide(hNAllTaus);
 Cout(p++, hMT2TauFound);

 //To add the syatematic uncertainty
 for(int j = 1; j < NumberOfMT2Bins; j++){
   float val = hMT2TauFound->GetBinContent(j);
   float Err = hMT2TauFound->GetBinError(j);
   hMT2TauFound->SetBinError(j, sqrt(Err * Err + sys * sys * val * val));
 }

 Cout(p++, hEtaPtTauFake);
 Cout(p++, hEtaPtJet);
 float TauFake  = hEtaPtTauFake->Integral();
  float Jet      = hEtaPtJet    ->Integral();
  float Fake     = TauFake/Jet;
  cout<<" Fake "<<Fake<<endl;
  hEtaPtTauFake->Divide(hEtaPtJet);
  Cout(p++, hEtaPtTauFake);

  Cout(p++, hMT2JetSig);

  //  hMT2JetSig->Scale(Fake);

  Cout(p++, hMT2JetSig);
  Cout(p++, hMT2TauSig);
  hMT2JetSig->Add(hMT2TauSig, -1);
  hMT2JetSig->Scale(-1.0);
  Cout(p++, hMT2JetSig);

  hMT2JetSig->Divide(hMT2TauFound);
  Cout(p++, hMT2JetSig);

  float TauFakeData  = hEtaPtTauFakeData->Integral();
  float JetData      = hEtaPtJetData    ->Integral();
  float FakeData     = TauFakeData/JetData;
  cout<<" FakeData "<<FakeData<<endl;
  hEtaPtTauFakeData->Divide(hEtaPtJetData);
  
  Cout(p++, hMT2JetSigData);
  //   hMT2JetSigData->Scale(FakeData);
  
  Cout(p++, hMT2JetSigData);
  Cout(p++, hMT2TauSigData);
  hMT2JetSigData->Add(hMT2TauSigData, -1);
  hMT2JetSigData->Scale(-1.0);
  Cout(p++, hMT2JetSigData);
    

  hMT2JetSigData->Divide(hMT2TauFound);
  
  Cout(p++, hMT2JetSigData);
  

  Cout(p++, hRatioRelaxed);
  Cout(p++, hRatioTight);

  hRatioRelaxed->Add(hRatioTight);
  
  hRatioTight->Divide(hRatioRelaxed);

  Cout(p++, hRatioTight);

//   hMT2JetSig->Multiply(hRatioTight);

//   hMT2JetSigData->Multiply(hRatioTight);

//   hNAllTaus->Multiply(hRatioTight);
 
  for(int j = 1; j < NumberOfMT2Bins; j++){
    
    float TotalErrSim = hMT2JetSig    ->GetBinError(j);
    float StatErrSim  = (hMT2TauSig    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
    float SysErrSim   = sqrt(fabs(TotalErrSim * TotalErrSim - StatErrSim * StatErrSim));

    float TotalErrData = hMT2JetSigData    ->GetBinError(j);
    float StatErrData  = (hMT2TauSigData    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
    float SysErrData   = sqrt(fabs(TotalErrData * TotalErrData - StatErrData * StatErrData));
    
    cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<"$: Sim  "<<hNAllTaus     ->GetBinContent(j)<<" +/- "<<hNAllTaus        ->GetBinError(j)*hRatioTight->GetBinContent(j)<<
                                                 " Pre  "<<hMT2JetSig    ->GetBinContent(j)<<" +/- "<<StatErrSim <<" +/- "<< SysErrSim <<                              
                                                 " Data "<<hMT2JetSigData->GetBinContent(j)<<" +/- "<<StatErrData<<" +/- "<< SysErrData<<endl;
  }

  for(int j = 0; j < NumberOfSamples; j++){
     AddOverAndUnderFlow(MT2[j]);
  }
  printYield();

THStack* h_stack     = new THStack(varname, "");
 double xbin[7] = {0.0, 35.0, 70.0, 110.0, 150.0, 250.0, 2000};
 TH1F *hnew[NumberOfSamples];
  for(int j = 0; j < NumberOfSamples; j++){
    hnew[j] = (TH1F*)MT2[j]->Rebin(6,"hnew",xbin);
    TH1F* mt2 = (TH1F*)MT2[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 2))
      //h_stack  -> Add(MT2[j]);
      h_stack  -> Add(hnew[j]);
    delete mt2;
  }  

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[5], "data", "l");
  
  //  TLegend *Legend1;
  printHisto(h_stack, MT2[5]->Rebin(6,"hnew1",xbin), MT2[4]->Rebin(6,"hnew4",xbin), Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0);

  plotRatioStack(h_stack,  MT2[4]->Rebin(6,"hnew2",xbin), MT2[5]->Rebin(6,"hnew3",xbin), true, false, "MT2_ratio", Legend1, "MT2", "Events", -4, 0);

  TCanvas *CanvasMT2 = new TCanvas("myMT2", "myMT2");
  hNAllTaus      ->SetLineColor(kBlue);
  hMT2JetSig     ->SetLineColor(kRed);
  hMT2JetSigData ->SetLineColor(kGreen);
  hNAllTaus      ->Draw();
  hMT2JetSig     ->Draw("sames");
  hMT2JetSigData ->Draw("sames");
}


void MassPlotter::NewTauContamination(int sample_index, unsigned int nevents, int flag){
  float sys = 0.00; //systematic uncertanty on the fake rate and efficiency
  float IsoCut = 1.5;//cut on the isolation 0.5:VLoose, 1.5:Loose...
  float IsoCutLoose = IsoCut - 1.0;
  float maxDPhiTauMET = 11.5;//To cut on the DeltaPhi between the tau and MET

  setFlags(flag); 

  TH1::SetDefaultSumw2();
  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    MT2[i] = new TH1F(varname+"_"+cnames[i], "", 150, 0, 750);
    MT2[i] -> SetFillColor  (ccolor[i]);
    MT2[i] -> SetLineColor  (ccolor[i]);
    MT2[i] -> SetLineWidth  (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[5] -> SetMarkerStyle(20);
  MT2[5] -> SetMarkerColor(kBlack);
  MT2[5] -> SetLineColor(kBlack);
  
  MT2[4] -> SetFillStyle(3004);
  MT2[4] -> SetFillColor(kBlack);

  TH2F* hEtaPtJet           = new TH2F("hEtaPtJet",           "hEtaPtJet",          1,0,1000.0,1,0,2.5);//tauable jets in MC in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFake       = new TH2F("hEtaPtTauFake",       "hEtaPtTauFake",      1,0,1000.0,1,0,2.5);//taus in MC in QCD dominated region  (MT2 < 70)
  TH2F* hEtaPtJetData       = new TH2F("hEtaPtJetData",       "hEtaPtJetData",      1,0,1000.0,1,0,2.5);//tauable jets in data in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFakeData   = new TH2F("hEtaPtTauFakeData",   "hEtaPtTauFakeData",  1,0,1000.0,1,0,2.5);//taus in data in QCD dominated region  (MT2 < 70)

  TH1F* hMT2TauTightFound        = new TH1F("hMT2TauTightFound",   "hMT2TauTightFound",   NumberOfMT2Bins - 1, xMT2bin);//Rec/IdT taus
  TH1F* hMT2TauLooseFound        = new TH1F("hMT2TauLooseFound",   "hMT2TauLooseFound",   NumberOfMT2Bins - 1, xMT2bin);//Rec/IdL taus

  TH1F* hMT2Tau0Sig         = new TH1F("hMT2Tau0Sig",     "hMT2Tau0Sig",     NumberOfMT2Bins - 1, xMT2bin);//events without tau(MC)
  TH1F* hMT2Tau1Sig         = new TH1F("hMT2Tau1Sig",     "hMT2Tau1Sig",     NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (MC)
  TH1F* hMT2Tau0SigData     = new TH1F("hMT2Tau0SigData", "hMT2Tau0SigData", NumberOfMT2Bins - 1, xMT2bin);//events without tau (data)
  TH1F* hMT2Tau1SigData     = new TH1F("hMT2Tau1SigData", "hMT2Tau1SigData", NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (data)
  TH1F* hRatioTight         = new TH1F("hRatioTight",       "hRatioTight",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut
  TH1F* hRatioRelaxed       = new TH1F("hRatioRelaxed",       "hRatioRelaxed",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut

  TH1F* hNAllTaus           = new TH1F("hNAllTaus",      "hNAllTaus",      NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus
  TH1F* hNAllTauJets        = new TH1F("hNAllTauJets",   "hNAllTauJets",   NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus matched with a rec jet

  TH1F* hMT2TauFake         = new TH1F("hMT2TauFake",   "hMT2TauFake",   NumberOfMT2Bins - 1, xMT2bin);//Fake rate (MC)
  TH1F* hMT2TauFakeData     = new TH1F("hMT2TauFakeData", "hMT2TauFakeData",   NumberOfMT2Bins - 1, xMT2bin);//Fake rate (data)
  TH1F* hMT2Unit            = new TH1F("hMT2Unit",      "hMT2Unit",   NumberOfMT2Bins - 1, xMT2bin);//Fake rate (MC)
 
  for(int ii = 0; ii < (int) fSamples.size(); ii++){

    if(sample_index != -1 && ii != sample_index)
      continue;

    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data")
      data = 1;

    if(Sample.type == "susy")// || Sample.sname == "QCD") 
      continue;

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    unsigned int nentries =  Sample.tree->GetEntries();

    float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents * Sample.PU_avg_weight );

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
 
      float weight = Weight;
      if(data == 0)
	weight *= fMT2tree->pileUp.Weight;
      else
	weight = 1.0;

      if (!(fMT2tree->misc.MET             >  30                                  &&
	    fMT2tree->NBJetsCSVT           > 0                                    &&
            fMT2tree->NJetsIDLoose40       >= 4                                   && 
           (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
 	   (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)                   &&
 	    fMT2tree->misc.Jet0Pass        == 1                                   &&
 	    fMT2tree->misc.Jet1Pass        == 1                                   &&
 	    fMT2tree->misc.PassJetID       == 1                                   &&
	    fMT2tree->misc.Vectorsumpt     <  70                                  &&
	    //          (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	    fMT2tree->misc.MinMetJetDPhi4  >  0.3                                  &&
	    //Six Or Quad
	    ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0)  || 
	     (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0)   || 
	     (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0)   || 
	     (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0))  &&
	    
	    (fMT2tree->misc.ProcessID!=0 || 
	     ((fMT2tree->misc.CrazyHCAL==0          && fMT2tree->misc.NegativeJEC==0   && 
	       fMT2tree->misc.CSCTightHaloIDFlag==0 && fMT2tree->misc.HBHENoiseFlag==0 && 
	       fMT2tree->misc.hcalLaserEventFlag==0 && fMT2tree->misc.trackingFailureFlag==0 && 
	       fMT2tree->misc.eeBadScFlag==0 && fMT2tree->misc.EcalDeadCellTriggerPrimitiveFlag==0 )&&
	      //Six Or Quad
	      ((((fMT2tree->trigger.HLT_QuadJet80_v1 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v2 ==1)||
		 (fMT2tree->trigger.HLT_QuadJet80_v3 ==1)) && 
		((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0 ) || 
		 (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0 )))|| 

	       (((fMT2tree->trigger.HLT_SixJet45_v1 ==1)||(fMT2tree->trigger.HLT_SixJet45_v2 ==1)||
		 (fMT2tree->trigger.HLT_SixJet45_v3 ==1)) && 
		((fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0 ) || 
		 (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0 ))))))))
	continue;

      float METPhi = fMT2tree->misc.METPhi; //To be used if a cone close to met is used to look for taus
	  
      for(int t = 0; t < fMT2tree->NTaus; t++){

 	if(fMT2tree->tau[t].MT > 100 || fMT2tree->tau[t].Isolation < IsoCut)
 	  continue;

 	float DeltaPhi_ = Util::DeltaPhi(fMT2tree->tau[t].lv.Phi(), fMT2tree->misc.METPhi);

	if(DeltaPhi_ > maxDPhiTauMET)
	  continue;

	if(data == 1){
	    
	  MT2[5]->Fill(fMT2tree->misc.MT2, weight);//data
	    
	}else{
	  MT2[4]->Fill(fMT2tree->misc.MT2, weight);
	    
	  if(Sample.sname == "Top")  
	    MT2[3]->Fill(fMT2tree->misc.MT2, weight);
	  else
	    if(Sample.sname == "DY")	
	      MT2[2]->Fill(fMT2tree->misc.MT2, weight);
	    else
	      if(Sample.sname == "Wtolnu")
		MT2[1]->Fill(fMT2tree->misc.MT2, weight);
	      else
		if(Sample.sname == "QCD")
		  MT2[0]->Fill(fMT2tree->misc.MT2, weight);
	  if(fMT2tree->misc.MT2 >= MT2Min)
	    cout<<"Tau:: "<<fMT2tree->misc.MT2<<endl;
	}}

      int Filled = 0;

      for(int j = 0; j < fMT2tree->NGenLepts; j++){
	if(abs(fMT2tree->genlept[j].ID) != 15)
	  continue;
	
	int HadTau = 1;
	    
	for(int i = 0; i < fMT2tree->NGenLepts; i++){
	  if((abs(fMT2tree->genlept[i].ID) == 11 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID) || 
	     (abs(fMT2tree->genlept[i].ID) == 13 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID)){
	    HadTau = 0;
	  }
	}

	if(HadTau == 0)
	  continue;
        
	if(Filled == 0){
	  //Filled = 1;
	  if((mT2 && fMT2tree->misc.MinMetJetDPhi   >  0.3) || (mT2b && fMT2tree->misc.MinMetJetDPhi4   >  0.3))
	    hRatioTight->Fill(fMT2tree->misc.MT2, weight);
	  else
	    hRatioRelaxed->Fill(fMT2tree->misc.MT2, weight);	  
	}
	
	hNAllTaus->Fill(fMT2tree->misc.MT2, weight);
	

	//To find the jet matched to this gen had tau
	int jetIndex = -1;
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();

	float genleptPhi = fMT2tree->genlept[j].lv.Phi();

	for(int k = 0; k < fMT2tree->NJets; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->jet[k].lv.Eta(), genleptPhi, fMT2tree->jet[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}	

	if(minDR < 0.1){
	  hNAllTauJets->Fill(fMT2tree->misc.MT2, weight);
	}

	jetIndex = -1;
	minDR = 100.0;

	//To find the rec/id tau matched to this gen had tau
	for(int k = 0; k < fMT2tree->NTaus; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->tau[k].lv.Eta(), genleptPhi, fMT2tree->tau[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}

	float DeltaPhi_ = Util::DeltaPhi(fMT2tree->tau[jetIndex].lv.Phi(), METPhi);

	if(DeltaPhi_ > maxDPhiTauMET)
	  continue;
	
	if(minDR < 0.1){
	  if(fMT2tree->tau[jetIndex].Isolation > IsoCutLoose && fMT2tree->tau[jetIndex].MT < 100)
	    hMT2TauLooseFound->Fill(fMT2tree->misc.MT2, weight);

	  if(fMT2tree->tau[jetIndex].Isolation > IsoCut && fMT2tree->tau[jetIndex].MT < 100)
	    hMT2TauTightFound->Fill(fMT2tree->misc.MT2, weight);
	}
      }
      
      //QCD dominated region, use these events to find the fake rate
      if(fMT2tree->misc.MT2 < 70){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    
	    if(fMT2tree->jet[j].isTauMatch < 0  || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;
	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCutLoose){
	      hEtaPtJetData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	      
	      if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){  
		hEtaPtTauFakeData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	      }}}}else{
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    
	    if(fMT2tree->jet[j].isTauMatch < 0  || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;	    

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;

	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCutLoose){
	      hEtaPtJet->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	      
	      if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){ 
		
		hEtaPtTauFake->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	      }}}}}
      
      //Fill The number of Rec/Id taus and tauable jets in the signal region
      if(fMT2tree->misc.MT2 >= MT2Min){
	if(data == 1){
	  int nTausLoose = 0;
	  int nTaus = 0;
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;
	    
	    if(fMT2tree->tau[j].MT < 100){
	      if(fMT2tree->tau[j].Isolation > IsoCutLoose){
		nTausLoose++;
		cout<<"TauSigLooseData:: "<<fMT2tree->misc.MT2<<endl;
		if(fMT2tree->tau[j].Isolation > IsoCut){
		  nTaus++;
		  cout<<"TauSigData:: "<<fMT2tree->misc.MT2<<endl;}}}}
	    
	    if(nTausLoose > 0){
	      if(nTaus == 0)
		//	nTaus=0;
		hMT2Tau0SigData->Fill(fMT2tree->misc.MT2, weight);
	      else
		hMT2Tau1SigData->Fill(fMT2tree->misc.MT2, weight);}	    
	  }else{
	  int nTausLoose = 0;
	  int nTaus = 0;
	  
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi_ = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi_ > maxDPhiTauMET)
	      continue;
	
	    if(fMT2tree->tau[j].MT < 100){
	      if(fMT2tree->tau[j].Isolation > IsoCutLoose){
		nTausLoose++;
		cout<<"TauSigLoose:: "<<fMT2tree->misc.MT2<<endl;
		if(fMT2tree->tau[j].Isolation > IsoCut){
		  nTaus++;
		  cout<<"TauSig:: "<<fMT2tree->misc.MT2<<endl;}}}}
	    
	  if(nTausLoose > 0){
	    if(nTaus == 0){
	      hMT2Tau0Sig->Fill(fMT2tree->misc.MT2, weight);}
	    else
	      hMT2Tau1Sig->Fill(fMT2tree->misc.MT2, weight);}}}
    }//for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) 
  }// for(int i = 0; i < fSamples.size(); i++){


 for(int j = 1; j < NumberOfMT2Bins; j++){
   float idT = hMT2TauTightFound  ->GetBinContent(j);
   float idL = hMT2TauLooseFound  ->GetBinContent(j);
   float jet= hNAllTauJets->GetBinContent(j);
   float gen= hNAllTaus   ->GetBinContent(j);
   float EffTotal = -1.0;
   if(gen != 0) EffTotal = idT/gen;
   float EffLtoT = -1.0;
   if(idL != 0) EffLtoT = idT/idL;
   float Eff = -1.0;
   if(jet != 0) Eff = idL/jet;
   float Acc = -1.0;
   if(gen != 0) Acc = jet/gen;

   cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<" Eff = "<<EffTotal<<", LtoT * Eff * Acc = "<<EffLtoT<<" * "<<Eff<<" * "<<Acc<<endl;
 }

 int p = 0;
 Cout(p++, hMT2TauTightFound);
 Cout(p++, hMT2TauLooseFound);
 Cout(p++, hNAllTauJets);
 Cout(p++, hNAllTaus);
 hMT2TauTightFound->Divide(hMT2TauLooseFound);
 cout<<"hMT2TauTightFound->Divide(hMT2TauLooseFound);"<<endl;
 Cout(p++, hMT2TauTightFound);

 //To add the syatematic uncertainty
 for(int j = 1; j < NumberOfMT2Bins; j++){
   float val = hMT2TauTightFound->GetBinContent(j);
   float Err = hMT2TauTightFound->GetBinError(j);
   hMT2TauTightFound->SetBinError(j, sqrt(Err * Err + sys * sys * val * val));
 }

 Cout(p++, hEtaPtTauFake);
 Cout(p++, hEtaPtJet);
 float TauFake  = hEtaPtTauFake->Integral();
 float Jet      = hEtaPtJet    ->Integral();
 float Fake     = TauFake/Jet;
 cout<<" Fake "<<Fake<<endl;
 hEtaPtTauFake->Divide(hEtaPtJet);
 Cout(p++, hEtaPtTauFake);
 for(int j = 1; j < NumberOfMT2Bins; j++){
   hMT2TauFake->SetBinContent(j,hEtaPtTauFake->GetBinContent(4));
   hMT2TauFake->SetBinError(j,hEtaPtTauFake->GetBinError(4));
   hMT2Unit->SetBinContent(j,1.0);
 }

  Cout(p++, hMT2TauFake);
  Cout(p++, hMT2Unit);

  Cout(p++, hMT2Tau0Sig);

  hMT2Tau0Sig->Multiply(hMT2TauFake);
  cout<<"hMT2Tau0Sig->Multiply(hMT2TauFake);"<<endl;
  Cout(p++, hMT2Tau0Sig);

  hMT2Unit->Add(hMT2TauFake, -1.0);
  cout<<"hMT2Unit->Add(hMT2TauFake, -1.0);"<<endl;
  Cout(p++, hMT2Unit);

  Cout(p++, hMT2Tau1Sig);

  hMT2Tau1Sig->Multiply(hMT2Unit);
  cout<<"hMT2Tau1Sig->Multiply(hMT2Unit);"<<endl;
  Cout(p++, hMT2Tau1Sig);

  hMT2Tau0Sig->Add(hMT2Tau1Sig, -1.0);
  cout<<"hMT2Tau0Sig->Add(hMT2Tau1Sig, -1.0);"<<endl;
  Cout(p++, hMT2Tau0Sig);

  hMT2TauFake->Add(hMT2TauTightFound, -1.0);
  cout<<"hMT2TauFake->Add(hMT2TauTightFound, -1.0);"<<endl;
  Cout(p++, hMT2TauFake);

  hMT2Tau0Sig->Divide(hMT2TauFake);
  cout<<"hMT2Tau0Sig->Divide(hMT2TauFake);"<<endl;
  Cout(p++, hMT2Tau0Sig);

  Cout(p++, hNAllTaus);

  hMT2TauLooseFound->Divide(hNAllTaus);
  cout<<"hMT2TauLooseFound->Divide(hNAllTaus);"<<endl;
  Cout(p++, hMT2TauLooseFound);
  
  hMT2Tau0Sig->Divide(hMT2TauLooseFound);
  cout<<"hMT2Tau0Sig->Divide(hMT2TauLooseFound);"<<endl;
  Cout(p++, hMT2Tau0Sig);



  float TauFakeData  = hEtaPtTauFakeData->Integral();
  float JetData      = hEtaPtJetData    ->Integral();
  float FakeData     = TauFakeData/JetData;
  Cout(p++, hEtaPtTauFakeData);
  Cout(p++, hEtaPtJetData);

  cout<<" FakeData "<<FakeData<<endl;
  hEtaPtTauFakeData->Divide(hEtaPtJetData);

  for(int j = 1; j < NumberOfMT2Bins; j++){
    hMT2TauFakeData->SetBinContent(j,hEtaPtTauFakeData->GetBinContent(4));
    hMT2TauFakeData->SetBinError(j,hEtaPtTauFakeData->GetBinError(4));
    hMT2Unit->SetBinContent(j,1.0);
    hMT2Unit->SetBinError(j,0.0);
  }
  Cout(p++, hMT2TauFakeData);
  Cout(p++, hMT2Unit);

  Cout(p++, hMT2Tau0SigData);

  hMT2Tau0SigData->Multiply(hMT2TauFakeData);
  cout<<"hMT2Tau0SigData->Multiply(hMT2TauFakeData);"<<endl;
  Cout(p++, hMT2Tau0SigData);

  hMT2Unit->Add(hMT2TauFakeData, -1.0);
  cout<<"hMT2Unit->Add(hMT2TauFakeData, -1.0);"<<endl;
  Cout(p++, hMT2Unit);

  Cout(p++, hMT2Tau1SigData);
  
  hMT2Tau1SigData->Multiply(hMT2Unit);
  cout<<"hMT2Tau1SigData->Multiply(hMT2Unit);"<<endl;
  Cout(p++, hMT2Tau1SigData);

  hMT2Tau0SigData->Add(hMT2Tau1SigData, -1.0);
  cout<<"hMT2Tau0SigData->Add(hMT2Tau1SigData, -1.0);"<<endl;
  Cout(p++, hMT2Tau0SigData);

  hMT2TauFakeData->Add(hMT2TauTightFound, -1.0);
  cout<<"hMT2TauFakeData->Add(hMT2TauTightFound, -1.0);"<<endl;
  Cout(p++, hMT2TauFakeData);

  hMT2Tau0SigData->Divide(hMT2TauFakeData);
  cout<<"hMT2Tau0SigData->Divide(hMT2TauFakeData);"<<endl;
  Cout(p++, hMT2Tau0SigData);

  Cout(p++, hMT2TauLooseFound);
  
  hMT2Tau0SigData->Divide(hMT2TauLooseFound);
  cout<<"hMT2Tau0SigData->Divide(hMT2TauLooseFound);"<<endl;
  Cout(p++, hMT2Tau0SigData);


  Cout(p++, hRatioRelaxed);
  Cout(p++, hRatioTight);

  hRatioRelaxed->Add(hRatioTight);
  
  hRatioTight->Divide(hRatioRelaxed);

  Cout(p++, hRatioTight);

//   hMT2Tau0Sig->Multiply(hRatioTight);

//   hMT2Tau0SigData->Multiply(hRatioTight);

//   hNAllTaus->Multiply(hRatioTight);
 
  for(int j = 1; j < NumberOfMT2Bins; j++){
    
    float TotalErrSim = hMT2Tau0Sig    ->GetBinError(j);
//     float StatErrSim  = (hMT2Tau1Sig    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
//     float SysErrSim   = sqrt(fabs(TotalErrSim * TotalErrSim - StatErrSim * StatErrSim));

     float TotalErrData = hMT2Tau0SigData    ->GetBinError(j);
//     float StatErrData  = (hMT2Tau1SigData    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
//     float SysErrData   = sqrt(fabs(TotalErrData * TotalErrData - StatErrData * StatErrData));
    
     cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<"$: Sim  "<<hNAllTaus     ->GetBinContent(j)<<" +/- "<<hNAllTaus        ->GetBinError(j)<<//*hRatioTight->GetBinContent(j)<<
       " Pre  "<<hMT2Tau0Sig    ->GetBinContent(j)<<" +/- "<<TotalErrSim<<//" +/- "<<StatErrSim <<" +/- "<< SysErrSim <<                              
      " Data "<<hMT2Tau0SigData->GetBinContent(j)<<" +/- "<<TotalErrData<<//" +/- "<<StatErrData<<" +/- "<< SysErrData<<
      endl;
  }

  for(int j = 0; j < NumberOfSamples; j++){
     AddOverAndUnderFlow(MT2[j]);
  }
  printYield();

THStack* h_stack     = new THStack(varname, "");
 double xbin[5] = {0.0, 35.0, 70.0, 110.0, 750.0};
 TH1F *hnew[NumberOfSamples];
  for(int j = 0; j < NumberOfSamples; j++){
    hnew[j] = (TH1F*)MT2[j]->Rebin(4,"hnew",xbin);
    TH1F* mt2 = (TH1F*)MT2[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 2))
      //h_stack  -> Add(MT2[j]);
      h_stack  -> Add(hnew[j]);
    delete mt2;
  }  

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[5], "data", "l");
  
  //  TLegend *Legend1;
  printHisto(h_stack, MT2[5]->Rebin(4,"hnew1",xbin), MT2[4]->Rebin(4,"hnew4",xbin), Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0);

  plotRatioStack(h_stack,  MT2[4]->Rebin(4,"hnew2",xbin), MT2[5]->Rebin(4,"hnew3",xbin), true, false, "MT2_ratio", Legend1, "MT2", "Events", -4, 0);

  TCanvas *CanvasMT2 = new TCanvas("myMT2", "myMT2");
  hNAllTaus       ->SetLineColor(kBlue);
  hMT2Tau0Sig     ->SetLineColor(kRed);
  hMT2Tau0SigData ->SetLineColor(kGreen);
  hNAllTaus       ->Draw();
  hMT2Tau0Sig     ->Draw("sames");
  hMT2Tau0SigData ->Draw("sames");
}

void MassPlotter::setFlags(int flag){
    mT2  = 0;
    mT2b = 0;
    High = 0;
    Low  = 0;
    BTagRelaxed = 0;

    if(flag < 10)  mT2 = 1;
    if(flag > 10)  mT2b= 1;
    if(flag == 5 || flag == 15)  Low = 1;
    if(flag == 7 || flag == 17)  High = 1;

    if(flag == 19){
      High = 1;
      BTagRelaxed = 1;}

    if(mT2){
      MT2Min = 150;

      xMT2bin[1] = 200.0;
      xMT2bin[2] = 275.0;
      xMT2bin[3] = 375.0;
      xMT2bin[4] = 500.0;
      xMT2bin[5] = 2000.0;
      
      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 30.0;
      newxMT2bin[2] = 90.0;
      newxMT2bin[3] = 120.0;
      newxMT2bin[4] = 150.0;
      newxMT2bin[5] = 200.0;
      newxMT2bin[6] = 275.0;
      newxMT2bin[7] = 375.0;
      newxMT2bin[8] = 500.0;
      newxMT2bin[9] = 750.0;
    }

    if(mT2b){
      MT2Min = 125;

      xMT2bin[1] = 150.0;
      xMT2bin[2] = 200.0;
      xMT2bin[3] = 300.0;
      xMT2bin[4] = 2000.0;
      xMT2bin[5] = 3000.0;
  
      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 25.0;
      newxMT2bin[2] = 50.0;
      newxMT2bin[3] = 75.0;
      newxMT2bin[4] = 100.0;
      newxMT2bin[5] = 125.0;
      newxMT2bin[6] = 150.0;
      newxMT2bin[7] = 200.0;
      newxMT2bin[8] = 300.0;
      newxMT2bin[9] = 750.0;
    }

    if(mT2b && High){
      MT2Min = 125;

      xMT2bin[1] = 150.0;
      xMT2bin[2] = 180.0;
      xMT2bin[3] = 260.0;
      xMT2bin[4] = 2000.0;
      xMT2bin[5] = 3000.0;

      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 25.0;
      newxMT2bin[2] = 50.0;
      newxMT2bin[3] = 75.0;
      newxMT2bin[4] = 100.0;
      newxMT2bin[5] = 125.0;
      newxMT2bin[6] = 150.0;
      newxMT2bin[7] = 180.0;
      newxMT2bin[8] = 260.0;
      newxMT2bin[9] = 750.0;
    }

    if(flag == 27){
      High = 1;
      MT2Min = 125;
      
      xMT2bin[1] = 1500.0;
      xMT2bin[2] = 2000.0;
      xMT2bin[3] = 2500.0;
      xMT2bin[4] = 20000.0;
      xMT2bin[5] = 30000.0;

      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 25.0;
      newxMT2bin[2] = 50.0;
      newxMT2bin[3] = 75.0;
      newxMT2bin[4] = 125.0;
      newxMT2bin[5] = 1500.0;
      newxMT2bin[6] = 2000.0;
      newxMT2bin[7] = 2500.0;
      newxMT2bin[8] = 20000.0;
      newxMT2bin[9] = 30000.0;
    }

    //For muTau 
    MT2Min = 20.0;

     xMT2bin[0] = MT2Min;
     xMT2bin[1] = 40.0;
     xMT2bin[2] = 50.0;
     xMT2bin[3] = 70.0;
     xMT2bin[4] = 90.0;
     xMT2bin[5] = 300.0;
/*    xMT2bin[0] = MT2Min;
    xMT2bin[1] = 60.0;
    xMT2bin[2] = 80.0;
    xMT2bin[3] = 100.0;
    xMT2bin[4] = 150.0;
    xMT2bin[5] = 300.0;*/

}






void MassPlotter::AddOverAndUnderFlow(TH1 * Histo, bool overflow, bool underflow){
  // Add underflow & overflow bins
  // This failed for older ROOT version when the first(last) bin is empty
  // and there are underflow (overflow) events --- must check whether this 
  // is still the case
  if(underflow){
  Histo->SetBinContent(1,
		       Histo->GetBinContent(0) + Histo->GetBinContent(1));
  Histo->SetBinError(1,
		     sqrt(Histo->GetBinError(0)*Histo->GetBinError(0)+
			  Histo->GetBinError(1)*Histo->GetBinError(1) ));
  Histo->SetBinContent(0, 0.0);
  } if(overflow){
  Histo->SetBinContent(Histo->GetNbinsX(),
		       Histo->GetBinContent(Histo->GetNbinsX()  )+ 
		       Histo->GetBinContent(Histo->GetNbinsX()+1) );
  Histo->SetBinError(Histo->GetNbinsX(),
		     sqrt(Histo->GetBinError(Histo->GetNbinsX()  )*
			  Histo->GetBinError(Histo->GetNbinsX()  )+
			  Histo->GetBinError(Histo->GetNbinsX()+1)*
			  Histo->GetBinError(Histo->GetNbinsX()+1)  ));
  Histo->SetBinContent(Histo->GetNbinsX() + 1, 0.0);
  }
}

void MassPlotter::Cout(int k, TH1 * Histo){
  cout<<k<<" "<<Histo->GetName()<<endl;
  for(int j = 1; j <= Histo->GetNbinsX(); j++)
    cout<<"Bin :: "<<j<<" "<<Histo->GetBinContent(j)<<" +/- "<<Histo->GetBinError(j)<<endl;
}

void MassPlotter::Cout(int k, TH2 * Histo){
  cout<<k<<" "<<Histo->GetName()<<endl;
  int nYbins = Histo->GetNbinsY();
  int nXbins = Histo->GetNbinsX();

  for(int i = 1; i <= nYbins; i++){
    for(int j = 1; j <= nXbins; j++){
      cout<<Histo->GetBinContent((nXbins + 2) * i + j)<<" +/- "<<Histo->GetBinError((nXbins + 2) * i + j)<<endl;
    }
  }
}

void MassPlotter::printYield(){

    cout<<"\n           & "<<MT2[0]->GetName()     << " & "<<MT2[2]->GetName()      <<" & "<<MT2[1]->GetName()      <<" & "<<MT2[3]->GetName()     <<" & "<<MT2[4]->GetName()<<" & "<<MT2[5]->GetName()<<" & "<<MT2[NumberOfSamples-2]->GetName()<<"& "<<MT2[NumberOfSamples]->GetName()     <<"& "<<MT2[NumberOfSamples-1]->GetName()     <<"\\"<<endl;
  
    
    double err;
    
    int binNumber1 = MT2[NumberOfSamples-2]->GetNbinsX();
    MT2[NumberOfSamples-2]->IntegralAndError(0,binNumber1,err);
 
     cout<<" full selection & "<<MT2[0]->Integral()     << " & "<<MT2[2]->Integral()      <<" & "<<MT2[1]->Integral()      <<" & "<<MT2[3]->Integral()     <<" & "<<MT2[4]->Integral()     << " & "<<MT2[5]->Integral()     << " & "<<MT2[NumberOfSamples-2]->Integral()<<" $pm$ "<<err<<" & "<<MT2[NumberOfSamples]->Integral()     <<" & "<<MT2[NumberOfSamples-1]->Integral()     <<"\\"<<endl;


    for(int i = 0; i < NumberOfMT2Bins-1; i++){
      binNumber1     = MT2[NumberOfSamples-2]->FindBin(xMT2bin[i] + 0.1);
      int binNumber2 = MT2[NumberOfSamples-2]->FindBin(xMT2bin[i+1] - 0.1);

      MT2[NumberOfSamples-2]->IntegralAndError(binNumber1, binNumber2,err);
      cout<<"$"<<xMT2bin[i]<<"-"<<xMT2bin[i+1]<<"$ &      "<<MT2[0]->Integral(binNumber1, binNumber2)<< " & "<<MT2[2]->Integral(binNumber1, binNumber2) <<" & "<<MT2[1]->Integral(binNumber1, binNumber2) <<" & "<<MT2[3]->Integral(binNumber1, binNumber2)<<" & "<<MT2[4]->Integral(binNumber1, binNumber2)<<" & "<<MT2[5]->Integral(binNumber1, binNumber2)<<" & "<<MT2[NumberOfSamples-2]->Integral(binNumber1, binNumber2)<<" $pm$ "<<err<<" & "<<MT2[NumberOfSamples]->Integral(binNumber1, binNumber2)<<" & "<<MT2[NumberOfSamples-1]->Integral(binNumber1, binNumber2)<<"\\"<<endl;
    }
}


void MassPlotter::GetGenInfo(MT2GenParticle *T1, MT2GenParticle *T2, MT2GenParticle *W1, MT2GenParticle *W2, MT2GenParticle *B1, MT2GenParticle *B2, MT2GenParticle *U1, MT2GenParticle *U2, MT2GenParticle *D1, MT2GenParticle *D2){
  T1->Reset();
  T2->Reset();
  W1->Reset();
  W2->Reset();
  B1->Reset();
  B2->Reset();
  U1->Reset();
  U2->Reset();
  D1->Reset();
  D2->Reset();
  
  MT2GenParticle *Mo;

  for(int i = 0; i < fMT2tree->NGenParticles; i++){
    // cout<<i<<" Id "<<fMT2tree->genparticle[i].ID<<" Mo "<<fMT2tree->genparticle[i].MIndex<<" GMo "<<fMT2tree->genparticle[i].GMIndex<<endl;  

    Mo  = fMT2tree->genparticle[i].GetMother(fMT2tree->genparticle, fMT2tree->NGenParticles);

    if(abs(fMT2tree->genparticle[i].ID) == 6){
      if(T1->Index == -9)
	*T1 = fMT2tree->genparticle[i];
      else
	*T2 = fMT2tree->genparticle[i];
    }

    if(abs(fMT2tree->genparticle[i].ID) == 24){
      if(Mo->Index == T1->Index)
	*W1 = fMT2tree->genparticle[i];
      else
	*W2 = fMT2tree->genparticle[i];
    }

    if(abs(fMT2tree->genparticle[i].ID) == 5 && abs(Mo->ID) == 6){
      if(Mo->Index == T1->Index)
	*B1 = fMT2tree->genparticle[i];
      else
	*B2 = fMT2tree->genparticle[i];
    }


    if(abs(fMT2tree->genparticle[i].ID) < 6 && abs(Mo->ID) == 24){
      if(Mo->Index == W1->Index){
	if(U1->Index == -9)
	  *U1 =  fMT2tree->genparticle[i];
	else 
	  *D1 =  fMT2tree->genparticle[i];
      }else{
	if(U2->Index == -9)
	  *U2 =  fMT2tree->genparticle[i];
	else 
	  *D2 =  fMT2tree->genparticle[i];
      }}
  }
  // cout<<" T1->ID "<<T1->ID<<endl;
  // cout<<" T2->ID "<<T2->ID<<endl;
  // cout<<" W1->ID "<<W1->ID<<endl;
  // cout<<" W2->ID "<<W2->ID<<endl;
  // cout<<" B1->ID "<<B1->ID<<endl;
  // cout<<" B2->ID "<<B2->ID<<endl;
  // cout<<" U1->ID "<<U1->ID<<endl;
  // cout<<" U2->ID "<<U2->ID<<endl;
  // cout<<" D1->ID "<<D1->ID<<endl;
  // cout<<" D2->ID "<<D2->ID<<endl;
}

float MassPlotter::MatchedJets(MT2GenParticle a, int &b){
  float minDR = 100.0;
  for(int i=0; i<fMT2tree->NJets; ++i){
    if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
    
    if(a.lv.DeltaR(fMT2tree->jet[i].lv) < minDR){
      minDR = a.lv.DeltaR(fMT2tree->jet[i].lv);
      b = i;
    }
  }
    return minDR;
}

void MassPlotter::InitializeTopSearch(TopSearch * myTopSearch, int nPart, TLorentzVector lv, int bjet){
	myTopSearch->SetEvt(nPart, 0, 1); 
	myTopSearch->SetEvt(nPart, 1, lv.Px()); 
	myTopSearch->SetEvt(nPart, 2, lv.Py()); 
	myTopSearch->SetEvt(nPart, 3, lv.Pz()); 
	myTopSearch->SetEvt(nPart, 4, lv.E()); 
	myTopSearch->SetEvt(nPart, 5, lv.M()); 
	myTopSearch->SetEvt(nPart, 6, lv.Pt()); 
	myTopSearch->SetEvt(nPart, 7, lv.Eta()); 
	myTopSearch->SetEvt(nPart, 8, lv.Phi()); 
	myTopSearch->SetEvt(nPart, 9, bjet);
}

void MassPlotter::InitializeHemisphereFinder(vector<float> &px, vector<float> &py, vector<float> &pz, vector<float> &E, TLorentzVector lv){
  	px.push_back(lv.Px());
	py.push_back(lv.Py());
	pz.push_back(lv.Pz());
	E .push_back(lv.E());
}

void MassPlotter::TopStudy(TString mySample, unsigned int nevents){

  TH1F* hWMass = new TH1F("W","W",100,0,300);
  TH1F* hTMass = new TH1F("T","T",50,0,300);
  TH2F* hTWMass = new TH2F("TW","TW",50,0,150,50,0,300);
  TH2F* hHemi12 = new TH2F("Hemi12","Hemi12",10,-0.5,9.5,10,-0.5,9.5);
  TH1F* hChi2T = new TH1F("Chi2T", "Chi2T",15,0,5);

  TH1F* hMinDR = new TH1F("MinDR","MinDR",1020,-1,101);
  TH1F* hDPtTop = new TH1F("DPtTop","DPtTop", 10000, -500, 500);

  double x[11] = {0.0, 30.0, 60.0, 90.0, 130.0, 170.0, 250.0, 350.0, 450.0, 1000.0, 1300.0}; 
  double xEta[6] = {0, 0.5, 1.0, 1.5, 2.5, 4.5};
 
  TH1F* hPtTop = new TH1F("PtTop","PtTop", 10, x);
  TH1F* hMT2Top = new TH1F("MT2Top","MT2Top", 10, x);
  TH1F* hEtaTop = new TH1F("EtaTop","EtaTop", 5, xEta);
  TH1F* hPtTopTogether = new TH1F("PtTopTogether","PtTopTogether", 10, x);
  TH1F* hPtTopEff = new TH1F("PtTopEff","PtTopEff", 10, x);
  TH1F* hEtaTopEff = new TH1F("EtaTopEff","EtaTopEff", 5, xEta);
  TH1F* hMT2TopEff = new TH1F("MT2TopEff","MT2TopEff", 10, x);

  TH1F* hPtW = new TH1F("PtW","PtW", 10, x);
  TH1F* hPtWTogether = new TH1F("PtWTogether","PtWTogether", 10, x);

  int NTopsW1b = 0;
  int NTopsW1j = 0;
  int NTopsRightB = 0;
  int NTopsRecTop = 0;


  int NTops0 = 0;

  int Together1 = 0;
  int SplittedW1 = 0;
  int SplittedT1 = 0;
  int NTops1 = 0;

  int Together2 = 0;
  int SplittedW2 = 0;
  int SplittedT2 = 0;
  int NTops2 = 0;

  int together = 0;
  int splittedW = 0;
  int splittedT = 0;
  int nTops = 0;

  int nGenFullHad = 0;
  int TogetherGenFullHad = 0;
  int SplittedBGenFullHad = 0;
  int SplittedWGenFullHad = 0;
  int CorrectTopGenFullHad = 0;

  int TogetherGenTBW = 0;

  int TogetherGenTBQQ = 0;
  int SplittedBGenTBQQ = 0;
  int SplittedWGenTBQQ = 0;

  MT2GenParticle T1, T2, W1, W2, B1, B2, U1, U2, D1, D2;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    sample Sample = fSamples[ii];
    // if(Sample.name != "Run2012A")
    if(Sample.name != mySample)
    //   //  if(Sample.sname != mySample)
      continue;

    // debug level for TopSearch:
    // = 0 if no printout at all
    // = 1 if only final statistics
    // = 2 if also short event per event output
    // = 3 if also longer event per event output
    const int debug = 3;
   
    TopSearch *myTopSearch = new TopSearch(); 
    myTopSearch->SetDebug(debug);
    myTopSearch->SetWithbtag(3); // = 1: btag preferred over chisq, = 2: correct b-jet requested for top
    myTopSearch->SetChisqTmax(2.0);
    myTopSearch->SetChisqWmax(2.0);
    myTopSearch->SetChisqTWlepmax(0.0);
    myTopSearch->SetChisqW1b(1.0);//chisqW1bpenlty
    myTopSearch->SetChisqWlept(1.0);//chisqWleptpenlty
    myTopSearch->SetRescaleW(1);// = 1: rescale W 4-vector to correct W mass, = 0: do not rescale

    //Changed from the default values
    myTopSearch->SetCkWoverlap(0);

    fMT2tree = new MT2tree();

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    unsigned int nentries =  Sample.tree->GetEntries();

    unsigned int minNentriesNevents = min(nentries, nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Events to read: " << minNentriesNevents << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (unsigned int jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 
      together = 0;
      splittedW = 0;
      splittedT = 0;
      nTops = 0;
 
      if ( fVerbose>2 && jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
      // if(fMT2tree->misc.TopDecayMode != 0)
      // 	continue;
      if((fMT2tree->NEles!=0  && fMT2tree->ele[0].lv.Pt() >= 10) || (fMT2tree->NMuons!=0 && fMT2tree->muo[0].lv.Pt() >= 10))
	continue;
      // if(fMT2tree->NJetsIDLoose < 3 )
      //  	continue;
  
     if (!(fMT2tree->misc.MET            >  30                                   &&
          fMT2tree->NBJets               > 0                                     &&
         (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                   &&
         (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                   &&
	  fMT2tree->misc.Jet0Pass        == 1                                    &&
	  fMT2tree->misc.Jet1Pass        == 1                                    &&
	  fMT2tree->misc.LeadingJPt      >  150                                  &&
	  fMT2tree->misc.SecondJPt       >  100                                  &&
	  fMT2tree->misc.PassJetID       == 1                                    &&
	  fMT2tree->misc.Vectorsumpt     <  70                                   &&
         (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	 ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0))
	   ))
	  continue;


     // if(fMT2tree->misc.MT2 < 170.0)
     //   continue;

     TLorentzVector w1, t1, w2, t2;

     // fill Pseudojets with selected objects
     vector<float> px, py, pz, E;
 
     int nPart = 0;
 
     // struct to keep track of indices in hemi association
     struct HemiObj {vector<TString> type; vector<int> index; vector<int> hemisphere;
     } hemiobjs;
      

     GetGenInfo(&T1, &T2, &W1, &W2, &B1, &B2, &U1, &U2, &D1, &D2);

     if(B1.Index == B2.Index){
       //Wjets
       for(int i=0; i<fMT2tree->NJets; ++i){
       if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
       //      cout<<nPart<<" jet: "<<i<<" ";
       //      jet[i].lv.Print();

       InitializeTopSearch(myTopSearch, nPart, fMT2tree->jet[i].lv, fMT2tree->jet[i].IsBJet(3));
       nPart++;

       InitializeHemisphereFinder(px, py, pz, E, fMT2tree->jet[i].lv);

       hemiobjs.index.push_back(i); 
       hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     }
     myTopSearch->SetNpart(nPart);
   
     myTopSearch->GoSearch();

     NTopsRecTop += myTopSearch->NpT();

     for(int i = 0; i < myTopSearch->NpT(); i++){
       TLorentzVector Top = myTopSearch->PT(i);

       hChi2T->Fill(myTopSearch->ChisqT(i));
       hPtTop->Fill(Top.Pt());
       hEtaTop->Fill(fabs(Top.Eta()));
       hMT2Top->Fill(fMT2tree->misc.MT2);

       int Wnumber = myTopSearch->IpT(i,0);
       int W_1 = myTopSearch->IpW(Wnumber, 1);

       if(W_1 == -990)
	 NTopsW1b++;
       else{
	 if(W_1 == -999)
	   NTopsW1j++;
	   
	 if(myTopSearch->heavyJet(i) > 0)
	   NTopsRightB++;
       }
     }
     }//Wjets

     if(U1.Index == D2.Index)      
       continue;
    
     MT2GenParticle TopSemiLep;
//      int TopSemiLepIndex = -1;
     if(abs(U1.Index + D1.Index + U2.Index + D2.Index) > 100){
       if(U1.Index == D1.Index){
	 TopSemiLep = T2;
// 	 if(T1.lv.Pt() > T2.lv.Pt())
// 	   TopSemiLepIndex = 1;
// 	 else
// 	   TopSemiLepIndex = 0;
       }
       else{
	 TopSemiLep = T1;
// 	 if(T1.lv.Pt() > T2.lv.Pt())
// 	   TopSemiLepIndex = 0;
// 	 else
// 	   TopSemiLepIndex = 1;
	
       }
       // cout<<" SemiLep MT2 :: "<<fMT2tree->misc.MT2<<" TopEta "<<TopSemiLep.lv.Eta()<<endl;
       // TopSemiLep.lv.Print();
       hMT2Top->Fill(fMT2tree->misc.MT2);
       hPtTop->Fill(TopSemiLep.lv.Pt());
       hEtaTop->Fill(fabs(TopSemiLep.lv.Eta()));


       hemiobjs.index.clear(); 
       hemiobjs.type.clear();
       hemiobjs.hemisphere.clear();
       nPart = 0;
       myTopSearch->SetNpart(nPart);   
       px.clear();
       py.clear();
       pz.clear();
       E.clear();

       for(int i=0; i<fMT2tree->NJets; ++i){
	 if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
	 //      cout<<nPart<<" jet: "<<i<<" ";
	 //      jet[i].lv.Print();

	 InitializeTopSearch(myTopSearch, nPart, fMT2tree->jet[i].lv, fMT2tree->jet[i].IsBJet(3));
	 nPart++;

	 InitializeHemisphereFinder(px, py, pz, E, fMT2tree->jet[i].lv);

	 hemiobjs.index.push_back(i); 
	 hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
       }

       myTopSearch->SetNpart(nPart);

       myTopSearch->GoSearch();
       
       float minDR = 100.0;
       int MatchedRecTopIndex = -1;

       NTopsRecTop += myTopSearch->NpT();

       for(int i = 0; i < myTopSearch->NpT(); i++){
	 TLorentzVector Top = myTopSearch->PT(i);
	 
	 if(TopSemiLep.lv.DeltaR(Top) < minDR){
	   minDR = TopSemiLep.lv.DeltaR(Top);
	   MatchedRecTopIndex = i;
	 }
       }
       
       hMinDR->Fill(minDR);
       if(minDR < 1.0){
	 int Wnumber = myTopSearch->IpT(MatchedRecTopIndex,0);
	 int W_1 = myTopSearch->IpW(Wnumber, 1);

	 if(W_1 == -990)
	   NTopsW1b++;
	 else{
	   if(W_1 == -999)
	     NTopsW1j++;
	   
	   if(myTopSearch->heavyJet(MatchedRecTopIndex) > 0)
	     NTopsRightB++;
	 }

	 // cout<<"Rec SemiLep MT2 :: "<<fMT2tree->misc.MT2<<" DR "<<minDR<<" MatchedRecTopIndex "<<MatchedRecTopIndex<<endl;
	 // myTopSearch->PT(MatchedRecTopIndex).Print();
	 hDPtTop->Fill((TopSemiLep.lv.Pt() - myTopSearch->PT(MatchedRecTopIndex).Pt())/TopSemiLep.lv.Pt());
	 hMT2TopEff->Fill(fMT2tree->misc.MT2);
	 hPtTopEff->Fill(TopSemiLep.lv.Pt());
	 hEtaTopEff->Fill(fabs(TopSemiLep.lv.Eta()));
       }
     }//SemiLept
     
     int b1, b2, u1, u2, d1, d2;

     float MaxDR = 0.0;
     
     if(B1.Index != -9 && B2.Index != -9 && U1.Index != -9 && U2.Index != -9 && D1.Index != -9 && D2.Index != -9){
       // hMinDR->Fill(MatchedJets(B1, b1));
       // hMinDR->Fill(MatchedJets(B2, b2));
       // hMinDR->Fill(MatchedJets(U1, u1));
       // hMinDR->Fill(MatchedJets(U2, u2));
       // hMinDR->Fill(MatchedJets(D1, d1));
       // hMinDR->Fill(MatchedJets(D2, d2));

       if(MatchedJets(B1, b1) > MaxDR)
	 MaxDR = MatchedJets(B1, b1);
       if(MatchedJets(B2, b2) > MaxDR)
	 MaxDR = MatchedJets(B2, b2);
       if(MatchedJets(U1, u1) > MaxDR)
	 MaxDR = MatchedJets(U1, u1);
       if(MatchedJets(U2, u2) > MaxDR)
	 MaxDR = MatchedJets(U2, u2);
       if(MatchedJets(D1, d1) > MaxDR)
	 MaxDR = MatchedJets(D1, d1);
       if(MatchedJets(D2, d2) > MaxDR)
	 MaxDR = MatchedJets(D2, d2);
     }else{
       // cout<<" B1 "<<B1.Index<<" B2 "<<B2.Index<<" D1 "<<D1.Index<<" D2 "<<D2.Index<<" U1 "<<U1.Index<<" U2 "<<U2.Index<<endl;
       continue;
     }

     // if(fabs(T1.lv.Eta()) > 2.4 || T1.lv.Pt() < 100.0 || fabs(T2.lv.Eta()) > 2.4 || T2.lv.Pt() < 100.0)
     // 	continue;

     // if(!(fabs(B1.lv.Eta()) < 2.4 && fabs(U1.lv.Eta()) < 2.4 && fabs(D1.lv.Eta()) < 2.4  && fabs(B2.lv.Eta()) < 2.4 && fabs(U2.lv.Eta()) < 2.4 && fabs(D2.lv.Eta()) < 2.4 && B1.lv.Pt() > 20.0 && U1.lv.Pt() > 20.0  && D1.lv.Pt() > 20.0  && B2.lv.Pt() > 20.0  && U2.lv.Pt() > 20.0  && D2.lv.Pt() > 20.0))
     //   continue;
     vector<MT2GenParticle> GenTops;

     if(T1.lv.Pt() > T2.lv.Pt() ){
       GenTops.push_back(T1);
       GenTops.push_back(T2);
     }else{
       GenTops.push_back(T2);
       GenTops.push_back(T1);
     }


     nGenFullHad++;
     hMT2Top->Fill(fMT2tree->misc.MT2);
     hMT2Top->Fill(fMT2tree->misc.MT2);
     // cout<<" Full Had MT2 :: "<<fMT2tree->misc.MT2<<" TopEta "<<T1.lv.Eta()<<endl;
     // T1.lv.Print();
     // cout<<" Full Had MT2 :: "<<fMT2tree->misc.MT2<<" TopEta "<<T2.lv.Eta()<<endl;
     // T2.lv.Print();
     hPtTop->Fill(T1.lv.Pt());
     hPtTop->Fill(T2.lv.Pt());
     hEtaTop->Fill(fabs(T1.lv.Eta()));
     hEtaTop->Fill(fabs(T2.lv.Eta()));
     hPtW->Fill(W1.lv.Pt());
     hPtW->Fill(W2.lv.Pt());

     if(T1.ID < T2.ID)
       for(int i = 0; i < fMT2tree->NGenParticles; i++)
	 cout<<i<<" Id "<<fMT2tree->genparticle[i].ID<<" Mo "<<fMT2tree->genparticle[i].MIndex<<" GMo "<<fMT2tree->genparticle[i].GMIndex<<endl; 

     int Overlap = 0;

     if(b1 == b2 || u1 == u2 || d1 == d2 || u1 == d2 || d1 == u2 || u1 == b2 || d1 == b2 || b1 == u2 || b1 == d2)
       Overlap += 1;

     if(b1 == u1 || b1 == d1 || u1 == d1)
       Overlap += 2;

     if(b2 == u2 || b2 == d2 || u2 == d2)
       Overlap += 4;

  
     nPart = 0;
     myTopSearch->SetNpart(nPart);   
  
     InitializeTopSearch(myTopSearch, nPart, B1.lv, 1);
     InitializeHemisphereFinder(px, py, pz, E, B1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeTopSearch(myTopSearch, nPart, U1.lv, 0);
     InitializeHemisphereFinder(px, py, pz, E, U1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
      
     InitializeTopSearch(myTopSearch, nPart, D1.lv, 0);
     InitializeHemisphereFinder(px, py, pz, E, D1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     InitializeTopSearch(myTopSearch, nPart, B2.lv, 1);
     InitializeHemisphereFinder(px, py, pz, E, B2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     InitializeTopSearch(myTopSearch, nPart, U2.lv, 0);
     InitializeHemisphereFinder(px, py, pz, E, U2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
      
     InitializeTopSearch(myTopSearch, nPart, D2.lv, 0);
     InitializeHemisphereFinder(px, py, pz, E, D2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     myTopSearch->SetNpart(nPart);

     myTopSearch->GoSearch();

     // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
     Hemisphere* hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     //      2,1  0.38
     //      1,3  0.33
     vector<int> grouping = hemisp->getGrouping();

     //    TLorentzVector pseudojet1(0.,0.,0.,0.);
     //     TLorentzVector pseudojet2(0.,0.,0.,0.);
	
     //     for(int i=0; i<px.size(); ++i){
     // 	if(grouping[i]==1){
     // 	  pseudojet1.SetPx(pseudojet1.Px() + px[i]);
     // 	  pseudojet1.SetPy(pseudojet1.Py() + py[i]);
     // 	  pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
     // 	  pseudojet1.SetE( pseudojet1.E()  + E[i]);	
     // 	}else if(grouping[i] == 2){
     // 	  pseudojet2.SetPx(pseudojet2.Px() + px[i]);
     // 	  pseudojet2.SetPy(pseudojet2.Py() + py[i]);
     // 	  pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
     // 	  pseudojet2.SetE( pseudojet2.E()  + E[i]);
     // 	}
     //     }

     // // define MET
     // TLorentzVector MET = B1.lv + B2.lv + U1.lv + U2.lv + D1.lv + D2.lv;

     // Float_t MT2 = fMT2tree->CalcMT2(0,0, pseudojet1, pseudojet2, MET);
     // hTMass->Fill(MT2);

 
     for(int i = 0; i < myTopSearch->NpT(); i++){
       // if(myTopSearch->NpT() == 1)
       //   continue;

       TLorentzVector Top = myTopSearch->PT(i);
       int b = myTopSearch->IpT(i,1);
       int Wnumber = myTopSearch->IpT(i,0);
       int W_0 = myTopSearch->IpW(Wnumber, 0);
       int W_1 = myTopSearch->IpW(Wnumber, 1);
       // cout<<" myTopSearch->NpT() "<<myTopSearch->NpT()<<" b "<<b<<" W_0 "<<W_0<<" W_1 "<<W_1<<endl;
       // cout<<" W "<<endl;
       // myTopSearch->PW(i).Print();
       // cout<<" W1.lv "<<endl;
       // W1.lv.Print();
       // cout<<" T2.lv "<<endl;
       // W2.lv.Print();
	
       if((b == 0 && fabs((myTopSearch->PW(i) - W1.lv).M()) < 1 && fabs((myTopSearch->PW(i) - W1.lv).E()) < 1 ) || 
	  (b == 3 && fabs((myTopSearch->PW(i) - W2.lv).M()) < 1 && fabs((myTopSearch->PW(i) - W2.lv).E()) < 1 ))
	 CorrectTopGenFullHad++;
       else
	 if((b + W_0 + W_1) == 3 || (b + W_0 + W_1) == 12 )
	   CorrectTopGenFullHad++;
	   

       //   else{
       //     //	    if(!((Top - T1.lv).M() < 0.1 || (Top - T2.lv).M() < 0.1)){
       //   if(!((b + W_0 + W_1) == 3 || (b + W_0 + W_1) == 12 )){
       //     cout<<" myTopSearch->NpT() "<<myTopSearch->NpT()<<" b "<<b<<" W_0 "<<W_0<<" W_1 "<<W_1<<endl;
       //     cout<<" Top "<<Top.M()<<endl;
       //     Top.Print();
       //     cout<<" T1.lv "<<T1.lv.M()<<endl;
       //     T1.lv.Print();
       //     cout<<" T2.lv "<<T2.lv.M()<<endl;
       //     T2.lv.Print();
       //     cout<<" W "<<myTopSearch->PW(i).M()<<endl;
       //     myTopSearch->PW(i).Print();
       //     cout<<" W1.lv "<<W1.lv.M()<<endl;
       //     W1.lv.Print();
       //     cout<<" W2.lv "<<W2.lv.M()<<endl;
       //     W2.lv.Print();
       //    }
       // }
     }

     delete hemisp;
     int B1Hemi = grouping[0];
     int B2Hemi = grouping[3];
     if(B1Hemi == grouping[1] && B1Hemi == grouping[2]){
       TogetherGenFullHad++;
       // hPtTopTogether->Fill(T1.lv.Pt());
       // hPtWTogether->Fill(W1.lv.Pt());
     }


     if(B2Hemi == grouping[4] && B2Hemi == grouping[5]){
       TogetherGenFullHad++;
       // hPtTopTogether->Fill(T2.lv.Pt());
       // hPtWTogether->Fill(W2.lv.Pt());
     }

     if(B1Hemi != grouping[1] && B1Hemi != grouping[2])
       SplittedBGenFullHad++;

     if(B2Hemi != grouping[4] && B2Hemi != grouping[5])
       SplittedBGenFullHad++;

     if(grouping[1] != grouping[2])
       SplittedWGenFullHad++;

     if(grouping[4] != grouping[5])
       SplittedWGenFullHad++;

     if(Overlap == 0 && MaxDR < 0.5){
       //if((b_0 == b1 && ((W_0 == u1 && W_1 == d1) || (W_0 == d1 && W_1 == u1))) || (b0 == b2 && ((W_0 == u2 && W_1 == d2) || (W_0 == d2 && W_1 == u2))))

       w1 = fMT2tree->jet[u1].lv + fMT2tree->jet[d1].lv;
       t1 = fMT2tree->jet[b1].lv + w1;

       w2 = fMT2tree->jet[u2].lv + fMT2tree->jet[d2].lv;
       t2 = fMT2tree->jet[b2].lv + w2;

       // hWMass->Fill(w1.M());
       // hWMass->Fill(w2.M());
       // hTMass->Fill(t1.M());
       // hTMass->Fill(t2.M());
       // hTWMass->Fill(w1.M(), t1.M());
       // hTWMass->Fill(w2.M(), t2.M());
     }

     // struct to keep track of indices in hemi association
     // struct HemiObj {vector<TString> type; vector<int> index; vector<int> hemisphere;
     //} hemiobjs;
     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     nPart = 0;
     InitializeHemisphereFinder(px, py, pz, E, B2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeHemisphereFinder(px, py, pz, E, W2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, T1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("Top"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;

     B1Hemi = grouping[0];
     if(B1Hemi == grouping[1]){
       TogetherGenTBW++;
       // hPtTopTogether->Fill(T2.lv.Pt());
       // hPtWTogether->Fill(W2.lv.Pt());
     }

     if(B1Hemi == grouping[1] && B1Hemi == grouping[2])
       cout<<"TogetherGenTBW::Everything in one hemi !!!!"<<endl;

     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     nPart = 0;
     InitializeHemisphereFinder(px, py, pz, E, B1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeHemisphereFinder(px, py, pz, E, W1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, T2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("Top"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;

     B1Hemi = grouping[0];
     if(B1Hemi == grouping[1]){
       TogetherGenTBW++;
       // hPtTopTogether->Fill(T1.lv.Pt());
       // hPtWTogether->Fill(W1.lv.Pt());
     }

     if(B1Hemi == grouping[1] && B1Hemi == grouping[2])
       cout<<"TogetherGenTBW::Everything in one hemi !!!!"<<endl;

     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     nPart = 0;
     InitializeHemisphereFinder(px, py, pz, E, B2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeHemisphereFinder(px, py, pz, E, U2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, D2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, T1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("Top"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;

     B1Hemi = grouping[0];
     if(B1Hemi == grouping[1] && B1Hemi == grouping[2]){
       TogetherGenTBQQ++;
       hPtTopTogether->Fill(T2.lv.Pt());
       hPtWTogether->Fill(W2.lv.Pt());
     }else
       if(B1Hemi != grouping[1] && B1Hemi != grouping[2])
	 SplittedBGenTBQQ++;
       else
	 if(grouping[1] != grouping[2])
	   SplittedWGenTBQQ++;
	 else
	   cout<<" B1Hemi  = "<<B1Hemi<<" grouping[1] "<<grouping[1]<<" grouping[2] "<<grouping[2]<<" grouping[3] "<<grouping[3]<<endl;

     if(B1Hemi == grouping[1] && B1Hemi == grouping[2] && B1Hemi == grouping[3])
       cout<<"TogetherGenTBQQ::Everything in one hemi !!!!"<<endl;


     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     nPart = 0;
     InitializeHemisphereFinder(px, py, pz, E, B1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeHemisphereFinder(px, py, pz, E, U1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, D1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, T2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("Top"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;

     B1Hemi = grouping[0];
     if(B1Hemi == grouping[1] && B1Hemi == grouping[2]){
       TogetherGenTBQQ++;
       hPtTopTogether->Fill(T1.lv.Pt());
       hPtWTogether->Fill(W1.lv.Pt());
     }else
       if(B1Hemi != grouping[1] && B1Hemi != grouping[2])
	 SplittedBGenTBQQ++;
       else
	 if(grouping[1] != grouping[2])
	   SplittedWGenTBQQ++;
	 else
	   cout<<" B1Hemi  = "<<B1Hemi<<" grouping[1] "<<grouping[1]<<" grouping[2] "<<grouping[2]<<" grouping[3] "<<grouping[3]<<endl;

     if(B1Hemi == grouping[1] && B1Hemi == grouping[2] && B1Hemi == grouping[3])
       cout<<"TogetherGenTBQQ::Everything in one hemi !!!!"<<endl;


     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     nPart = 0;
     myTopSearch->SetNpart(nPart);   
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     for(int i=0; i<fMT2tree->NJets; ++i){
       if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
       //      cout<<nPart<<" jet: "<<i<<" ";
       //      jet[i].lv.Print();

       InitializeTopSearch(myTopSearch, nPart, fMT2tree->jet[i].lv, fMT2tree->jet[i].IsBJet(3));
       nPart++;

       InitializeHemisphereFinder(px, py, pz, E, fMT2tree->jet[i].lv);

       hemiobjs.index.push_back(i); 
       hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     }

     myTopSearch->SetNpart(nPart);

     myTopSearch->GoSearch();

     // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;
   
     int hemi1 = 0;
     int hemi2 = 0;
     for(int i = 0; i <(int) px.size(); i++){
       if(grouping[i] == 1)
	 hemi1++;
       if(grouping[i] == 2)
	 hemi2++;
     }
     if((hemi1 + hemi2) != (int) px.size())
       cout<<" hemi1 "<<hemi1<<" hemi2 "<<hemi2<<" px.size() "<<px.size()<<endl;
      
     hHemi12->Fill(hemi1, hemi2);

     for(int i = 0; i < myTopSearch->NpT(); i++){
       TLorentzVector Top = myTopSearch->PT(i);


       // if(abs(myTopSearch->heavyJet(i)) > 200 || nTops > 2)
       // if(nTops > 2)
       //   continue;

       int b_0 = myTopSearch->IpT(i,1);
       int WNumber = myTopSearch->IpT(i,0);
       int W_0 = myTopSearch->IpW(WNumber, 0);
       int W_1 = myTopSearch->IpW(WNumber, 1);

       // if(fabs(myTopSearch->Mw(WNumber) - 80.4) > 0.5 || fabs((myTopSearch->PT(i)).M() - 172.0) > 1.0)
       //   continue;

       nTops++;
       // cout<<" grouping size "<<grouping.size()<<" E.size() "<<E.size()<<" fMT2tree->NJets "<<fMT2tree->NJets<<endl;
       int bHemi = grouping[b_0];
       // cout<<" nTops "<<nTops<<endl;
       // cout<<" heavyJet(i) "<<myTopSearch->heavyJet(i)<<" mW "<<myTopSearch->Mw(WNumber)<<" mT "<<(myTopSearch->PT(i)).M()<<endl;
       hWMass->Fill(myTopSearch->Mw(WNumber));
       hTMass->Fill((myTopSearch->PT(i)).M());
       hTWMass->Fill(myTopSearch->Mw(WNumber),(myTopSearch->PT(i)).M()); 
       hChi2T->Fill(myTopSearch->ChisqT(i));

       // cout<<" b_0, hemi "<<b_0<<" ,"<<bHemi<<endl;
       // cout<<" W_0, hemi "<<W_0<<" ,"<<grouping[W_0]<<endl;
       // cout<<" W_1, hemi "<<W_1<<" ,"<<grouping[W_1]<<endl;
	
       if(abs(myTopSearch->heavyJet(i)) < 200){
	 if(grouping[W_0]==bHemi && grouping[W_1]==bHemi)
	   together++;
	 else
	   if(grouping[W_0]!=bHemi && grouping[W_1]!=bHemi)
	     splittedT++;
	   else
	     if(grouping[W_0]!= grouping[W_1])
	       splittedW++;
	     else{
	       cout<<"What is this? "<<endl;
	       cout<<"b_0, hemi "<<b_0<<", "<<bHemi<<endl;
	       cout<<"W_0, hemi "<<W_0<<", "<<grouping[W_0]<<endl;
	       cout<<"W_1, hemi "<<W_1<<", "<<grouping[W_1]<<endl;
	     }
       }else{
	 if(grouping[W_0]==bHemi)
	   together++;
	 else
	   splittedT++;
       }
     }
     // for(int i = 0; i < myTopSearch->NpT(); i++){
     //   cout<<"Rec Full Had  "<<myTopSearch->PT(i).Eta()<<endl;}
     NTopsRecTop += myTopSearch->NpT();
     vector<int> MatchedGenTop;
     for(int j = 0; j < 2; j++){
       float minDR = 100.0;
       int MatchedRecTopIndex = -1;

       for(int i = 0; i < myTopSearch->NpT(); i++){
	 
	 if(MatchedGenTop.size() > 0 && i == MatchedGenTop[0])
	   continue;

	 TLorentzVector Top = myTopSearch->PT(i);
	 
	 if(GenTops[j].lv.DeltaR(Top) < minDR){
	   minDR = GenTops[j].lv.DeltaR(Top);
	   MatchedRecTopIndex = i;
	 }
       }

       hMinDR->Fill(minDR);
       if(minDR < 1.0){
	 MatchedGenTop.push_back(MatchedRecTopIndex);
	 int Wnumber = myTopSearch->IpT(MatchedRecTopIndex,0);
	 int W_1 = myTopSearch->IpW(Wnumber, 1);
	 if(W_1 == -990)
	   NTopsW1b++;
	 else{
	   if(W_1 == -999)
	     NTopsW1j++;
	   
	   if(myTopSearch->heavyJet(MatchedRecTopIndex) > 0)
	     NTopsRightB++;
	 }

	 hDPtTop->Fill((GenTops[j].lv.Pt() - myTopSearch->PT(MatchedGenTop[j]).Pt())/GenTops[j].lv.Pt());
	 hMT2TopEff->Fill(fMT2tree->misc.MT2);
	 // cout<<"Rec Full Had MT2 :: "<<fMT2tree->misc.MT2<<" DR "<<minDR<<" MatchedRecTopIndex "<<MatchedRecTopIndex<<endl;
	 // myTopSearch->PT(MatchedGenTop[j]).Print();
	 hPtTopEff->Fill(GenTops[j].lv.Pt());
	 hEtaTopEff->Fill(fabs(GenTops[j].lv.Eta()));
       }
     }
     
     if(nTops == 0)
       NTops0++;

     if(nTops == 1){
       NTops1 += nTops;
       Together1  += together;
       SplittedW1 += splittedW;
       SplittedT1 += splittedT;
     }
   
     if(nTops >= 2){
       NTops2 += nTops;
       Together2  += together;
       SplittedW2 += splittedW;
       SplittedT2 += splittedT;
     }

    }

  }



  cout<<" Total number of events with GenFullHad = "<<nGenFullHad<<endl;
  cout<<"TogetherGenTBW = "<<TogetherGenTBW<<", "<<(TogetherGenTBW/(2.0 * nGenFullHad))<<endl;
  cout<<"==================="<<endl;

  cout<<"TogetherGenTBQQ = "<<TogetherGenTBQQ<<", "<<(TogetherGenTBQQ/(2.0 * nGenFullHad))<<endl;
  cout<<"SplittedTGenTBQQ = "<<SplittedBGenTBQQ<<", "<<(SplittedBGenTBQQ/(2.0 * nGenFullHad))<<endl;
  cout<<"SplittedWGenTBQQ = "<<SplittedWGenTBQQ<<", "<<(SplittedWGenTBQQ/(2.0 * nGenFullHad))<<endl;
  cout<<"==================="<<endl;
 
  cout<<"CorrectTopGenFullHad = "<<CorrectTopGenFullHad<<", "<<(CorrectTopGenFullHad/(2.0 * nGenFullHad))<<endl;

  cout<<"==================="<<endl;
  cout<<"TogetherGenFullHad = "<<TogetherGenFullHad<<", "<<(TogetherGenFullHad/(2.0 * nGenFullHad))<<endl;
  cout<<"SplittedTGenFullHad = "<<SplittedBGenFullHad<<", "<<(SplittedBGenFullHad/(2.0 * nGenFullHad))<<endl;
  cout<<"SplittedWGenFullHad = "<<SplittedWGenFullHad<<", "<<(SplittedWGenFullHad/(2.0 * nGenFullHad))<<endl;
  cout<<"==================="<<endl;

  // cout<<"Total number of events = "<<(NTops0 + NTops1 + NTops2/2)<<endl;  

  cout<<"Events with 0 Top = "<<NTops0<<endl;
  cout<<"==================="<<endl;
  cout<<"Events with 1 Top = "<<NTops1<<endl;
  cout<<"Together = "<<Together1<<", "<<(Together1/(1.0 * NTops1))<<endl;
  cout<<"SplittedT = "<<SplittedT1<<", "<<(SplittedT1/(1.0 * NTops1))<<endl;
  cout<<"SplittedW = "<<SplittedW1<<", "<<(SplittedW1/(1.0 * NTops1))<<endl;
  cout<<"==================="<<endl;
  cout<<"Tops in events with >= 2 Top = "<<NTops2<<endl;
  cout<<"Together = "<<Together2<<", "<<(Together2/(1.0 * NTops2))<<endl;
  cout<<"SplittedT = "<<SplittedT2<<", "<<(SplittedT2/(1.0 * NTops2))<<endl;
  cout<<"SplittedW = "<<SplittedW2<<", "<<(SplittedW2/(1.0 * NTops2))<<endl;


  cout<< " W1b "<<NTopsW1b<< " W1j "<<NTopsW1j<<" RightB "<<(NTopsRightB + NTopsW1b)<<" NTopsRecTop "<<NTopsRecTop<<endl;


  TCanvas *MyC = new TCanvas("MyC","MyC");
  MyC->Divide(4,4);
  MyC->cd(1);
  hWMass->GetXaxis()->SetRangeUser(0.0,300.0);
  hWMass->Draw();
  hTMass->Draw("sames");
  MyC->cd(2);
  hTWMass->Draw();
  MyC->cd(3);
  hHemi12->Draw();
  MyC->cd(4);
  hChi2T->GetXaxis()->SetRangeUser(0.0,5.0);
  hChi2T->Draw();

  MyC->cd(5);
  hPtTop->Draw();
  MyC->cd(6);
  //hPtTopTogether->Divide(hPtTop);
  TGraphAsymmErrors *gPtTopTogether = new TGraphAsymmErrors(hPtTopTogether, hPtTop);
  gPtTopTogether->GetXaxis()->SetTitle("Top Pt (GeV)");
  gPtTopTogether->GetYaxis()->SetTitle("Together Fraction");
  gPtTopTogether->Draw("AP");
  MyC->cd(7);
  //hPtWTogether->Divide(hPtW);
  TGraphAsymmErrors *gPtWTogether = new TGraphAsymmErrors(hPtWTogether, hPtW);
  gPtWTogether->GetXaxis()->SetTitle("W Pt (GeV)");
  gPtWTogether->GetYaxis()->SetTitle("Together Fraction");
  gPtWTogether->Draw("AP");
  
  MyC->cd(8);

  cout<<"Number of Tops = "<<hPtTop->Integral()<<endl;
  cout<<"Number of Rec Tops = "<<hPtTopEff->Integral()<<endl;
  cout<<"Eff = "<<hPtTopEff->Integral()/hPtTop->Integral()<<endl;
  
  //hPtTopEff->Divide(hPtTop);
  TGraphAsymmErrors *gPtTopEff = new TGraphAsymmErrors(hPtTopEff, hPtTop);
  gPtTopEff->GetXaxis()->SetTitle("Top Pt (GeV)");
  gPtTopEff->GetYaxis()->SetTitle("Top Finding Eff");
  gPtTopEff->Draw("AP");

  MyC->cd(9);
  hEtaTop->Draw();
  MyC->cd(10);
  //hEtaTopEff->Divide(hEtaTop);
  TGraphAsymmErrors *gEtaTopEff = new TGraphAsymmErrors(hEtaTopEff, hEtaTop);
  gEtaTopEff->GetXaxis()->SetTitle("Top Eta");
  gEtaTopEff->GetYaxis()->SetTitle("Top Finding Eff");
  gEtaTopEff->Draw("AP");

  MyC->cd(11);
  hMT2Top->Draw();
  MyC->cd(12);
  cout<<"Number of Tops (MT2) = "<<hMT2Top->Integral()<<endl;
  cout<<"Number of Rec Tops = "<<hMT2TopEff->Integral()<<endl;
  cout<<"Eff = "<<hMT2TopEff->Integral()/hMT2Top->Integral()<<endl;

  // hMT2TopEff->Divide(hMT2Top);
  TGraphAsymmErrors *gMT2TopEff = new TGraphAsymmErrors(hMT2TopEff, hMT2Top);
  gMT2TopEff->GetXaxis()->SetTitle("Event MT2");
  gMT2TopEff->GetYaxis()->SetTitle("Top Finding Eff");
  gMT2TopEff->Draw("AP");

  MyC->cd(13);
  hDPtTop->Draw();
  MyC->cd(14);
  hMinDR->Draw();
}



void MassPlotter::MyMakePlot(unsigned int nevents){
  static const int NumberOfSamples = 7;
  TH1F* myHisto[NumberOfSamples];


  TH1::SetDefaultSumw2();
  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "susy", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,      0,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    myHisto[i] = new TH1F(varname+"_"+cnames[i], "", 35, 0, 3.5);
    myHisto[i] -> SetFillColor  (ccolor[i]);
    myHisto[i] -> SetLineColor  (ccolor[i]);
    myHisto[i] -> SetLineWidth  (2);
    myHisto[i] -> SetMarkerColor(ccolor[i]);
    myHisto[i] -> SetStats(false);
  }

  myHisto[6] -> SetMarkerStyle(20);
  myHisto[6] -> SetMarkerColor(kBlack);
  myHisto[6] -> SetLineColor(kBlack);
  
  myHisto[4] -> SetFillStyle(3004);
  myHisto[4] -> SetFillColor(kBlack);
  
  myHisto[5] -> SetLineColor(kBlack);
  myHisto[5] -> SetLineStyle(kDotted);


  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    sample Sample = fSamples[ii];
 
    fMT2tree = new MT2tree();
    int counter = 0;
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    unsigned int nentries =  Sample.tree->GetEntries();
    float weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents * Sample.PU_avg_weight );
    unsigned int minNentriesNevents = min(nentries, nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (unsigned int jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

    //   if(jentry ==64){
// 	//	if(((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0) || (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0) && (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5) ))

// 	cout<<"NEle "<<fMT2tree->NEles<<" pt0 "<<fMT2tree->ele[0].lv.Pt()<<endl;
// 	if(fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                
// 	  cout<<jentry<<endl;
// 	cout<<"nJ "<<fMT2tree->NJetsIDLoose40<<" j4pt "<<fMT2tree->jet[3].lv.Pt()<<" j5pt "<<fMT2tree->jet[4].lv.Pt()<<" j6pt "<<fMT2tree->jet[5].lv.Pt()<<" j7pt "<<fMT2tree->jet[6].lv.Pt()<<endl;

//       }

      if (!(fMT2tree->misc.MET             >  30                                  &&
	    fMT2tree->misc.MT2             >  100                                 &&
	    fMT2tree->NBJetsCSVT           > 0                                     &&
            fMT2tree->NJetsIDLoose40       >= 4                                   && 
           (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
 	   (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)                   &&
 	    fMT2tree->misc.Jet0Pass        == 1                                    &&
 	    fMT2tree->misc.Jet1Pass        == 1                                    &&
// 	    // 	  fMT2tree->misc.SecondJPt       >  100                                  &&
// 	    // 	  fMT2tree->misc.LeadingJPt      >  150                                  &&
 	    fMT2tree->misc.PassJetID       == 1                                    &&
	    fMT2tree->misc.Vectorsumpt     <  70                                   &&
	    //          (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	    fMT2tree->misc.MinMetJetDPhi4  >  0.3                                  &&
	    //Six Or Quad
	   ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0) || (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0))                                  &&
	    //HT Trigger
	    // 	   fMT2tree->misc.HT             > 750                                   && 

	    (fMT2tree->misc.ProcessID!=0 || ((fMT2tree->misc.CrazyHCAL==0 &&     fMT2tree->misc.NegativeJEC==0 && fMT2tree->misc.CSCTightHaloIDFlag==0 && fMT2tree->misc.HBHENoiseFlag==0 && fMT2tree->misc.hcalLaserEventFlag==0 &&  fMT2tree->misc.trackingFailureFlag==0 && fMT2tree->misc.eeBadScFlag==0 && fMT2tree->misc.EcalDeadCellTriggerPrimitiveFlag==0 )&&
	     //Six Or Quad
	     ((((fMT2tree->trigger.HLT_QuadJet80_v1 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v2 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v3 ==1)) && ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0 ) || (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0 )))|| (((fMT2tree->trigger.HLT_SixJet45_v1 ==1)||(fMT2tree->trigger.HLT_SixJet45_v2 ==1)||(fMT2tree->trigger.HLT_SixJet45_v3 ==1)) && ((fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0 ) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0 )))))
	     //HT Trigger
	     //    ((fMT2tree->trigger.HLT_PFHT650_v5 ==1) || (fMT2tree->trigger.HLT_PFHT650_v6 ==1) || (fMT2tree->trigger.HLT_PFHT650_v7 ==1) || (fMT2tree->trigger.HLT_PFHT650_v8 ==1) || (fMT2tree->trigger.HLT_PFHT650_v9 ==1))
	        )
))
	continue;
      cout<<++counter<<" : "<<jentry<<endl;
      double myVar = fMT2tree->misc.MinMetBJetDPhi;

      if(Sample.type == "data"){
	    
	myHisto[6]->Fill(myVar, weight);//data
	    
      }else{
	if(Sample.sname == "SMS")
	  myHisto[5]->Fill(myVar, weight);
	else{
	  myHisto[4]->Fill(myVar, weight);
	    
	  if(Sample.sname == "Top")  
	    myHisto[3]->Fill(myVar, weight);
	  else
	    if(Sample.sname == "DY")	
	      myHisto[2]->Fill(myVar, weight);
	    else
	      if(Sample.name == "WJetsToLNu")
		myHisto[1]->Fill(myVar, weight);
	      else
		if(Sample.sname == "QCD")
		  myHisto[0]->Fill(myVar, weight);
	}
      }
    }}
  for(int j = 0; j < NumberOfSamples; j++){
    AddOverAndUnderFlow(myHisto[j]);
  }
  THStack* h_stack     = new THStack(varname, "");
  for(int j = 0; j < NumberOfSamples; j++){

    TH1F* mt2 = (TH1F*)myHisto[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 3))
      h_stack  -> Add(myHisto[j]);
    delete mt2;
  }
    
  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(myHisto[0], "QCD", "f");
  Legend1->AddEntry(myHisto[1], "W+jets", "f");
  Legend1->AddEntry(myHisto[2], "Z+jets", "f");
  Legend1->AddEntry(myHisto[3], "Top", "f");
  Legend1->AddEntry(myHisto[5], "SUSY", "l");
  Legend1->AddEntry(myHisto[6], "data", "l");
  
  //  TLegend *Legend1;
  printHisto(h_stack, myHisto[6], myHisto[4], Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0, 1);

  plotRatioStack(h_stack,  myHisto[4], myHisto[6], true, false, "_ratio", Legend1, "MT2", "Events", -4, 0);
}
  



void MassPlotter::SortHighMT2(float MT2cut, unsigned int nevents){

  vector<int> RunNo;
  vector<int> LSNo;
  vector<int> EventNo;
  vector<int> counterNo;
  vector<float> vecMT2;
  int counter = 0;
  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    sample Sample = fSamples[ii];
    if(Sample.type != "data")
      continue;

    fMT2tree = new MT2tree();

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    unsigned int nentries =  Sample.tree->GetEntries();

    unsigned int minNentriesNevents = min(nentries, nevents);
 
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Events to read: " << minNentriesNevents << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (unsigned int jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

      if (!(fMT2tree->misc.MET            >  30                                   &&
	    fMT2tree->NBJets               > 0                                     &&
	    (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                   &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                   &&
	    fMT2tree->misc.Jet0Pass        == 1                                    &&
	    fMT2tree->misc.Jet1Pass        == 1                                    &&
	    fMT2tree->misc.SecondJPt       >  100                                  &&
	    fMT2tree->misc.LeadingJPt      >  150                                  &&
	    fMT2tree->misc.PassJetID       == 1                                    &&
	    fMT2tree->misc.Vectorsumpt     <  70                                   &&
	    (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	    //Six Or Quad
	    //((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0))                                  &&
	    //HT Trigger
	    fMT2tree->misc.HT             > 750                                   && 

	    ((  fMT2tree->misc.ProcessID==10 || (fMT2tree->misc.CrazyHCAL==0 &&     fMT2tree->misc.NegativeJEC==0 && fMT2tree->misc.CSCTightHaloIDFlag==0 && fMT2tree->misc.HBHENoiseFlag==0 && fMT2tree->misc.hcalLaserEventFlag==0 &&  fMT2tree->misc.trackingFailureFlag==0 && fMT2tree->misc.eeBadScFlag==0 && fMT2tree->misc.EcalDeadCellTriggerPrimitiveFlag==0 ))&&
	     //Six Or Quad
	     //	   ((((fMT2tree->trigger.HLT_QuadJet80_v1 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v2 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v3 ==1)) && fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0 ) || (((fMT2tree->trigger.HLT_SixJet45_v1 ==1)||(fMT2tree->trigger.HLT_SixJet45_v2 ==1)||(fMT2tree->trigger.HLT_SixJet45_v3 ==1)) && fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0 ))
	     //HT Trigger
	     ((fMT2tree->trigger.HLT_PFHT650_v5 ==1) || (fMT2tree->trigger.HLT_PFHT650_v6 ==1) || (fMT2tree->trigger.HLT_PFHT650_v7 ==1) || (fMT2tree->trigger.HLT_PFHT650_v8 ==1) || (fMT2tree->trigger.HLT_PFHT650_v9 ==1))
	

	     )
	    ))
	continue;

      if(fMT2tree->misc.MT2 > MT2cut){
	RunNo.push_back(fMT2tree->misc.Run);
	LSNo.push_back(fMT2tree->misc.LumiSection);
	EventNo.push_back(fMT2tree->misc.Event);
	counterNo.push_back(counter);
	vecMT2.push_back(fMT2tree->misc.MT2);
	counter++;
      }
    }}
  
  // cout<<" MT2 > "<< MT2cut<<endl;
  // for(int i = 0; i < counter; i++){
  //   cout<<counterNo[i]<<" MT2 = "<< MT2[counterNo[i]]<<" Run::LS::Event "<<RunNo[counterNo[i]]<<"::"<<LSNo[counterNo[i]]<<"::"<<EventNo[counterNo[i]]<<endl;
  // }
    
  counterNo = Util::VSort(counterNo, vecMT2, true);
  cout<<" MT2 > "<< MT2cut<<endl;
  for(int i = 0; i < counter; i++){
    cout<<counterNo[i]<<" MT2 = "<< vecMT2[counterNo[i]]<<" Run::LS::Event "<<RunNo[counterNo[i]]<<"::"<<LSNo[counterNo[i]]<<"::"<<EventNo[counterNo[i]]<<endl;
  }
}


void MassPlotter::Efficiency(TString mySample){
  
  TString s = ";m_{Glu};m_{LSP}";
  int nBins = 50;
  float min = 0;
  float max = 2500;
  
  if(mySample == "SMST2"){
    s = ";m_{Stop};m_{LSP}";
    nBins = 100;
  }

  TH2F *SMSEventsAllCuts4or6 = new TH2F("SMSEventsAllCuts4or6", "AllCuts4or6" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsAllCuts4or6MT2_100DelPhi = new TH2F("SMSEventsAllCuts4or6MT2_100DelPhi", "4or6MT2_100DelPhi" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsAllCutsHT = new TH2F("SMSEventsAllCutsHT", "AllCutsHT" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsAllCutsHTMT2_100DelPhi = new TH2F("SMSEventsAllCutsHTMT2_100DelPhi", "HTMT2_100DelPhi" + s, nBins, min, max, nBins, min, max);

  TH2F *SMSEvents = new TH2F("SMSEvents", "PreSelection" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsHT = new TH2F("SMSEventsHT", "HT" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEvents4or6 = new TH2F("SMSEvents4or6", "4 or 6 Jets" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMT2 = new TH2F("SMSEventsMT2", "MT2" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMT2100 = new TH2F("SMSEventsMT2100", "MT2" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMT2Top = new TH2F("SMSEventsMT2Top", "MT2Top" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsNTop1 = new TH2F("SMSEventsNTop1", "NTop1" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsNTop2 = new TH2F("SMSEventsNTop2", "NTop2" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsLeptonVeto = new TH2F("SMSEventsLeptonVeto", "LeptonVeto" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsJetPass = new TH2F("SMSEventsJetPass", "JetPass" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsSecondJet = new TH2F("SMSEventsSecondJet", "SecondJet" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsLeadingJet = new TH2F("SMSEventsLeadingJet", "LeadingJet" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMetJetDPhi = new TH2F("SMSEventsMetJetDPhi", "MetJetDPhi" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMetJetDPhi5 = new TH2F("SMSEventsMetJetDPhi5", "MetJetDPhi 0.5 " + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMetJetDPhi10 = new TH2F("SMSEventsMetJetDPhi10", "MetJetDPhi 1.0 " + s, nBins, min, max, nBins, min, max);
  TH1F *hDeltaPhiMetTop = new TH1F("DeltaPhiMetTop", "DeltaPhiMetTop", 32, 0, 3.14);
  TH1F *hDeltaPhiMetTopTTBar = new TH1F("DeltaPhiMetTopTTBar", "DeltaPhiMetTopTTBar", 32, 0, 3.14);

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    sample Sample = fSamples[ii];

    if(Sample.name != mySample && Sample.sname != "Top")
      continue;

    unsigned int nentries = Sample.tree->GetEntries();
 
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
 

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    if(Sample.sname == "Top"){
      for (unsigned int jentry=0; jentry<nentries;jentry++) {
	Sample.tree->GetEntry(jentry); 
	if ( jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	}
	if (!(fMT2tree->misc.MET                > 30                                    &&
	      fMT2tree->misc.MT2                > 130                                   &&
	      fMT2tree->NBJetsCSVM              > 0                                     &&
	      (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
	      (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)                   &&
	      fMT2tree->misc.Jet0Pass        == 1                                    &&
	      fMT2tree->misc.Jet1Pass        == 1                                    &&
	      // fMT2tree->misc.SecondJPt       >  100                                  &&
	      // fMT2tree->misc.LeadingJPt      >  150                                  &&
	      fMT2tree->misc.PassJetID       == 1                                      &&
	      //   	     fMT2tree->misc.Vectorsumpt        < 70                                    &&
	      // fMT2tree->misc.MinMetJetDPhi4            >  0.6                                  &&
	      ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0) ||
	       (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0)) &&
	      fMT2tree->NJetsIDLoose40          >=4))
	  continue;

	if(fMT2tree->NTops                   >  0){
	  hDeltaPhiMetTopTTBar->Fill(Util::DeltaPhi(fMT2tree->fitTop[0].lv.Phi(), fMT2tree->misc.METPhi));
	}}
      continue;
    }
    h_SMSEvents    = (TH2F*) Sample.file->Get("h_SMSEvents");
	
    for (unsigned int jentry=0; jentry<nentries;jentry++) {
      Sample.tree->GetEntry(jentry); 

      if ( jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);}

      if (!(fMT2tree->misc.MET                > 30                                    &&
	    fMT2tree->NBJets                  > 0                                     &&
	    fMT2tree->misc.Vectorsumpt        < 70                                    &&
	    fMT2tree->NJetsIDLoose40          >=4))
	continue;
      

      if ((fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                   &&
	  (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                   &&
	  fMT2tree->misc.Jet0Pass        == 1                                    &&
	  fMT2tree->misc.Jet1Pass        == 1                                    &&
	  fMT2tree->misc.SecondJPt       >  100                                  &&
	  fMT2tree->misc.LeadingJPt      >  150                                  &&
	  fMT2tree->misc.PassJetID       == 1                                    
	  ){
	//Six Or Quad
	if((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0)){ 
	  if((fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) && fMT2tree->misc.MT2  >  125)
	    SMSEventsAllCuts4or6->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
      
	  if((fMT2tree->misc.MinMetJetDPhi >0.5||fMT2tree->misc.MinMetJetDPhiIndex>3) && fMT2tree->misc.MT2  >  100)
	    SMSEventsAllCuts4or6MT2_100DelPhi->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
	}
	//HT Trigger
	if(fMT2tree->misc.HT             > 750){
	  if((fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) && fMT2tree->misc.MT2  >  125)
	    SMSEventsAllCutsHT->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
      
	  if((fMT2tree->misc.MinMetJetDPhi >0.5||fMT2tree->misc.MinMetJetDPhiIndex>3) && fMT2tree->misc.MT2  >  100)
	    SMSEventsAllCutsHTMT2_100DelPhi->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
	}
      }

      SMSEvents->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.HT              > 750)
	SMSEventsHT->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.MT2             > 125)
	SMSEventsMT2->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if(fMT2tree->misc.MT2             > 100)
	SMSEventsMT2100->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if((fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10) && (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10))
	SMSEventsLeptonVeto->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if((fMT2tree->misc.Jet0Pass        == 1) && (fMT2tree->misc.Jet1Pass        == 1) && (fMT2tree->misc.PassJetID       == 1 ))
	SMSEventsJetPass->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.SecondJPt       >  100) 
	SMSEventsSecondJet->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.LeadingJPt      >  150)
	SMSEventsLeadingJet->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3)
	SMSEventsMetJetDPhi->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if(fMT2tree->misc.MinMetJetDPhi >0.5||fMT2tree->misc.MinMetJetDPhiIndex>3)
	SMSEventsMetJetDPhi5->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if(fMT2tree->misc.MinMetJetDPhi >1.0||fMT2tree->misc.MinMetJetDPhiIndex>3)
	SMSEventsMetJetDPhi10->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0))
	SMSEvents4or6->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.MT2Top             > 125)
	SMSEventsMT2Top->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if(fMT2tree->NTops                   >  0){
	SMSEventsNTop1->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
    
      }
  
      if(fMT2tree->NTops                   >  1)
	SMSEventsNTop2->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
    

      if (!(fMT2tree->misc.MET                > 30                                    &&
	    fMT2tree->misc.MT2                > 130                                   &&
	    fMT2tree->NBJetsCSVM              > 0                                     &&
	    (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)                   &&
	    fMT2tree->misc.Jet0Pass        == 1                                    &&
	    fMT2tree->misc.Jet1Pass        == 1                                    &&
	    // fMT2tree->misc.SecondJPt       >  100                                  &&
	    // fMT2tree->misc.LeadingJPt      >  150                                  &&
	    fMT2tree->misc.PassJetID       == 1                                      &&
	    //   	     fMT2tree->misc.Vectorsumpt        < 70                                    &&
	    // fMT2tree->misc.MinMetJetDPhi4            >  0.6                                  &&
	    ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0) ||
	     (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0)) &&
	    fMT2tree->NJetsIDLoose40          >=4))
	continue;
  
      if(fMT2tree->NTops                   >  0){
	hDeltaPhiMetTop->Fill(Util::DeltaPhi(fMT2tree->fitTop[0].lv.Phi(), fMT2tree->misc.METPhi));
      }
    }}
  
  if(mySample == "SMST2")
    h_SMSEvents->Rebin2D(5,5);
  else
    h_SMSEvents->Rebin2D(10,10);

  TH2F* PreSelected = (TH2F*)SMSEvents->Clone();
  PreSelected->SetName("PreSelected");
  
  TCanvas *MyC = new TCanvas("MyC","MyC");
  MyC->Divide(4,3);
  MyC->cd(1);
  SMSEvents->Divide(h_SMSEvents);
  SMSEvents->Draw("COLZ");
  MyC->cd(2);
  SMSEventsHT->Divide(PreSelected);
  SMSEventsHT->Draw("COLZ");
  MyC->cd(3);
  SMSEvents4or6->Divide(PreSelected);
  SMSEvents4or6->Draw("COLZ");
  MyC->cd(4);
  SMSEventsMT2->Divide(PreSelected);
  SMSEventsMT2->Draw("COLZ");
  MyC->cd(5);
  // SMSEventsMT2Top->Divide(PreSelected);
  // SMSEventsMT2Top->Draw("COLZ");
  SMSEventsMT2100->Divide(PreSelected);
  SMSEventsMT2100->Draw("COLZ");
  MyC->cd(6);
  SMSEventsSecondJet->Divide(PreSelected);
  SMSEventsSecondJet->Draw("COLZ");
  MyC->cd(7);
  SMSEventsLeadingJet->Divide(PreSelected);
  SMSEventsLeadingJet->Draw("COLZ");
  MyC->cd(8);
  SMSEventsMetJetDPhi->Divide(PreSelected);
  SMSEventsMetJetDPhi->Draw("COLZ");
  MyC->cd(9);
  SMSEventsMetJetDPhi5->Divide(PreSelected);
  SMSEventsMetJetDPhi5->Draw("COLZ");
  MyC->cd(10);
  SMSEventsMetJetDPhi10->Divide(PreSelected);
  SMSEventsMetJetDPhi10->Draw("COLZ");
  // SMSEventsLeptonVeto->Divide(PreSelected);
  // SMSEventsLeptonVeto->Draw("COLZ");
  // SMSEventsJetPass->Divide(PreSelected);
  // SMSEventsJetPass->Draw("COLZ");
  MyC->cd(11);
  SMSEventsNTop1->Divide(PreSelected);
  SMSEventsNTop1->Draw("COLZ");
  MyC->cd(12);
  SMSEventsNTop2->Divide(PreSelected);
  SMSEventsNTop2->Draw("COLZ");

  TCanvas *MyC3 = new TCanvas("MyC3","MyC3");
  hDeltaPhiMetTopTTBar->Scale(1.0/hDeltaPhiMetTopTTBar->Integral());
  hDeltaPhiMetTop->Scale(1.0/hDeltaPhiMetTop->Integral());
  hDeltaPhiMetTopTTBar->SetLineColor(4);
  hDeltaPhiMetTop->Draw();
  hDeltaPhiMetTopTTBar->Draw("sames");

  TCanvas *MyC2 = new TCanvas("MyC2","MyC2");
  MyC2->Divide(2,2);

  MyC2->cd(1);
  SMSEventsAllCuts4or6->Divide(h_SMSEvents);
  SMSEventsAllCuts4or6->Draw("COLZ");

  MyC2->cd(2);
  SMSEventsAllCuts4or6MT2_100DelPhi->Divide(h_SMSEvents);
  SMSEventsAllCuts4or6MT2_100DelPhi->Draw("COLZ");

  MyC2->cd(3);
  SMSEventsAllCutsHT->Divide(h_SMSEvents);
  SMSEventsAllCutsHT->Draw("COLZ");

  MyC2->cd(4);
  SMSEventsAllCutsHTMT2_100DelPhi->Divide(h_SMSEvents);
  SMSEventsAllCutsHTMT2_100DelPhi->Draw("COLZ");
}


void MassPlotter::TopStudy2(TString mySample, unsigned int nevents){
  double x[11] = {0.0, 30.0, 60.0, 90.0, 130.0, 170.0, 250.0, 350.0, 450.0, 1000.0, 1300.0}; 
  TH1F* hPtTop = new TH1F("PtTop","PtTop", 10, x);
  TH1F* hPtTopEff = new TH1F("PtTopEff","PtTopEff", 10, x);

  MT2GenParticle T1, T2, W1, W2, B1, B2, U1, U2, D1, D2;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    sample Sample = fSamples[ii];
    if(Sample.name != mySample)
      //   //  if(Sample.sname != mySample)
      continue;

    // debug level for TopSearch:
    // = 0 if no printout at all
    // = 1 if only final statistics
    // = 2 if also short event per event output
    // = 3 if also longer event per event output
    const int debug = 0;
   
    TopSearch *myTopSearch = new TopSearch(); 
    myTopSearch->SetDebug(debug);
    myTopSearch->SetWithbtag(1); // = 1: btag preferred over chisq, = 2: correct b-jet requested for top
    myTopSearch->SetChisqTmax(2.0);
    myTopSearch->SetChisqWmax(2.0);
    myTopSearch->SetChisqTWlepmax(0.0);
    myTopSearch->SetChisqW1b(1.0);//chisqW1bpenlty
    myTopSearch->SetChisqWlept(1.0);//chisqWleptpenlty
    myTopSearch->SetRescaleW(1);// = 1: rescale W 4-vector to correct W mass, = 0: do not rescale

    fMT2tree = new MT2tree();

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    unsigned int nentries =  Sample.tree->GetEntries();

    unsigned int minNentriesNevents = min(nentries, nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Events to read: " << minNentriesNevents << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (unsigned int jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);}

      if (!(fMT2tree->misc.MET            >  30                                   &&
	    fMT2tree->NBJets               > 0                                     &&
	    (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                   &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                   &&
	    fMT2tree->misc.Jet0Pass        == 1                                    &&
	    fMT2tree->misc.Jet1Pass        == 1                                    &&
	    fMT2tree->misc.LeadingJPt      >  150                                  &&
	    fMT2tree->misc.SecondJPt       >  100                                  &&
	    fMT2tree->misc.PassJetID       == 1                                    &&
	    fMT2tree->misc.Vectorsumpt     <  70                                   &&
	    (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	    ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0))
	    ))
	continue;
      int nPart = 0;
      //Read thegen info and find the particle indices 
      GetGenInfo(&T1, &T2, &W1, &W2, &B1, &B2, &U1, &U2, &D1, &D2);

      //Reject Fully Hadronic (quarks do not have mother)
      if(U1.Index == D2.Index)      
	continue;

      MT2GenParticle TopSemiLep;
      //To select semileptonic ttbar 
      if(abs(U1.Index + D1.Index + U2.Index + D2.Index) > 100){
	if(U1.Index == D1.Index){
	  TopSemiLep = T2;
	}else{
	  TopSemiLep = T1;}

	hPtTop->Fill(TopSemiLep.lv.Pt());

	//clear the memory
	nPart = 0;
	myTopSearch->SetNpart(nPart);   


	for(int i=0; i<fMT2tree->NJets; ++i){
	  if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
	  InitializeTopSearch(myTopSearch, nPart, fMT2tree->jet[i].lv, fMT2tree->jet[i].IsBJet(3));
	  nPart++;
	}

	myTopSearch->SetNpart(nPart);

	myTopSearch->GoSearch();
       
	float minDR = 100.0;
  
	for(int i = 0; i < myTopSearch->NpT(); i++){
	  TLorentzVector Top = myTopSearch->PT(i);
	 
	  if(TopSemiLep.lv.DeltaR(Top) < minDR){
	    minDR = TopSemiLep.lv.DeltaR(Top);
	  }
	}

	if(minDR < 1.0){
	  hPtTopEff->Fill(TopSemiLep.lv.Pt());
	}
      }}}

  TCanvas *MyC = new TCanvas("MyC", "MyC");
  // hPtTopEff->Divide(hPtTop);
  //hPtTopEff->Draw();

  TGraphAsymmErrors *gPtTopEff = new TGraphAsymmErrors(hPtTopEff, hPtTop);
  gPtTopEff->GetXaxis()->SetTitle("Top Pt (GeV)");
  gPtTopEff->GetYaxis()->SetTitle("Top Finding Eff");
  gPtTopEff->Draw("AP");
}

void MassPlotter::vs(unsigned int  nevents, TString cuts, TString trigger){ 
   
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  TH2F *AvsBSignal   = new TH2F("AvsBSignal"  , "AvsBSignal"  , 10,    -1.0, 1.0, 10, -1.5, 1.0);
  TH2F *AvsBSignal2  = new TH2F("AvsBSignal2" , "AvsBSignal2" , 10,    -1.0, 1.0, 10, -1.5, 1.0);
  TH2F *AvsBBkg      = new TH2F("AvsBBkg"     , "AvsBBkg"     , 10,    -1.0, 1.0, 10, -1.5, 1.0);
  TH2F *AvsB2Signal  = new TH2F("AvsB2Signal" , "AvsB2Signal" , 60,    0, 600, 60,    0, 150);
  TH2F *AvsB2Signal2 = new TH2F("AvsB2Signal2", "AvsB2Signal2", 60,    0, 600, 60,    0, 150);
  TH2F *AvsB2Bkg     = new TH2F("AvsB2Bkg"    , "AvsB2Bkg"    , 60,    0, 600, 60,    0, 150);
  TH2F *AvsB3Signal  = new TH2F("AvsB3Signal" , "AvsB3Signal" , 60,    0, 300, 60,    0, 150);
  TH2F *AvsB3Signal2 = new TH2F("AvsB3Signal2", "AvsB3Signal2", 60,    0, 300, 60,    0, 150);
  TH2F *AvsB3Bkg     = new TH2F("AvsB3Bkg"    , "AvsB3Bkg"    , 60,    0, 300, 60,    0, 150);

  for(unsigned int i = 0; i < fSamples.size(); i++){

    float OldIntegral = AvsBBkg->Integral() + AvsBSignal->Integral() + AvsBSignal2->Integral();

    sample Sample = fSamples[i];
   
    TString myCuts = cuts;
    if( Sample.type!="data") continue; //
    myCuts += " && " + trigger;//just for completeness

//     if((Sample.sname != "Wtolnu") && (Sample.type != "susy"))
//       continue;
	   
    fMT2tree = new MT2tree();
 
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    
    float Weight;

    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);
    
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;


    if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){

      Weight = Sample.lumi;
      if(fPUReweight) Weight /= Sample.PU_avg_weight;
       
    }

    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();
    
    for (unsigned int jentry = 0; jentry < min(nentries, nevents); jentry++) {
     
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

      float weight = Weight;
     
      weight *= fMT2tree->SFWeight.BTagCSV40eq0 * fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();
	
      if(fPUReweight)
	weight *= fMT2tree->pileUp.Weight;
	
      if(Sample.sname == "Wtolnu")
	weight *= fMT2tree->muTau[0].GetTauWjetsSF();
     
      cout<<fMT2tree->misc.Run<<":"<<fMT2tree->misc.LumiSection<<":"<<fMT2tree->misc.Event<<endl;

      float Q_MT = 1 - 80 * 80/(2 * fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Pt() * fMT2tree->misc.MET);

      float DeltaPhi_L_MET = DeltaPhi(fMT2tree->pfmet[0].Phi(),fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Phi());

      float CosDeltaPhi_L_MET = TMath::Cos(DeltaPhi_L_MET);
      if(Sample.type == "susy"){
	if((fMT2tree->Susy.MassGlu - fMT2tree->Susy.MassLSP) < 150.0 && (fMT2tree->Susy.MassGlu - fMT2tree->Susy.MassLSP) > 100.0){
	  // 	cout << fMT2tree->DeltaREleEle();
	  AvsBSignal->Fill(CosDeltaPhi_L_MET, Q_MT, weight);
	  AvsB2Signal->Fill(fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].MT + fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT, fMT2tree->muTau[0].GetMT2(), weight);
	  AvsB3Signal->Fill(fMT2tree->muTau[0].GetLV().Pt(), fMT2tree->muTau[0].GetMT2(), weight);
	}
	//     if(Sample.sname == "Wtolnu")
	else   if((fMT2tree->Susy.MassGlu - fMT2tree->Susy.MassLSP) < 400.0 && (fMT2tree->Susy.MassGlu - fMT2tree->Susy.MassLSP) > 300.0){
	  // 	cout << fMT2tree->DeltaREleEle();
	  AvsBSignal2->Fill(CosDeltaPhi_L_MET, Q_MT, weight);
	  AvsB2Signal2->Fill(fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].MT + fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT, fMT2tree->muTau[0].GetMT2(), weight);
	  AvsB3Signal2->Fill(fMT2tree->muTau[0].GetLV().Pt(), fMT2tree->muTau[0].GetMT2(), weight);
	}
      }
      else
	{
	  AvsBBkg->Fill(CosDeltaPhi_L_MET, Q_MT, weight);
	  AvsB2Bkg->Fill(fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].MT + fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT, fMT2tree->muTau[0].GetMT2(), weight);
	  AvsB3Bkg->Fill(fMT2tree->muTau[0].GetLV().Pt(), fMT2tree->muTau[0].GetMT2(), weight);
	}
    }
    
    float NewIntegral = AvsBBkg->Integral() + AvsBSignal->Integral() + AvsBSignal2->Integral();
    cout<<" Integral "<<(NewIntegral - OldIntegral) <<endl;
    
  }

  TCanvas *MyC = new TCanvas("MyC", "MyC");

  MyC->Divide(3,3);
  MyC->SetLogz();
  MyC->cd(1);
  AvsBBkg->Draw("colz");
  MyC->cd(2);
  AvsBSignal->Draw("colz");
  MyC->cd(3);
  AvsBSignal2->Draw("colz");
  MyC->cd(4);
  AvsB2Bkg->Draw("colz");
  MyC->cd(5);
  AvsB2Signal->Draw("colz");
  MyC->cd(6);
  AvsB2Signal2->Draw("colz");
  MyC->cd(7);
  AvsB3Bkg->Draw("colz");
  MyC->cd(8);
  AvsB3Signal->Draw("colz");
  MyC->cd(9);
  AvsB3Signal2->Draw("colz");
  
}




void MassPlotter::SpecialMakePlot(unsigned int nevents, TString cuts, TString trigger){
  TH1::SetDefaultSumw2();

  TH1F *NBJets[NumberOfSamples+1];
  TString  cnames[NumberOfSamples+1] = {"QCD", "Wjets", "Zjets", "Top", "MC", "data", "SMS"};
  int      ccolor[NumberOfSamples+1] = {  401,     417,     419,   855,  500,    632,     0};
  TString varname = "NBJetsCSVT";
  for (int i=0; i<(NumberOfSamples+1); i++){
    NBJets[i] = new TH1F(varname+"_"+cnames[i], "", 4, 0, 4);
    NBJets[i] -> SetFillColor  (ccolor[i]);
    NBJets[i] -> SetLineColor  (ccolor[i]);
    NBJets[i] -> SetLineWidth  (2);
    NBJets[i] -> SetMarkerColor(ccolor[i]);
    NBJets[i] -> SetStats(false);
  }

  NBJets[5] -> SetMarkerStyle(20);
  NBJets[5] -> SetMarkerColor(kBlack);
  NBJets[5] -> SetLineColor(kBlack);
  
  NBJets[4] -> SetFillStyle(3004);
  NBJets[4] -> SetFillColor(kBlack);

  NBJets[6] -> SetLineColor(kBlack);
  NBJets[6] -> SetLineStyle(kDotted);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){

//     if(fSamples[ii].type != "data" || fSamples[ii].name != "Run2012D1F")
//       continue;
    float Weight = fSamples[ii].xsection * fSamples[ii].kfact * fSamples[ii].lumi / (fSamples[ii].nevents);

     std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << fSamples[ii].name << endl;
    cout << "   File:           " << (fSamples[ii].file)->GetName() << endl;
    cout << "   Events:         " << fSamples[ii].nevents  << endl;
    cout << "   Events in tree: " << fSamples[ii].tree->GetEntries() << endl; 
    cout << "   Xsection:       " << fSamples[ii].xsection << endl;
    cout << "   kfactor:        " << fSamples[ii].kfact << endl;
    cout << "   avg PU weight:  " << fSamples[ii].PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    fMT2tree = new MT2tree();
    fSamples[ii].tree->SetBranchAddress("MT2tree", &fMT2tree);

    TString myCuts = cuts;

    if(fSamples[ii].type == "data"){
      myCuts += " && " + trigger;
    }


    fSamples[ii].tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();

    for (unsigned int i=0;i<min(nentries, nevents); i++) {
      if(i > nevents) continue;

      if ( i % 100000 == 0 ){ 	
	fprintf(stdout, "\rProcessed events: %6d of %6d ", i + 1, nentries);
	fflush(stdout);
      }
 
      fSamples[ii].tree->GetEntry(myEvtList->GetEntry(i));
 
     float weight = Weight;
     if(fSamples[ii].type == "data")
       weight = 1.0;
     else{
       if(fSamples[ii].type == "susy")
	 weight *= 1.0;
       else{
	 weight *= 1.0;//(fMT2tree->pileUp.Weight/fSamples[ii].PU_avg_weight);//
	 if(fMT2tree->NBJetsCSVT == 0)
	   weight *=fMT2tree->SFWeight.BTagCSV40eq0;
	 if(fMT2tree->NBJetsCSVT == 1)
	   weight *=fMT2tree->SFWeight.BTagCSV40eq1;
	 if(fMT2tree->NBJetsCSVT == 2)
	   weight *=fMT2tree->SFWeight.BTagCSV40eq2;
	 if(fMT2tree->NBJetsCSVT >= 3)
	   weight *=fMT2tree->SFWeight.BTagCSV40ge3;
       }}
   
     if(fSamples[ii].type == "data"){
       NBJets[5]->Fill(fMT2tree->NBJetsCSVT, weight);//data
     }else{
       if(fSamples[ii].type == "susy")
	 NBJets[6]->Fill(fMT2tree->NBJetsCSVT, weight);
       else
	 NBJets[4]->Fill(fMT2tree->NBJetsCSVT, weight);
	    
       if(fSamples[ii].sname == "Top")  
	 NBJets[3]->Fill(fMT2tree->NBJetsCSVT, weight);
       else
	 if(fSamples[ii].sname == "DY")	
	   NBJets[2]->Fill(fMT2tree->NBJetsCSVT, weight);
	 else
	   if(fSamples[ii].sname == "Wtolnu")
	     NBJets[1]->Fill(fMT2tree->NBJetsCSVT, weight);
	   else
	     if(fSamples[ii].sname == "QCD")
	       NBJets[0]->Fill(fMT2tree->NBJetsCSVT, weight);
     }
     
   
    }
    cout << endl;
  }
  
  THStack* h_stack     = new THStack(varname, "");
  for(int j = 0; j < (NumberOfSamples+ 1); j++){
    AddOverAndUnderFlow(NBJets[j]);

    if(j < (NumberOfSamples - 2))
      h_stack  -> Add(NBJets[j]);
  }  

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(NBJets[0], "QCD", "f");
  Legend1->AddEntry(NBJets[1], "W+jets", "f");
  Legend1->AddEntry(NBJets[2], "Z+jets", "f");
  Legend1->AddEntry(NBJets[3], "Top", "f");
  Legend1->AddEntry(NBJets[5], "data", "l");
  Legend1->AddEntry(NBJets[6], "SMS", "l");

  cout << "------------------------------------"                << endl
       << "QCD Integral:      " << NBJets[0]->Integral()  << endl
       << "W+jets Integral:   " << NBJets[1]->Integral()  << endl
       << "Z+jets Integral:   " << NBJets[2]->Integral()  << endl
       << "Top Integral:      " << NBJets[3]->Integral()  << endl
       << "TOTAL BG:          " << NBJets[4]->Integral()  << endl
       << "SUSY:              " << NBJets[6]->Integral()  << endl
       << "Data:              " << NBJets[5]->Integral()  << endl;
        

  cout<< "SUSY1:              " << NBJets[6]->Integral()  << endl;

  //  TLegend *Legend1;
  printHisto(h_stack, NBJets[5], NBJets[4], NBJets[6],Legend1 , "MTC", "hist", true, varname, "Events", -4, 0, 1,1);
  cout<< "SUSY2:              " << NBJets[6]->Integral()  << endl;

  plotRatioStack(h_stack,  NBJets[4], NBJets[5], NBJets[6],true, false, "MT2_ratio", Legend1, varname, "Events", -4, -10, 0, 1);
  cout<< "SUSY3:              " << NBJets[6]->Integral()  << endl;

}


void MassPlotter::MakeCutFlowTable( std::vector<std::string> all_cuts ){
  bool SetfPUReweight = fPUReweight; // To be able to control the pile weight from the run_MassPlot.C

  for( vector<string>::const_iterator cut = all_cuts.begin(); cut != all_cuts.end() ; cut++){
    cout<<*cut<<endl;
  }

  TH1::SetDefaultSumw2();
  std::vector< TH1* > CutFlowHistos;
  string TotalBkgtitle = "Total Bkg";
  TH1 *hTotalBkg = new TH1D(TotalBkgtitle.c_str(), TotalBkgtitle.c_str(), all_cuts.size() , 0 , all_cuts.size() );
  for(unsigned int ii = 0; ii < fSamples.size(); ii++){

    string hName_ = "h" + string(fSamples[ii].name.Data()) + "_cutflow" ;
    TH1* hCutFlowSample = new TH1D( hName_.c_str() , fSamples[ii].sname  , all_cuts.size() , 0 , all_cuts.size() );

    if (fSamples[ii].name == "QCD-Pt-15-20-MuEnriched"){
      fPUReweight = false;}
    else fPUReweight = SetfPUReweight;

    Double_t weight=0;
    if(fPUReweight) weight = fSamples[ii].xsection * fSamples[ii].kfact * fSamples[ii].lumi / (fSamples[ii].nevents*fSamples[ii].PU_avg_weight);
    else            weight = fSamples[ii].xsection * fSamples[ii].kfact * fSamples[ii].lumi / (fSamples[ii].nevents);
    
    TTree* theTree = fSamples[ii].tree ;
    cout << "|" << endl << fSamples[ii].name ;
    TEventList* latest_event_list = NULL;

    string full_cut = "1 ";
    int cut_number = 0 ; 

    for( vector<string>::const_iterator cut = all_cuts.begin(); cut != all_cuts.end() ; cut++){

      full_cut += ("&&" + *cut);

      TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
      if(full_cut.find("NBJetsCSVM") != string::npos ) btagweight = TString::Format("SFWeight.BTagCSV40eq%d",abs(0));

      TString ChannelSpecificSF = "1.00";
		  
      if(myChannel == "muTau"){
	if(fMuIdSF)
	  ChannelSpecificSF += "*muTau[0].muIdSF";
	if(fMuIsoSF)
	  ChannelSpecificSF += "*muTau[0].muIsoSF";
	if(fMuTrgSF)
	  ChannelSpecificSF += "*muTau[0].muTrgSF";
	if(fTauTrgSF)
	  ChannelSpecificSF += "*muTau[0].tauTrgSF";
	if(fTauWjetsSF && (fSamples[ii].sname == "Wtolnu"))
	  ChannelSpecificSF += "*muTau[0].tauWjetsSF";	  
	if(fTauEnergySF )
	  ChannelSpecificSF += "*muTau[0].tauEnergySF";	  
      }
      else
	if(myChannel == "ee"){
	  if(fE0IdIsoSF)
	    ChannelSpecificSF += "*doubleEle[0].Ele0IdIsoSF";
	  if(fE1IdIsoSF)
	    ChannelSpecificSF += "*doubleEle[0].Ele1IdIsoSF";
	  if(fETrgSF)
	    ChannelSpecificSF += "*doubleEle[0].DiEleTrgSF";
	}
	else if (myChannel == "eleMu"){          
	  if (fMuTrgSFeleMu)
	    ChannelSpecificSF += "*eleMu[0].MuTrgSF";
	  if (fMuIdIsoSF)
	    ChannelSpecificSF += "*eleMu[0].MuIdIsoSF";
	  if (fEleTrgSF)
	    ChannelSpecificSF += "*eleMu[0].EleTrgSF";
	  if ( fEleIdIsoSF )
	    ChannelSpecificSF += "*eleMu[0].EleIdIsoSF";
	}else if(myChannel == "eleTau"){
	  ChannelSpecificSF = "eleTau[0].tauTrgSF * eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF";
	  if( fSamples[ii].sname == "Wtolnu" )
	    ChannelSpecificSF += " * ( ( (eleTau[0].tauWjetsSF > 0) * eleTau[0].tauWjetsSF ) + ( (eleTau[0].tauWjetsSF <= 0) * 1.0) )";
	}
      
      if(fStitching && (fSamples[ii].sname == "Wtolnu" || (fSamples[ii].shapename == "ZJetsToLL" && fSamples[ii].name != "DYToLL_M10To50"))){
	weight = fSamples[ii].lumi;
	if(fPUReweight) weight /= fSamples[ii].PU_avg_weight;
      }
       	  
      TString selection;
      if(fSamples[ii].type != "data"){
	if(fPUReweight && fbSFReWeight) selection = TString::Format("(%.15f*pileUp.Weight*%s*%s) * (%s)",weight, btagweight.Data(), ChannelSpecificSF.Data(), full_cut.c_str());
	else if(fPUReweight)  selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",   weight,                    ChannelSpecificSF.Data(), full_cut.c_str());
	else if(fbSFReWeight) selection = TString::Format("(%.15f*%s*%s) * (%s)",              weight, btagweight.Data(), ChannelSpecificSF.Data(), full_cut.c_str());
	else                  selection = TString::Format("(%.15f*%s) * (%s)",                 weight,                    ChannelSpecificSF.Data(), full_cut.c_str()); 
      }else{
	weight = 1;
	selection = TString::Format("(%.15f) * (%s)", weight,  full_cut.c_str()); 
      }

      string full_weight = selection.Data();

      if( cut_number == 0 )
	theTree->Draw( ">>event_list" ,   full_cut.c_str() );
      else{
	theTree->SetEventList( latest_event_list );
	theTree->Draw( ">>event_list" ,  full_cut.c_str() );
      }
      latest_event_list = (TEventList*)( gDirectory->Get("event_list")->Clone("tmp_event_list") );
      theTree->SetEventList( latest_event_list );
      theTree->Draw( "1>>temp_hist" , full_weight.c_str() );
      TH1* hTemp = (TH1*)( gDirectory->Get("temp_hist") );
      double err_yields;
      double yields = hTemp->IntegralAndError( 0 , hTemp->GetNbinsX() , err_yields );
      hCutFlowSample->SetBinContent( cut_number+1 , yields );
      hCutFlowSample->SetBinError( cut_number+1 , err_yields );
      cout << "|" << yields << "\\pm" << err_yields << flush;
      
      cut_number++;
    }

 
    bool was_found_in_snames = false;
    for( std::vector< TH1* >::iterator sname_histo = CutFlowHistos.begin() ; sname_histo != CutFlowHistos.end() ; sname_histo++){
      if( string((*sname_histo)->GetTitle()) == string(fSamples[ii].sname.Data()) ){
	(*sname_histo)->Add( hCutFlowSample );
	was_found_in_snames = true;
      }
    }
    if( !was_found_in_snames ){
      string hName = "h" + string(fSamples[ii].sname.Data()) + "_cutflow" ;
      gROOT->cd();
      CutFlowHistos.push_back( (TH1*)(hCutFlowSample->Clone( hName.c_str() )) );
    }
    if(!(fSamples[ii].sname.Contains("ata") || fSamples[ii].sname.Contains("SUSY")))
      hTotalBkg->Add(hCutFlowSample);

    delete hCutFlowSample;
  }
  std::vector< TH1* >::iterator TotalBkgPlace = CutFlowHistos.end();
  TotalBkgPlace--;

  CutFlowHistos.insert(TotalBkgPlace, hTotalBkg);

  for( std::vector< TH1* >::iterator sname_histo = CutFlowHistos.begin() ; sname_histo != CutFlowHistos.end() ; sname_histo++){
    cout << endl << (*sname_histo)->GetTitle() << "&" ;
    for( int i = 1 ; i <= (*sname_histo)->GetNbinsX() ; i++){
      cout << (*sname_histo)->GetBinContent(i) << "\\pm" << (*sname_histo)->GetBinError(i) << "&" ;
    }
    cout << endl;
  }

  for(unsigned int i = 1 ; i <=  all_cuts.size(); i++){
    if(i == 1)
      cout <<  "          " << "&" ;
    else   
      cout << (all_cuts[i-1]) << " &" ;
    for( std::vector< TH1* >::iterator sname_histo = CutFlowHistos.begin() ; sname_histo != CutFlowHistos.end() ; sname_histo++){
      if(i == 1)
	cout <<  (*sname_histo)->GetTitle() << "&" ;
      else{
	if((TString) (*sname_histo)->GetTitle() != (TString) TotalBkgtitle.c_str())
	  cout << (*sname_histo)->GetBinContent(i) << "&" ;
	else
	  cout << (*sname_histo)->GetBinContent(i) << "\\pm" << (*sname_histo)->GetBinError(i) << "&" ;
      }}
    cout << endl;
  }
}


double MassPlotter::DeltaPhi(double phi1, double phi2){
  // From cmssw reco::deltaPhi()
  double result = phi1 - phi2;
  while( result >   TMath::Pi() ) result -= TMath::TwoPi();
  while( result <= -TMath::Pi() ) result += TMath::TwoPi();
  return TMath::Abs(result);
}
 

#include "TEfficiency.h"
void MassPlotter::TauEfficiency(TString cuts, unsigned int nevents, TString myfileName , TString SampleSName){
  TH1::SetDefaultSumw2();
  double ptBins[] {20,50,100,200,500};
  double etaBins[] {0 , 1.479 , 2.3 };
  double metBins[] {30 , 70 , 110 , 200 , 300 };
  double mt2Bins[] {0, 15 , 40 , 65 , 90 , 500 };
  
  TEfficiency *hPtEta  = new TEfficiency("hPtEta", "PtEta",  4 , ptBins , 2 , etaBins );
  TEfficiency *hPtMET  = new TEfficiency("hPtMET", "PtMET",  4 , ptBins , 4 , metBins );
  TEfficiency *hPt = new TEfficiency("hPt" , "Pt" , 4 , ptBins ); 
  TEfficiency *hMT2 = new TEfficiency("hMT2" , "MT2" , 5 , mt2Bins ); 
  TEfficiency *hOne = new TEfficiency("hOne" , "One" , 1 , 0 , 2);
  
  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
    TString myCuts = cuts;
    sample Sample = fSamples[ii];
    
    if(Sample.sname != SampleSName )
      continue;
    if( SampleSName == "DY" )
      if(! (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50") )
	continue;

    float Weight = 1.0;
    if( (Sample.sname == "Wtolnu"  || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")) ){
      Weight = (Sample.lumi / Sample.PU_avg_weight);
    }


    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    Sample.tree->SetBranchStatus("*", 0);
    Sample.tree->SetBranchStatus("*NTaus" , 1 );
    Sample.tree->SetBranchStatus("*tau*" , 1 );
    Sample.tree->SetBranchStatus("*misc*MET" , 1 );

    Sample.tree->SetBranchStatus("*SFWeight*BTagCSV40eq0" , 1 );

    Sample.tree->SetBranchStatus("*NGenLepts" , 1 );
    Sample.tree->SetBranchStatus("*genlept*" , 1 );
    Sample.tree->SetBranchStatus("muTau*" , 1 );
    Sample.tree->SetBranchStatus("eleTau*" , 1 );
    Sample.tree->SetBranchStatus("NBJetsCSVM" , 1 );
    Sample.tree->SetBranchStatus("*MET*" , 1 );
    Sample.tree->SetBranchStatus("*MinMetJetDPhiPt40*" , 1 );

    Sample.tree->SetBranchStatus("pileUp*Weight" , 1);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();//Sample.tree->GetEntries();
    
    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));
      
      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

      float weight = Weight;
     
      weight *= (fMT2tree->pileUp.Weight * fMT2tree->SFWeight.BTagCSV40eq0);
 
      int GenTau = 0;
     
      for(int j = 0; j < fMT2tree->NGenLepts; j++){
	if(abs(fMT2tree->genlept[j].ID) != 15)
	  continue;
		

	int HadTau = 1;
	    
	for(int i = 0; i < fMT2tree->NGenLepts; i++){
	  if((abs(fMT2tree->genlept[i].ID) == 11 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID) || 
	     (abs(fMT2tree->genlept[i].ID) == 13 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID)){
	    HadTau = 0;
	  }
	}

 	if(HadTau == 0){
 	  continue;
	}

	GenTau++;

	//To find the jet matched to this gen had tau
	int jetIndex = -1;
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();
	float genleptPhi = fMT2tree->genlept[j].lv.Phi();
	
	for(int k = 0; k < fMT2tree->NTaus; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->tau[k].lv.Eta(), genleptPhi, fMT2tree->tau[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}

	if(minDR < 0.1){
	  if(fMT2tree->tau[jetIndex].PassTau_ElTau == 1){

	    //bool pass = ((fMT2tree->tau[jetIndex].PassTau_MuTau == 1 ) && (fMT2tree->tau[jetIndex].Isolation3Hits <= 1.));
	    bool pass = ((fMT2tree->tau[jetIndex].PassTau_ElTau == 1 ) && (fMT2tree->tau[jetIndex].Isolation3Hits <= 1.));
	    hPtEta->FillWeighted(pass, weight , fMT2tree->tau[jetIndex].lv.Pt() , fabs( fMT2tree->tau[jetIndex].lv.Eta() ) );
	    hPtMET->FillWeighted(pass, weight , fMT2tree->tau[jetIndex].lv.Pt() , fMT2tree->misc.MET );
	    hPt->FillWeighted(pass, weight , fMT2tree->tau[jetIndex].lv.Pt() );
	    hOne->FillWeighted(pass, weight , 1.0 );

	    if( fMT2tree->eleTau[0].GetTauIndex0() == jetIndex && 
		fMT2tree->eleTau[0].GetEleIndex0() >= 0 && 
		fMT2tree->eleTau[0].GetSumCharge() == 0 )
	      hMT2->FillWeighted(pass , weight , fMT2tree->eleTau[0].GetMT2() );
	  }
	}
      }
    }
  }

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_" + SampleSName + "_PRHistos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  hPtEta->Write();
  hPtMET->Write();
  hPt->Write();
  hOne->Write();
  hMT2->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" cuts "<<cuts<<endl;
}

void MassPlotter::TauFakeRate(TString cuts, TString trigger, unsigned int nevents, TString myfileName){
 TH1::SetDefaultSumw2();

 
 TH2F *hPtEtaAll  = new TH2F("hPtEtaAll", "hPtEtaAll",  60, -3.0, 3.0, 1000, 0, 1000);
 TH2F *hPtEtaPass = new TH2F("hPtEtaPass","hPtEtaPass", 60, -3.0, 3.0, 1000, 0, 1000);

 TH2F *hTauPtMuPt  = new TH2F("hTauPtMuPt", "hTauPtMuPt",  40, 20, 220, 40, 20, 220);
 TH2F *hTauPtMET   = new TH2F("hTauPtMET",  "hTauPtMET",   40, 20, 220, 40, 20, 220);
 
 float xBins[7] = {30, 50, 75, 100, 150, 250, 500};

 TH1F *hMETAll   = new TH1F("hMETAll"  ,"hMETAll"  ,6, xBins);
 TH1F *hMETPass  = new TH1F("hMETPass" ,"hMETPass" ,6, xBins);
 TH1F *hMuPtAll  = new TH1F("hMuPtAll" ,"hMuPtAll" ,6, xBins);
 TH1F *hMuPtPass = new TH1F("hMuPtPass","hMuPtPass",6, xBins);
 TH1F *hMT2All   = new TH1F("hMT2All"  ,"hMT2All"  ,6, xBins);
 TH1F *hMT2Pass  = new TH1F("hMT2Pass" ,"hMT2Pass" ,6, xBins);

 for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
   TString myCuts = cuts;
 
   int data = 0;
   sample Sample = fSamples[ii];
    
   if(Sample.type == "data"){
     data = 1;
     myCuts += " && " + trigger;
   } else
     if(Sample.sname != "Wtolnu")
       continue;

//   if(Sample.type != "mc")
//     continue;

   fMT2tree = new MT2tree();
   Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

   float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

   std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry); 
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

      if(Sample.sname == "Wtolnu"){
	int genMuTauStatus = fMT2tree->GenLeptonAnalysis(0,1,1);
	
	if(genMuTauStatus != 613 && genMuTauStatus != 623)
	  continue;       
      }
      
      float weight = Weight;

      if(data == 1)
 	weight = 1.0;
      else{
	if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){
	  weight = Sample.lumi;
	  if(fPUReweight) weight /= Sample.PU_avg_weight;
      	}
      
      weight *= fMT2tree->SFWeight.BTagCSV40eq0 * fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

      if(fPUReweight)
	weight *= fMT2tree->pileUp.Weight;

      if(Sample.sname == "Wtolnu")
	weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	}
     
      int chargeMu = 0;      
      TLorentzVector LeadingMuon;
       
      for(int i=0; i<fMT2tree->NMuons; ++i){ 
//	if(fMT2tree->muo[i].PassMu0_TauMu == 1 && fMT2tree->muo[i].lv.Pt() > 27.0 && chargeMu == 0){
	if(fMT2tree->muo[i].PassMu0_TauMu == 1 && fMT2tree->muo[i].lv.Pt() > 20.0 && chargeMu == 0){
	  LeadingMuon = fMT2tree->muo[i].lv;
	  chargeMu = fMT2tree->muo[i].Charge;
	}
	else
	  if(fMT2tree->muo[i].RejMu_TauTau == 1){
	    chargeMu = 0;
	    break;
	  }
      }
      
      if(chargeMu == 0)
	continue;

      for(int i = 0; i < fMT2tree->NEles; i++){
	if(fMT2tree->ele[i].IDVetoMuTau){
	  chargeMu = 0;
	}
      }
      
      if(chargeMu == 0)
	continue;

      float MaxPtMuTau = -1;

      int TauInd = -1;

      for(int t=0; t<fMT2tree->NTaus; t++){ 
	  

// 	if(fMT2tree->tau[t].PassQCDTau_MuTau == 0)
 	if(fMT2tree->tau[t].PassTau_MuTau == 0)
//	if(fMT2tree->tau[t].PassQCDTau0_TauTau == 0)
	  continue;

	float deltaR = Util::GetDeltaR(fMT2tree->tau[t].lv.Eta(), LeadingMuon.Eta(), fMT2tree->tau[t].lv.Phi(), LeadingMuon.Phi());
	
	if(deltaR < 0.2)
	  continue;

	float maxPtMuTau = fMT2tree->tau[t].lv.Pt() + LeadingMuon.Pt();

	if(maxPtMuTau > MaxPtMuTau){
	  MaxPtMuTau = maxPtMuTau;
	  TauInd = t;
	}
      }

      if(TauInd == -1)
	continue;
      
      if(chargeMu == fMT2tree->tau[TauInd].Charge)
	 continue;

      float myMT2 = fMT2tree->CalcMT2(0, false, fMT2tree->tau[TauInd].lv, LeadingMuon, fMT2tree->pfmet[0]);

//       if(fMT2tree->tau[TauInd].PassQCDTau_MuTau == 1){
       if(fMT2tree->tau[TauInd].PassTau_MuTau == 1){
//      if(fMT2tree->tau[TauInd].PassQCDTau0_TauTau == 1){
	hPtEtaAll->Fill(fMT2tree->tau[TauInd].lv.Eta(), fMT2tree->tau[TauInd].lv.Pt(), weight); 
	hMETAll->Fill(fMT2tree->misc.MET, weight);
	hMT2All->Fill(myMT2, weight);
	hMuPtAll->Fill(LeadingMuon.Pt(), weight);
      }
      
       if(fMT2tree->tau[TauInd].PassTau_MuTau == 1 && fMT2tree->tau[TauInd].Isolation3Hits == 1.){
//      if(fMT2tree->tau[TauInd].PassQCDTau0_TauTau == 1 && fMT2tree->tau[TauInd].Isolation3Hits == 1.){
	hPtEtaPass->Fill(fMT2tree->tau[TauInd].lv.Eta(), fMT2tree->tau[TauInd].lv.Pt(), weight); 
	hMETPass->Fill(fMT2tree->misc.MET, weight);
	hMT2Pass->Fill(myMT2, weight);
	hMuPtPass->Fill(LeadingMuon.Pt(), weight);      
	hTauPtMuPt->Fill(fMT2tree->tau[TauInd].lv.Pt(), LeadingMuon.Pt());
	hTauPtMET ->Fill(fMT2tree->tau[TauInd].lv.Pt(), fMT2tree->misc.MET);
      }
    }
 }
 //First save the files, then plot the histos, otherwise the result is lost when run in screen
 TString fileName = fOutputDir;
 if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_FRHistos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  hPtEtaAll->Write();
  hPtEtaPass->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

 TCanvas *MyC = new TCanvas("Fake","Fake");
 MyC->Divide(4,2);
 MyC->cd(1);
 hMuPtAll->Draw();
 MyC->cd(2);
 hMETAll->Draw();
 MyC->cd(3);
 hMT2All->Draw();
 MyC->cd(4);
 hTauPtMuPt->Draw();
 MyC->cd(5);
 hMuPtPass->Divide(hMuPtAll);
 hMuPtPass->Draw();
 MyC->cd(6);
 hMETPass->Divide(hMETAll);
 hMETPass->Draw();
 MyC->cd(7);
 hMT2Pass->Divide(hMT2All);
 hMT2Pass->Draw();
 MyC->cd(8);
 hTauPtMET->Draw();
}



void MassPlotter::muTauAnalysis(TString cuts, TString trigger, unsigned int nevents, TString myfileName, int type){

  TH1::SetDefaultSumw2();
  setFlags(10);

  TString  cnames[NumberOfSamples+1] = {"QCD",     "W",    "ZX",  "Top",       "WW",  "Higgs",   "MC", "susy","data"};
  int      ccolor[NumberOfSamples+1] = { 401,      417,     419,   855,         603,     kRed,    603,      1, 632};
  TString varname = "MT2";
  for (int i=0; i<=(NumberOfSamples); i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i], "", 50, 0, 5);
    //following to check the promptness and faking of the events
    //MT2[i] = new TH1D(varname+"_"+cnames[i], "", 8, 0, 8);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);

    //following to check the promptness and faking of the events: 
//     TString  genStatus[8] = {"pp", "pf", "fp", "ff", "tp", "tf", "nothing", "wrong"};
//     for(int k = 0; k < MT2[i]->GetNbinsX(); k++)
//       MT2[i] ->GetXaxis()->SetBinLabel(k+1,genStatus[k]);
    
}

  MT2[NumberOfSamples] -> SetMarkerStyle(20);
  MT2[NumberOfSamples] -> SetMarkerColor(kBlack);
  MT2[NumberOfSamples] -> SetLineColor(kBlack);
  
  MT2[NumberOfSamples-2] -> SetFillStyle(3004);
  MT2[NumberOfSamples-2] -> SetFillColor(kBlack);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data"){
      data = 1;
      myCuts += " && " + trigger;
    }

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    
    float Weight;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry);
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

      if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){

	Weight = Sample.lumi;
	if(fPUReweight) Weight /= Sample.PU_avg_weight;
      
      }
 
      float weight = Weight;

      if(data == 1)
	weight = 1.0;
      else{
	weight *= fMT2tree->SFWeight.BTagCSV40eq0 * fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();
	
	if(fPUReweight)
	  weight *= fMT2tree->pileUp.Weight;
	
	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();
      }
      
      //float myQuantity = fMT2tree->muTau[0].GetMT2();
        
 //      int nExtraTaus = 0;
//       int tauIndex = -1;//fMT2tree->muTau[0].GetTauIndex0();


//       for(int i=0; i<fMT2tree->NTaus; ++i){
// 	if(fMT2tree->tau[i].PassQCDTau_MuTau && fMT2tree->tau[i].Isolation > 0 && i != fMT2tree->muTau[0].GetTauIndex0())//
// 	  nExtraTaus++;
//       }
        cout<<fMT2tree->misc.Run<<" "<<fMT2tree->misc.LumiSection<<" "<< fMT2tree->misc.Event<<" "<<fMT2tree->muTau[0].GetMT2() <<" " <<fMT2tree->muTau[0].GetLV().M()<<" "<<fMT2tree->muTau[0].GetLV().Pt()<<" "<<fMT2tree->muTau[0].GetLV().Eta()<<" "<<fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pt()<<" "<<fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Eta()<<" "<<fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Pt()<<" "<<" "<<fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Eta()<<" "<<fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT<<" "<<fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].MT<<" "<<fMT2tree->muTau[0].GetDPhi()<<" "<<fMT2tree->misc.MinMetJetDPhiPt40<<endl;

      
      int chargeMu = 0;      
      TLorentzVector LeadingMuon;
       
      for(int i=0; i<fMT2tree->NMuons; ++i){ 
	if(fMT2tree->muo[i].PassMu0_TauMu == 1 && fMT2tree->muo[i].lv.Pt() > 27.0 && chargeMu == 0){
	  LeadingMuon = fMT2tree->muo[i].lv;
	  chargeMu = fMT2tree->muo[i].Charge;
	}
	else
	  if(fMT2tree->muo[i].RejMu_TauTau == 1){
	    chargeMu = 0;
	    break;
	  }
      }
      
      if(chargeMu == 0)
	continue;

      for(int i = 0; i < fMT2tree->NEles; i++){
	if(fMT2tree->ele[i].IDVetoMuTau){
	  chargeMu = 0;
	}
      }

      if(chargeMu == 0)
	continue;

      float MaxPtMuTau = -1;

      int TauInd = -1;

      for(int t=0; t<fMT2tree->NTaus; t++){ 
	  

	if(fMT2tree->tau[t].PassQCDTau_MuTau == 0)
	  continue;

	float deltaR = Util::GetDeltaR(fMT2tree->tau[t].lv.Eta(), LeadingMuon.Eta(), fMT2tree->tau[t].lv.Phi(), LeadingMuon.Phi());
	
	if(deltaR < 0.2)
	  continue;

	float maxPtMuTau = fMT2tree->tau[t].lv.Pt() + LeadingMuon.Pt();
	if(maxPtMuTau > MaxPtMuTau){
	  MaxPtMuTau = maxPtMuTau;
	  TauInd = t;
	}
      }

      if(TauInd == -1)
	continue;

      if((chargeMu + fMT2tree->tau[TauInd].Charge) != 0)
	continue;
      
      float myMT2 = fMT2tree->CalcMT2(0, false, fMT2tree->tau[TauInd].lv, LeadingMuon, fMT2tree->pfmet[0]);

      if(fMT2tree->tau[TauInd].PassTau_MuTau == 0)
 	continue;
      
      //myQuantity = fMT2tree->pileUp.NVertices;//nExtraTaus;

      //following to check the promptness and faking of the events:
 //      TString myQuantity = fMT2tree->GenLeptonAnalysisInterpretation(0,1,1, false);
      float myQuantity = 100.0;
      
      for(int j=0; j<fMT2tree->NJets; ++j){ 
	if (!((fMT2tree->jet[j].lv.Pt() > 20)  &&  fabs(fMT2tree->jet[j].lv.Eta()<2.3)))  
	  continue;

	float deltaR = Util::GetDeltaR(LeadingMuon.Eta(), fMT2tree->jet[j].lv.Eta(), LeadingMuon.Phi(), fMT2tree->jet[j].lv.Phi());
	if(deltaR < myQuantity)
	  myQuantity = deltaR;
	
      }
	
    

      //if(fMT2tree->tau[TauInd].PassTau_MuTau != 1 || fMT2tree->tau[TauInd].Isolation3Hits != 1.)
      //  continue;

     /*
	std::vector<int> Tau0;
	std::vector<int> Mu0;

	for(int i=0; i<fMT2tree->NTaus; ++i){
	if(fMT2tree->tau[i].PassTau_MuTau)
	Tau0.push_back(i);
	}
	for(int i=0; i<fMT2tree->NMuons; ++i){
	if(fMT2tree->muo[i].PassQCDMu0_MuMu)
	Mu0.push_back(i);
	}
	std::pair<int,int> indecies = fMT2tree->MuTauParing(Tau0,Mu0);

	int selected = 0;

	if(indecies.first != -1 && indecies.second != -1){
	float pairCharge = fMT2tree->tau[indecies.first].Charge + fMT2tree->muo[indecies.second].Charge;
	float Mass = (fMT2tree->tau[indecies.first].lv + fMT2tree->muo[indecies.second].lv).M();
	if(pairCharge != 0 && (Mass > 15.0 && (Mass < 45.0 || Mass > 75))){
	myQuantity = fMT2tree->CalcMT2(0, false, fMT2tree->tau[indecies.first].lv, fMT2tree->muo[indecies.second].lv, fMT2tree->pfmet[0]);
	tauIndex = indecies.first;
	selected = 1;
	}
	}

  //    int jetCounter = 0;
      
//       for(int j=0; j<fMT2tree->NJets; ++j){ 
// 	if(fMT2tree->jet[j].isPFIDLoose==false) continue;
// 	if (!((fMT2tree->jet[j].lv.Pt() > 20)  &&  fabs(fMT2tree->jet[j].lv.Eta()<2.3)))  
// 	  continue;
	
// 	jetCounter++;

// 	if (jetCounter==1)
// 	  continue;

// 	if(fMT2tree->jet[j].isTauMatch < 0)
// 	  continue;
// 	myQuantity = fMT2tree->tau[fMT2tree->jet[j].isTauMatch].lv.Pt();


	if(selected == 0)
	continue;
      */
      if(data == 1){
      
	MT2[NumberOfSamples]->Fill(myQuantity, weight);//data
      
      }else{
	if(Sample.sname == "SUSY")
	  MT2[NumberOfSamples-1]->Fill(myQuantity, weight);
	else
	  MT2[NumberOfSamples-2]->Fill(myQuantity, weight);
      
	if(Sample.sname == "Top")
	  MT2[3]->Fill(myQuantity, weight);
	else
	  if(Sample.sname == "DY")	
	    MT2[2]->Fill(myQuantity, weight);
	  else
	    if(Sample.sname == "Wtolnu"){
	      MT2[1]->Fill(myQuantity, weight);}
	    else
	      if(Sample.sname == "QCD")
		MT2[0]->Fill(myQuantity, weight);
	      else
		if(Sample.sname == "VV")
		  MT2[4]->Fill(myQuantity, weight);
		else
		  if(Sample.sname == "Higgs")
		    MT2[5]->Fill(myQuantity, weight);
      }

}
  
  }//for(ii<fSamples)


  for(int j = 0; j <= (NumberOfSamples); j++){
    AddOverAndUnderFlow(MT2[j], true, true);
  }
   printYield();

  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j <= (NumberOfSamples); j++){
    // MT2[j]->Rebin(3);
    if(j <= (NumberOfSamples - 3))
      h_stack -> Add(MT2[j]);
  }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W", "f");
  Legend1->AddEntry(MT2[2], "ZX", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW", "f");
  Legend1->AddEntry(MT2[5], "Higgs", "f"); 
  Legend1->AddEntry(MT2[NumberOfSamples-1], "SMS", "l");
  Legend1->AddEntry(MT2[NumberOfSamples], "data", "l");


  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + myfileName +"_Histos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  h_stack->Write();

  for(int j = 0; j <= (NumberOfSamples); j++)
    MT2[j]->Write();
  Legend1->Write();

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  printHisto(h_stack, MT2[NumberOfSamples], MT2[NumberOfSamples-2], MT2[NumberOfSamples-1], Legend1 , "MTC", "hist", true, "MT2", "Events", 0, -10, 2, true);

  plotRatioStack(h_stack, MT2[NumberOfSamples-2], MT2[NumberOfSamples], MT2[NumberOfSamples-1], true, false, "MT2_ratio", Legend1, "MT2", "Events", 0, -10, 2, true, "C");


  TCanvas *MyC = new TCanvas("MyC", "MyC");
  MyC->Divide(2,1);
  MyC->cd(1);
  TGraphAsymmErrors* sig1 = plotSig(MT2[NumberOfSamples-1], MT2[NumberOfSamples-2], "MT2", "Lower Cut", /*type = */ 2, /*sys = */ 0.1);
  sig1->Draw("ACP");

  MyC->cd(2);
  TGraphAsymmErrors* sig2 = plotSig(MT2[NumberOfSamples-1], MT2[NumberOfSamples-2], "MT2", "Upper Cut", /*type = */ 2, /*sys = */ 0.1);
  sig2->Draw("ACP");
}

void MassPlotter::elemuAnalysis(TString cuts, TString trigger, unsigned int nevents, TString myfileName ,int type){

  TH1::SetDefaultSumw2();

  TString cnames[NumberOfSamples+1] = {"QCD", "Wjets", "Zjets", "Top","WWjets", "MC", "susy","data"};
  int ccolor[NumberOfSamples+1] = { 401, 417, 419, 600, 603,500, 1, 632};
  TString varname = "MT2";
  
    for (int i=0; i<(NumberOfSamples+1); i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i], "", 1000, 0, 1000);
  //MT2[i] = new TH1D(varname+"_"+cnames[i], "", 11, 0, 11);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }
  //  TString  genStatus[11] = {"pp", "pf", "fp", "ff", "tp", "pt","tf","ft","tt", "nothing", "wrong"};
  //  for(int k = 0; k < MT2[i]->GetNbinsX(); k++)
   //   MT2[i] ->GetXaxis()->SetBinLabel(k+1,genStatus[k]);




  MT2[7] -> SetMarkerStyle(20);
  MT2[7] -> SetMarkerColor(kBlack);
  MT2[7] -> SetLineColor(kBlack);
  
  MT2[5] -> SetFillStyle(3004);
  MT2[5] -> SetFillColor(kBlack);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data"){
      data = 1;
      myCuts += " && " + trigger;
    }

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over : " <<endl;	
    cout << " Name: " << Sample.name << endl;
    cout << " File: " << (Sample.file)->GetName() << endl;
    cout << " Events: " << Sample.nevents << endl;
    cout << " Events in tree: " << Sample.tree->GetEntries() << endl;
    cout << " Xsection: " << Sample.xsection << endl;
    cout << " kfactor: " << Sample.kfact << endl;
    cout << " avg PU weight: " << Sample.PU_avg_weight << endl;
    cout << " Weight: " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry);
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){
fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
fflush(stdout);
      }
 
      float weight = Weight;

      if(data == 1)
       weight = 1.0;
      else{

      float muIdSFeleMu = fMT2tree->muo[fMT2tree->eleMu[0].mu0Ind].GetMuIDISOSFelemu();
	

	
	float muTrgSFeleMu = fMT2tree->muo[fMT2tree->eleMu[0].mu0Ind].GetMuTrgSFelemu();
	

      
	float eleIdSFeleMu = fMT2tree->ele[fMT2tree->eleMu[0].ele0Ind].GetEleIDISOSFelemu();
	


       float eleTrgSFeleMu = fMT2tree->ele[fMT2tree->eleMu[0].ele0Ind].GetEleTrgSFelemu();
  

       
	weight = Weight * muIdSFeleMu *  muTrgSFeleMu * eleIdSFeleMu * eleTrgSFeleMu;   

if(Sample.type != "susy")
weight *= (fMT2tree->pileUp.Weight * fMT2tree->SFWeight.BTagCSV40eq0/Sample.PU_avg_weight);// * fMT2tree->SFWeight.TauTagge1/Sample.PU_avg_weight);//
      }
      
      float myQuantity = fMT2tree->eleMu[0].MT2;
    //TString myQuantity = fMT2tree->GenLeptonAnalysisInterpretation( 1, 1, 0, false );

if(data == 1){
      
      MT2[7]->Fill(myQuantity, weight);//data
      
    }else{
      if(Sample.sname == "SUSY")
MT2[6]->Fill(myQuantity, weight);
      else
MT2[5]->Fill(myQuantity, weight);
  
if(Sample.sname == "VV")
MT2[4]->Fill(myQuantity, weight);
    
    else if(Sample.sname == "Top")
MT2[3]->Fill(myQuantity, weight);
      else
if(Sample.sname == "DY")	
MT2[2]->Fill(myQuantity, weight);
else
if(Sample.sname == "Wtolnu")
MT2[1]->Fill(myQuantity, weight);
else
if(Sample.sname == "QCD")
MT2[0]->Fill(myQuantity, weight);
    }}
  
  }



  for(int j = 0; j < (NumberOfSamples); j++){
    AddOverAndUnderFlow(MT2[j], true, true);
  }
  printYield();

  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j < (NumberOfSamples); j++){
    // MT2[j]->Rebin(3);
    TH1F* mt2 = (TH1F*)MT2[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 3))
      h_stack -> Add(MT2[j]);
    delete mt2;
  }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WWjets", "f");
  Legend1->AddEntry(MT2[6], "SMS", "l");
   Legend1->AddEntry(MT2[7], "data", "l");
  // vector<TH1D*> h_signals;
  // h_signals.push_back(&(*MT2[5]));
  // TLegend *Legend1;

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + myfileName +"_Histos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  h_stack->Write();
  MT2[0]->Write();
  MT2[1]->Write();
  MT2[2]->Write();
  MT2[3]->Write();
  MT2[4]->Write();
  MT2[5]->Write();
  MT2[6]->Write();
  MT2[7]->Write();
  Legend1->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  printHisto(h_stack, MT2[7], MT2[5], MT2[6], Legend1 , "MTC", "hist", true, "MT2", "Events", 0, -10, 2, true);

  plotRatioStack(h_stack, MT2[5], MT2[7], MT2[6], true, false, "MT2_ratio", Legend1, "MT2", "Events", 0, -10, 2, true);

  TCanvas *MyC = new TCanvas("MyC", "MyC");
  MyC->Divide(2,1);
  MyC->cd(1);
  TGraphAsymmErrors* sig1 = plotSig(MT2[6], MT2[5], "MT2", "Lower Cut", /*type = */ 2, /*sys = */ 0.1);
  sig1->Draw("ACP");

  MyC->cd(2);
  TGraphAsymmErrors* sig2 = plotSig(MT2[6], MT2[5], "MT2", "Upper Cut", /*type = */ 2, /*sys = */ 0.1);
  sig2->Draw("ACP");


}

 void MassPlotter::EleFakeRateforeleMu(TString cuts, TString trigger, unsigned int nevents, TString myfileName){
 TH1::SetDefaultSumw2();

 
 TH2F *hPtEtaAll  = new TH2F("hPtEtaAll", "hPtEtaAll",  60, -3.0, 3.0, 1000, 0, 1000);
 TH2F *hPtEtaPass = new TH2F("hPtEtaPass","hPtEtaPass", 60, -3.0, 3.0, 1000, 0, 1000);
 TH2F *hPtEtaFakeRate = new TH2F("hPtEtaFakeRate","hPtEtaFakeRate", 60, -3.0, 3.0, 1000, 0, 1000);
  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
   TString myCuts = cuts;
 
   int data = 0;
   sample Sample = fSamples[ii];
    
   if(Sample.type == "data"){
     data = 1;
     myCuts += " && " + trigger;
   }else 
     continue;
   

   fMT2tree = new MT2tree();
   Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

   float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

   std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "   looping over :  " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry); 
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

     float weight = Weight;

      if(data == 1)
 	weight = 1.0;
      else{
	if(Sample.type != "susy")
	  weight *= (fMT2tree->pileUp.Weight * fMT2tree->SFWeight.BTagCSV40eq0/Sample.PU_avg_weight);
      }      

      int hasMuon = 0;
      TLorentzVector LeadingMuon;
      
        for(int i=0; i<fMT2tree->NMuons; ++i){ 
	if(fMT2tree->muo[i].PassMu0_EleMu== 1){
	  hasMuon = 1;
	  LeadingMuon = fMT2tree->muo[i].lv;
	  break;
	}
      }
      
        if(hasMuon == 0)
	continue;

      float MaxPtEleMu = -1;

      int EleInd = -1;

      for(int t=0; t<fMT2tree->NEles; t++){ 
	  float deltaR = Util::GetDeltaR(fMT2tree->ele[t].lv.Eta(), LeadingMuon.Eta(), fMT2tree->ele[t].lv.Phi(), LeadingMuon.Phi());
	
	if(deltaR < 0.1)
	  continue;

	float maxPtEleMu= fMT2tree->ele[t].lv.Pt() + LeadingMuon.Pt();

	if(maxPtEleMu > MaxPtEleMu){
	  MaxPtEleMu = maxPtEleMu;
	  EleInd = t;
	}
      }

      if(EleInd == -1)
	continue;

      

      
	hPtEtaAll->Fill(fMT2tree->ele[EleInd].lv.Eta(), fMT2tree->ele[EleInd].lv.Pt()); 
	
     
      
      if(fMT2tree->ele[EleInd].IDSelEMU== 1){
	hPtEtaPass->Fill(fMT2tree->ele[EleInd].lv.Eta(), fMT2tree->ele[EleInd].lv.Pt()); 
	
      }
    }
 }
 TString fileName = fOutputDir;
 if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_FRHistos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  hPtEtaAll->Write();
  hPtEtaPass->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

 TCanvas *MyC = new TCanvas("Fake","Fake");
 MyC->Divide(2,2);
 MyC->cd(1);
 hPtEtaPass->Draw();
 MyC->cd(2);
 hPtEtaAll->Draw();
 MyC->cd(3);
 hPtEtaFakeRate->Divide(hPtEtaPass,hPtEtaAll);

 }

void MassPlotter::EleEfficiencyforelemu(TString cuts, unsigned int nevents, TString myfileName){
  TH1::SetDefaultSumw2();
 
//   TH2F *hPtEtaAll  = new TH2F("hPtEtaAll", "hPtEtaAll",  60, -3.0, 3.0, 1000, 0, 1000);
//   TH2F *hPtEtaPass = new TH2F("hPtEtaPass","hPtEtaPass", 60, -3.0, 3.0, 1000, 0, 1000);
     TH1F *hPtEtaAll  = new TH1F("hPtEtaAll", "hPtEtaAll",  200, 0, 200);
     TH1F *hPtEtaPass = new TH1F("hPtEtaPass","hPtEtaPass", 200, 0, 200);
     TH1F *hPtEtaPromptRate = new TH1F("hPtEtaPromptRate","hPtEtaPromptRate", 200, 0, 200);


     for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
     TString myCuts = cuts;
 

     sample Sample = fSamples[ii];
    
     if(Sample.sname != "DY")
     continue;

     fMT2tree = new MT2tree();
     Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

     float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry); 
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

     float weight = Weight;
     
     weight *= (fMT2tree->pileUp.Weight * fMT2tree->SFWeight.BTagCSV40eq0/Sample.PU_avg_weight);//* fMT2tree->SFWeight.TauTagge1/Sample.PU_avg_weight);//
 
     int GenEle = 0;
     
     for(int j = 0; j < fMT2tree->NGenLepts; j++){
     if(abs(fMT2tree->genlept[j].ID) != 11)
     continue;
     GenEle++;

	
	
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();

	float genleptPhi = fMT2tree->genlept[j].lv.Phi();
	
	
	for(int t = 0;t < fMT2tree->NEles; t++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->ele[t].lv.Eta(), genleptPhi, fMT2tree->ele[t].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    
	  }
	

	if(minDR < 0.1){
	  
	    hPtEtaAll->Fill(fMT2tree->eleMu[0].MT2, weight);
	     
	 

	  if(fMT2tree->ele[t].IDSelEMU== 1){
	    hPtEtaPass->Fill(fMT2tree->eleMu[0].MT2, weight);
	        
	  }
	}
     }
    }
  }
 }
  
     TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_PRHistos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  hPtEtaAll->Write();
  hPtEtaPass->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" cuts "<<cuts<<endl;
  TCanvas *MyC = new TCanvas("MyC", "MyC");
  MyC->Divide(2,2);
  MyC->cd(1);
  hPtEtaPass->Draw();
  MyC->cd(2);
  hPtEtaAll->Draw();
  MyC->cd(3);
  hPtEtaPromptRate->Divide(hPtEtaPass,hPtEtaAll);
 }

float FRWeight(float fakeRate, float promptRate, int isolated){

  float fRweight = fakeRate * (1 - promptRate);
  
  if(isolated != 1)
    fRweight -= fakeRate;

  fRweight /= (fakeRate - promptRate);
  
  return fRweight;
}


void MassPlotter::muTauWJetsEstimation(TString cuts, TString trigger, TString myfileName, bool calculateTauMTPass, float sysFR, float sysPR, float sysMT2W, float sysTauMT){

  TH1::SetDefaultSumw2();
  setFlags(10);
  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  fileName = fileName + myfileName;
  TFile *file = new TFile(fileName.Data(), "READ");

  TH2* hPtEtaAll  = (TH2*) file->Get("hPtEtaAll");
  TH2* hPtEtaPass = (TH2*) file->Get("hPtEtaPass");

  float xbin1[4] = {0.0, 50.0, 100.0, 1000.0}; //Barrel
  TH1F* hPt1All  = new TH1F("hPt1All",  "Barrel",  3,xbin1);
  TH1F* hPt1Pass = new TH1F("hPt1Pass", "Barrel", 3,xbin1);

  float xbin2[4] = {0.0, 50.0, 100.0, 1000.0}; //Endcap
  TH1F* hPt2All  = new TH1F("hPt2All",  "Endcap",  3,xbin2);
  TH1F* hPt2Pass = new TH1F("hPt2Pass", "Endcap", 3,xbin2);
//   float xbin1[2] = {0.0, 1000.0}; //Barrel
//   TH1F* hPt1All  = new TH1F("hPt1All",  "Barrel",  1,xbin1);
//   TH1F* hPt1Pass = new TH1F("hPt1Pass", "Barrel", 1,xbin1);

//   float xbin2[2] = {0.0, 1000.0}; //Endcap
//   TH1F* hPt2All  = new TH1F("hPt2All",  "Endcap",  1,xbin2);
//   TH1F* hPt2Pass = new TH1F("hPt2Pass", "Endcap", 1,xbin2);
 

  for(int i = 1; i < 61; i++){
    for(int j = 1; j < 1001; j++){

      double binIJAll  = hPtEtaAll ->GetBinContent(i,j);
      double binIJPass = hPtEtaPass->GetBinContent(i,j);
      
      if(i < 17 || i > 44){//Endcap
	hPt2All ->Fill(j - 0.5, binIJAll);
	hPt2Pass->Fill(j - 0.5, binIJPass);
      }else{
	hPt1All ->Fill(j - 0.5, binIJAll);
	hPt1Pass->Fill(j - 0.5, binIJPass);
      }
    }
  }

  hPt1Pass->Divide(hPt1All);
  
  hPt2Pass->Divide(hPt2All);

  static const int nbins = 6;//3;//11
  //  double bins[nbins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,125.0,150.0,175.0,200.0,250.0,300.0,400.0};      //MT2
  //  double bins[nbins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0,120.0,200.0};      //MT2
    double bins[nbins+1] = {0.0,20.0,40.0,50.0,70.0,90.0,300.0};      //MT2

  //double bins[nbins+1] = {40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0};      //MT2

  //double bins[nbins+1] = {0.0,18.0,36.0,54.0,72.0,90.0,300.0};      //MT2
  //double bins[nbins+1] = {0.0, 50.0, 100.0, 1000.0};
  //    double binsTauPt[nbins+1] = {-20.0,40.0,60.0,80.0,100.0,150,1300.0};      //tauPt
  /*
  fileName = fOutputDir + "/MT2_MuTau_Over_QCDMuTau_SignalSelectionNoZVeto_PRHistos.root";
  
  TFile *file2 = new TFile(fileName.Data(), "READ");

  TH1* hOldMT2All  = (TH1*) file2->Get("hPtEtaAll");
  TH1* hOldMT2Pass = (TH1*) file2->Get("hPtEtaPass");

//   TH1* hMT2All   = hOldMT2All ->Rebin(nbins, "hMT2All" , bins);
//   TH1* hMT2Pass  = hOldMT2Pass->Rebin(nbins, "hMT2Pass", bins);
  TH1* hMT2All   = hOldMT2All->Rebin(10, "hMT2All");
  TH1* hMT2Pass  = hOldMT2Pass->Rebin(10, "hMT2Pass");
  
//   for(int j = 1; j < hMT2All->GetNbinsX(); j++){
//     cout<<" bin "<<j<<" hMT2All  "<<hMT2All ->GetBinContent(j)<<" +- "<<hMT2All ->GetBinError(j)<<endl;
//     cout<<" bin "<<j<<" hMT2Pass "<<hMT2Pass->GetBinContent(j)<<" +- "<<hMT2Pass->GetBinError(j)<<endl;
//   }

//   TH1F* hMT2All  = new TH1F("hMT2All",  "All",  nbins, bins);
//   TH1F* hMT2Pass = new TH1F("hMT2Pass", "Pass", nbins, bins);

//   for(int j = 1; j < 201; j++){
    
//     double binJAll  = hOldMT2All ->GetBinContent(j);
//     double binJPass = hOldMT2Pass->GetBinContent(j);
    
//     hMT2All ->Fill(j - 0.5, binJAll);
//     hMT2Pass->Fill(j - 0.5, binJPass);
             
//   }
  
  hMT2Pass->Divide(hMT2All);
//   for(int j = 1; j < hMT2All->GetNbinsX(); j++){
//     cout<<" bin "<<j<<" hMT2Pass "<<hMT2Pass->GetBinContent(j)<<" +- "<<hMT2Pass->GetBinError(j)<<endl;
//   }
*/

  TH1F *myWeights       = new TH1F("myWeights",       "myWeights",       30, -0.1, 0.2);
  TH2F *WeightsFakeRate = new TH2F("WeightsFakeRate", "WeightsFakeRate", 30, 0, 0.025, 30, -0.1, 0.2);
 
  TString  cnames[NumberOfSamples+1] = {"QCD",     "W",    "ZX",  "Top",       "WW",  "Higgs",   "MC", "susy","data"};
  int      ccolor[NumberOfSamples+1] = { 401,      417,     419,   855,         603,     kRed,    603,      1, 632};

  TH1D *LooseNonTight = new TH1D("LooseNonTight", "LooseNonTight", nbins, bins);
  TH1D *Tight = new TH1D("Tight", "Tight", nbins, bins);


  TString varname = "MT2";
  for (int i=0; i<=(NumberOfSamples); i++){
    if(calculateTauMTPass || i == NumberOfSamples){
      MT2[i] = new TH1D(varname+"_"+cnames[i], "", nbins, bins);
      MT2[i] -> SetFillColor (ccolor[i]);
      MT2[i] -> SetLineColor (ccolor[i]);
      MT2[i] -> SetLineWidth (2);
      MT2[i] -> SetMarkerColor(ccolor[i]);
      MT2[i] -> SetStats(false);
    }
  }

  MT2[NumberOfSamples] -> SetMarkerStyle(20);
  MT2[NumberOfSamples] -> SetMarkerColor(kBlack);
  MT2[NumberOfSamples] -> SetLineColor(kBlack);

  if(calculateTauMTPass){
    MT2[NumberOfSamples-2] -> SetFillStyle(3004);
    MT2[NumberOfSamples-2] -> SetFillColor(kBlack);
  }

  TH1D *EstimPRUp   =  new TH1D("EstimPRUp",   "EstimPRUp",   nbins, bins);
  TH1D *EstimPRDown =  new TH1D("EstimPRDown", "EstimPRDown", nbins, bins);
  TH1D *EstimFRUp   =  new TH1D("EstimFRUp",   "EstimFRUp",   nbins, bins);
  TH1D *EstimFRDown =  new TH1D("EstimFRDown", "EstimFRDown", nbins, bins);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  if(calculateTauMTPass)
    tauMTPass = new TH1F("tauMTPass", "tauMTPass", nbins, bins);

  TH1F *tauMTAll  = new TH1F("tauMTAll",  "tauMTAll",  nbins, bins);

  TH1F *MT2Weight = new TH1F("MT2Weight", "MT2Weight", nbins, bins);
 
  TH1F *Sigma_wi2 = new TH1F("Sigma_wi2", "Sigma_wi2", nbins, bins);

  TH1F *MT2EstSys = new TH1F("MT2EstSys", "MT2EstSys", nbins, bins);
 
  //QCDTau_MuTau -> TightTau_MuTau  
//   MT2Weight->SetBinContent(1,1.332621);
//   MT2Weight->SetBinContent(2,1.205256);
//   MT2Weight->SetBinContent(3,1.040942);
//   MT2Weight->SetBinContent(4,0.93485);
//   MT2Weight->SetBinContent(5,0.7379218);
//   MT2Weight->SetBinContent(6,1.576364);
//   MT2Weight->SetBinError(1,0.3319932);
//   MT2Weight->SetBinError(2,0.4059702);
//   MT2Weight->SetBinError(3,0.3366078);
//   MT2Weight->SetBinError(4,0.298416);
//   MT2Weight->SetBinError(5,0.2766776);
//   MT2Weight->SetBinError(6,1.403669);

  //LooseTau_MuTau -> TightTau_MuTau 
  MT2Weight->SetBinContent(1,0.937382);
  MT2Weight->SetBinContent(2,0.9769556);
  MT2Weight->SetBinContent(3,0.9696749);
  MT2Weight->SetBinContent(4,0.9930092);
  MT2Weight->SetBinContent(5,1.048872);
  MT2Weight->SetBinContent(6,1.609784);
  MT2Weight->SetBinError(1,0.234823);
  MT2Weight->SetBinError(2,0.329825);
  MT2Weight->SetBinError(3,0.3149124);
  MT2Weight->SetBinError(4,0.3181093);
  MT2Weight->SetBinError(5,0.4041598);
  MT2Weight->SetBinError(6,1.581523);


    for(int j = 0; j < MT2Weight->GetNbinsX(); j++){
      float value = MT2Weight->GetBinContent(j+1);

      float error = MT2Weight->GetBinError(j+1);
      
      float newError = sqrt(error * error + sysMT2W * sysMT2W * value * value);

      MT2Weight->SetBinError(j+1, newError);
    }

 
  // vector of all histos
  vector<TH1D*> h_samples;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){

    TString myCuts = cuts;
    
    sample Sample = fSamples[ii];

    h_samples.push_back(new TH1D(varname+"_"+Sample.name, "", nbins, bins));
 
    h_samples[ii] -> Sumw2();
    h_samples[ii] -> SetLineColor(Sample.color);
   
    h_samples[ii] -> SetMarkerColor(Sample.color);
    h_samples[ii] -> SetStats(false);
    if(Sample.type == "susy" ){
      h_samples[ii] -> SetLineColor(kBlack);
      h_samples[ii] -> SetLineStyle(kDotted);
    }

//     if(!calculateTauMTPass && Sample.type != "data")
//       continue;

//      if(Sample.sname != "Wtolnu")//Closure
//        continue;

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    float Weight;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    if(Sample.type == "data"){
   // if(Sample.type == "mc"){//To run on MC

      //myCuts += "&& (muTau[0].Isolated == 1) && (tau[muTau[0].tau0Ind].Isolation3Hits == 1.)";
      //myCuts += "&& (muTau[0].Isolated >= 0)";//QCD->Tight
      myCuts += "&& (muTau[0].Isolated > 0)";//Loose->Tight
      if(Sample.type == "data" && Sample.sname != "Wtolnu")    myCuts += "&&" + trigger;

      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
      
      unsigned int nentries = myEvtList->GetN();
      
      for (unsigned int jentry=0; jentry<nentries; jentry++) {
	
	Sample.tree->GetEntry(myEvtList->GetEntry(jentry));
	
	if ( fVerbose>2 && jentry % 100000 == 0 ){
	  fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	  fflush(stdout);
	}
	
	//int muIndex = fMT2tree->muTau[0].GetMuIndex0();

	//if((fMT2tree->muTau[0].GetIsolated() != 1) || (fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].Isolation3Hits > 1.))
	//continue;

	if(Sample.type == "data" && Sample.sname == "Wtolnu"){
	  int genMuTauStatus = fMT2tree->GenLeptonAnalysis(0,1,1);//Closue
	  
	  if(genMuTauStatus != 613 && genMuTauStatus != 623)// && genMuTauStatus != 633 && genMuTauStatus != -1)//Closure
	    continue;       
	}
	
	if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){
	  
	  Weight = Sample.lumi;
	  if(fPUReweight) Weight /= Sample.PU_avg_weight;
	  
	}

	float weight = Weight;
	
	weight *= fMT2tree->SFWeight.BTagCSV40eq0 * fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

	if(fPUReweight)
	  weight *= fMT2tree->pileUp.Weight;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	
	if(Sample.type == "data" && Sample.sname != "Wtolnu")
	  weight = 1;

	float myQuantity = fMT2tree->muTau[0].GetMT2();//fMT2tree->muo[muIndex].MT;

	int tauIndex = fMT2tree->muTau[0].GetTauIndex0();

	float tauEta = fMT2tree->tau[tauIndex].lv.Eta();

	float tauPt  = fMT2tree->tau[tauIndex].lv.Pt();
	
// 	myQuantity = tauPt;
	float fakeRate    = 0.0;
 	float fakeRateErr = 0.0;

	if(fabs(tauEta) < 1.4){
	  int binNumber = hPt1Pass->FindBin(tauPt);
	  fakeRate      = hPt1Pass->GetBinContent(binNumber);
 	  fakeRateErr   = hPt1Pass->GetBinError(binNumber);
	}else{
	  int binNumber = hPt2Pass->FindBin(tauPt);
	  fakeRate      = hPt2Pass->GetBinContent(binNumber);
 	  fakeRateErr   = hPt2Pass->GetBinError(binNumber);
	}

	TH1F fakeRateOS_SSCorrection("", "", 1, 0, 1);
	TH1F fakeRateUnCorrected("", "", 1, 0, 1);

	//PtEta_MuTauTight_Over_Loose_pfOnly_WJets_OppSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root  0.489688 +/- 0.003775
	//PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root 0.428217 +/- 0.00647227
	fakeRateOS_SSCorrection.SetBinContent(1, 1.14355);
	fakeRateOS_SSCorrection.SetBinError(1, 0.0194);

	//PtEta_MuTauTight_Over_Loose_pfOnly_WJets_OppSign_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root 0.467944 +/- 0.00564925
	//PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root 0.41091 +/- 0.0090992
// 	fakeRateOS_SSCorrection.SetBinContent(1, 1.138799);
// 	fakeRateOS_SSCorrection.SetBinError(1, 0.02872);


// PtEta_MuTauTight_Over_Loose_SingleMu_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root
// 0.4464 +/- 0.0062
  	fakeRate      = 0.4464;
	fakeRateErr   = 0.0062;


// PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_ExtraLepVeto_MET_NBJets_Weighted_FRHistos.root
// 0.4282 +/- 0.0065
// 	fakeRate      = 0.4282;
// 	fakeRateErr   = 0.0065;
// PtEta_MuTauTight_Over_Loose_pfOnly_WJets_SameSign_ExtraLepVeto_MT2gt40_ZVeto_minDPhi_MET_NBJets_Weighted_FRHistos.root
// 0.41091 +/- 0.0090992
//	fakeRate      = 0.41091;  
//	fakeRateErr   = 0.0090992;

	fakeRateUnCorrected.SetBinContent(1, fakeRate);
	fakeRateUnCorrected.SetBinError(1, fakeRateErr);

 	fakeRateUnCorrected.Multiply(&fakeRateOS_SSCorrection);

	fakeRate = fakeRateUnCorrected.GetBinContent(1);
	fakeRateErr = fakeRateUnCorrected.GetBinError(1);

	float SysFR = sqrt(fakeRateErr * fakeRateErr + sysFR * sysFR * fakeRate * fakeRate);

	// P +  F = Loose
	//pP + fF = Loose - LooseNonTight
	//FakeContribution = f * F 
	//F * (f - p) = (1 - p) Loose - LooseNonTight
	//F * (f - p) = (1 - p) Tight - p * LooseNonTight

	float promptRate = 0.52;// 0.5349 +- 0.0031 QCD To Tight MT2_MuTau_Over_QCDMuTau_SignalSelectionNoZVeto_DYToLL_M50_TauTightIso_DY_PRHistos.root
	promptRate = 0.77; //0.7659 +- 0.0032  Loose to Tight MT2_MuTauTight_Over_Loose_SignalSelectionNoZVeto_DYToLL_M50_DY_PRHistos.root

	float promptRateErr = 0.0032;

	float SysPR = sqrt(promptRateErr * promptRateErr + sysPR * sysPR * promptRate * promptRate);

	int isolated = 0;//fMT2tree->muTau[0].GetIsolated();

	if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].Isolation3Hits == 1.)
	  isolated = 1;

	float fRweight = FRWeight(fakeRate, promptRate, isolated);

// 	myWeights->Fill(fRweight);
	
	WeightsFakeRate->Fill(fakeRate, fRweight);

 	MT2[NumberOfSamples]->Fill(myQuantity, fRweight * weight);
	
	if(isolated == 1)
	  Tight->Fill(myQuantity, weight);
	else
	  LooseNonTight->Fill(myQuantity, weight);
	
	float fRweightFRUp = FRWeight((fakeRate + SysFR), promptRate, isolated);

	myWeights->Fill(fRweightFRUp);

	EstimFRUp->Fill(myQuantity, fRweightFRUp * weight);
 
// 	fRweightFRUp -= fRweight;

// 	float fRweightFRUp2 = fRweightFRUp * fRweightFRUp;

 	float fRweightFRDown = FRWeight((fakeRate - SysFR), promptRate, isolated);

// 	fRweightFRDown -= fRweight;

// 	float fRweightFRDown2 = fRweightFRDown * fRweightFRDown;
	
// 	EstimFRDown->Fill(myQuantity, 0.5 * (fRweightFRUp2 + fRweightFRDown2) * weight);
 	EstimFRDown->Fill(myQuantity, fRweightFRDown * weight);
 
 	float fRweightPRUp = FRWeight(fakeRate, (promptRate + SysPR), isolated);

	EstimPRUp->Fill(myQuantity, fRweightPRUp * weight);
 
// 	fRweightPRUp -= fRweight;

// 	float fRweightPRUp2 = fRweightPRUp * fRweightPRUp;

 	float fRweightPRDown = FRWeight(fakeRate, (promptRate - SysPR), isolated);

// 	fRweightPRDown -= fRweight;

// 	float fRweightPRDown2 = fRweightPRDown * fRweightPRDown;
	
//	EstimPRDown->Fill(myQuantity, 0.5 * (fRweightPRUp2 + fRweightPRDown2) * weight);
	EstimPRDown->Fill(myQuantity, fRweightPRDown * weight);
      }
    }else{
// 	    continue;//To run on MC
      if(Sample.tree->GetEntries() == 0)
	continue;
      
      myCuts += "&& (muTau[0].Isolated == 1) && (tau[muTau[0].tau0Ind].Isolation3Hits == 1.)";
  
      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

      unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();

      for (unsigned int jentry=0; jentry<nentries;jentry++) {
    
	Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

	if ( fVerbose>2 && jentry % 100000 == 0 ){
	  fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	  fflush(stdout);
	}

	int genMuTauStatus = fMT2tree->GenLeptonAnalysis(0,1,1);

	if(genMuTauStatus != 613 && genMuTauStatus != 623 && genMuTauStatus != 633 && genMuTauStatus != -1)
	  continue;       

	if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){
	  
	  Weight = Sample.lumi;
	  if(fPUReweight) Weight /= Sample.PU_avg_weight;
	  
	}
	
	float weight = Weight;
	
	weight *= fMT2tree->SFWeight.BTagCSV40eq0 * fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

	if(fPUReweight)
	  weight *= fMT2tree->pileUp.Weight;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	float myQuantity = fMT2tree->muTau[0].GetMT2();

	int tauIndex = fMT2tree->muTau[0].GetTauIndex0();

	float tauEta = fMT2tree->tau[tauIndex].lv.Eta();

	float tauPt  = fMT2tree->tau[tauIndex].lv.Pt();

// 	myQuantity = tauPt;

//  	if(fMT2tree->tau[tauIndex].MT > 200)
	  h_samples[ii]->Fill(myQuantity, weight);

	tauMTAll->Fill(myQuantity, weight);
	
  	if(fMT2tree->tau[tauIndex].MT > 200)
	  tauMTPass->Fill(myQuantity, weight);
	
      }

      AddOverAndUnderFlow(h_samples[ii]);

      if(Sample.type == "susy"){
	MT2[NumberOfSamples-1]->Add(h_samples[ii]);
      }else
	MT2[NumberOfSamples-2]->Add(h_samples[ii]);

      if (Sample.sname == "QCD") {
	MT2[0]->Add(h_samples[ii]);
      }
      else if(Sample.sname == "Top" ){
	MT2[3]->Add(h_samples[ii]);
      }
      else if(Sample.sname=="Wtolnu") {
	MT2[1]->Add(h_samples[ii]);
      }
      else if(Sample.sname=="DY") {
	MT2[2]->Add(h_samples[ii]);
      }
      else if(Sample.sname == "VV") {
	MT2[4]->Add(h_samples[ii]);
      }
      else if(Sample.sname == "Higgs") {
	MT2[5]->Add(h_samples[ii]);
      }
  

    }
  }
  //Error Caculations
  //S = Sigma(w_i)* MT2Weight * tauMTEff
  //(DeltaS/S)2  = (Delta_w/w)2 
  //(Delta_w/w)2 = Sigma(w_i2) + (Delta_w_fr)2 + (Delta_w_pr)2

  for(int j = 0; j < MT2[NumberOfSamples]->GetNbinsX(); j++){
    Sigma_wi2->SetBinContent(j+1, MT2[NumberOfSamples]->GetBinError(j+1));
  }
  
  int p = 0;
  Cout(p++, EstimPRDown);
  Cout(p++, EstimPRUp);
  Cout(p++, EstimFRDown);
  Cout(p++, EstimFRUp);
  /*
  EstimPRDown->Divide(MT2Weight);
  EstimPRUp->Divide(MT2Weight);

  EstimFRDown->Divide(MT2Weight);
  EstimFRUp->Divide(MT2Weight);

  Sigma_wi2->Divide(MT2Weight);
  */
  Cout(p++, EstimPRDown);
  Cout(p++, EstimPRUp);
  Cout(p++, EstimFRDown);
  Cout(p++, EstimFRUp);


  Cout(p++, MT2[NumberOfSamples]); 
  //MT2[NumberOfSamples]->Divide(MT2Weight);
  cout<<"estimation corrected by MT2 Weights "<<endl;

  Cout(p++, MT2[NumberOfSamples]); 
  Cout(p++, Tight);
  Cout(p++, LooseNonTight);

  if(calculateTauMTPass){

    for(int j = 0; j <= (NumberOfSamples); j++){
      AddOverAndUnderFlow(MT2[j], true, true);

  //     for(int i = 1; i <= MT2[j]->GetNbinsX(); i++){
// 	float val = MT2[j]->GetBinContent(i);
// 	float err = MT2[j]->GetBinError(i);
	
// 	MT2[j]->SetBinError(i, sqrt(err * err + 0.25 * 0.25 * val * val));
//       }
    }

    AddOverAndUnderFlow(tauMTAll,  true, true);
    
    Cout(p++, tauMTAll);
    
    AddOverAndUnderFlow(tauMTPass, true, true);
    
    Cout(p++, tauMTPass);
    
    tauMTPass->Divide(tauMTAll);
    
    tauMTPass->SetBinContent(0,0.00541489);
    tauMTPass->SetBinContent(1,0.002377618);
    tauMTPass->SetBinContent(2,0.002377618);
    tauMTPass->SetBinContent(3,0.002897278);
    tauMTPass->SetBinContent(4,0.008198583);
    tauMTPass->SetBinContent(5,0.04965259);
    tauMTPass->SetBinContent(6,0.06136292);
    tauMTPass->SetBinError(0,0.001346403);
    tauMTPass->SetBinError(1,0.0009614156);
    tauMTPass->SetBinError(2,0.0009614156);
    tauMTPass->SetBinError(3,0.000607802);
    tauMTPass->SetBinError(4,0.001363476);
    tauMTPass->SetBinError(5,0.010027);
    tauMTPass->SetBinError(6,0.02594956);

    Cout(p++, tauMTPass);

    for(int j = 0; j < tauMTPass->GetNbinsX(); j++){
      float value = tauMTPass->GetBinContent(j+1);

      float error = tauMTPass->GetBinError(j+1);
      
      float newError = sqrt(error * error + sysTauMT * sysTauMT * value * value);

      tauMTPass->SetBinError(j+1, newError);
    }

    TString fileName = fOutputDir;
    if(!fileName.EndsWith("/")) fileName += "/";
    Util::MakeOutputDir(fileName);
    fileName = fileName  + "tauMTPassCompared.root";
    TFile *savefile = new TFile(fileName.Data(), "RECREATE");
    savefile ->cd();
  
    tauMTPass->Write();

    for(int j = 0; j < NumberOfSamples; j++)
      MT2[j]->Write();
  
    savefile->Close();
    std::cout << "Saved histograms in " << savefile->GetName() << std::endl;

    Cout(p++, MT2[NumberOfSamples]); 
    /*
    MT2[NumberOfSamples]->Multiply(tauMTPass);
    EstimPRDown->Multiply(tauMTPass);
    EstimPRUp->Multiply(tauMTPass);
    EstimFRDown->Multiply(tauMTPass);
    EstimFRUp->Multiply(tauMTPass);
    Sigma_wi2->Multiply(tauMTPass);
    */
  }else{
    TString fileName = fOutputDir;
    if(!fileName.EndsWith("/")) fileName += "/";
    fileName = fileName + "tauMTPass.root";
    TFile *file = new TFile(fileName.Data(), "READ");
 
    tauMTPass = (TH1F*) file->Get("tauMTPass");

    MT2[0] = (TH1*) file->Get("MT2_QCD");
    MT2[1] = (TH1*) file->Get("MT2_W");
    MT2[2] = (TH1*) file->Get("MT2_ZX");
    MT2[3] = (TH1*) file->Get("MT2_Top");
    MT2[4] = (TH1*) file->Get("MT2_WW");
    MT2[5] = (TH1*) file->Get("MT2_Higgs");
    MT2[NumberOfSamples-2] = (TH1*) file->Get("MT2_MC");
    MT2[NumberOfSamples-1] = (TH1*) file->Get("MT2_susy");

    Cout(p++, tauMTPass);

    for(int j = 0; j < tauMTPass->GetNbinsX(); j++){
      float value = tauMTPass->GetBinContent(j+1);

      float error = tauMTPass->GetBinError(j+1);
      
      float newError = sqrt(error * error + sysTauMT * sysTauMT * value * value);

      tauMTPass->SetBinError(j+1, newError);
    }
    Cout(p++, tauMTPass);
    
    AddOverAndUnderFlow(MT2[NumberOfSamples], true, true);

    Cout(p++, MT2[NumberOfSamples]); 
    
    MT2[NumberOfSamples]->Multiply(tauMTPass);
    EstimPRDown->Multiply(tauMTPass);
    EstimPRUp->Multiply(tauMTPass);
    EstimFRDown->Multiply(tauMTPass);
    EstimFRUp->Multiply(tauMTPass);
    Sigma_wi2->Multiply(tauMTPass);
  }

  Cout(p++, MT2[NumberOfSamples]); 
  Cout(p++, EstimPRDown);
  Cout(p++, EstimFRUp);
  Cout(p++, EstimPRDown);
  Cout(p++, EstimFRUp);

  for(int j = 0; j < MT2[NumberOfSamples]->GetNbinsX(); j++){
    float totalError = MT2[NumberOfSamples]->GetBinError(j+1);
    float centralVal = MT2[NumberOfSamples]->GetBinContent(j+1);
    float PRDown = EstimPRDown->GetBinContent(j+1);
    float FRDown = EstimFRDown->GetBinContent(j+1);
    float PRUp = EstimPRUp->GetBinContent(j+1);
    float FRUp = EstimFRUp->GetBinContent(j+1);
    
    float PR = fabs(PRUp - PRDown)/2.0;//
    
    EstimPRDown->SetBinContent(j+1, 1.0);
    EstimPRDown->SetBinError(j+1, PR/centralVal);
 
    float FR = fabs(FRUp - FRDown)/2.0;//
    
    EstimFRDown->SetBinContent(j+1, 1.0);
    EstimFRDown->SetBinError(j+1, FR/centralVal);
    

    MT2EstSys->SetBinContent(j+1, MT2[NumberOfSamples]->GetBinContent(j+1));
    float StatisticalSystematicError = sqrt(totalError * totalError + PR * PR + FR * FR);//adding the error from FR and PR
    MT2[NumberOfSamples]->SetBinError(j+1, StatisticalSystematicError);
    
    float StatisticalError = Sigma_wi2->GetBinContent(j+1);
    
    MT2EstSys->SetBinError(j+1, StatisticalError);

  }

  Cout(p++, MT2EstSys);

  Cout(p++, MT2[NumberOfSamples]); 
  Cout(p++, EstimPRDown);
  Cout(p++, EstimFRDown);

  printYield();

  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j <= (NumberOfSamples); j++){
    if(j <= (NumberOfSamples - 3))
      h_stack -> Add(MT2[j]);
  }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W", "f");
  Legend1->AddEntry(MT2[2], "ZX", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW", "f");
  Legend1->AddEntry(MT2[5], "Higgs", "f");
  Legend1->AddEntry(MT2[NumberOfSamples-1], "SMS", "l");
  Legend1->AddEntry(MT2[NumberOfSamples], "data", "l");
 
  std::cout << " Fake Rate read from " << myfileName << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;
  
  TCanvas *MyC = new TCanvas("MyC","MyC");
  MyC->Divide(2,2);
  MyC->cd(1);
  hPt1Pass->Draw();
  MyC->cd(2);
  hPt2Pass->Draw();
  MyC->cd(3);
  AddOverAndUnderFlow(WeightsFakeRate);
  WeightsFakeRate->Draw();
  MyC->cd(4);
  AddOverAndUnderFlow(myWeights);
  myWeights->Draw();

  printHisto(h_stack, MT2[NumberOfSamples], MT2[NumberOfSamples-2], MT2[NumberOfSamples-1], Legend1 , "MTC", "hist", true, "MT2", "Events", 0, -10, 2, true);

  plotRatioStack(h_stack, MT2[NumberOfSamples-2], MT2[NumberOfSamples], MT2[NumberOfSamples-1], true, false, "MT2_ratio", Legend1, "MT2", "Events", 0, -10, 2, true, "C");


}


void MassPlotter::DrawMyPlots(TString myfileName, double *xbin, int NumberOfBins){

 TH1::SetDefaultSumw2();

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  fileName = fileName + myfileName;
  TFile *file = new TFile(fileName.Data(), "READ");
  
  TH1* MT2_MC = (TH1*) file->Get("MT2_MC");
  TH1* MT2_susy = (TH1*) file->Get("MT2_susy");
  TH1* MT2_data = (TH1*) file->Get("MT2_data");

  if(xbin != 0){
    TH1* MT2_Wjets = (TH1*) file->Get("MT2_Wjets");
    TH1* MT2_WWjets = (TH1*) file->Get("MT2_WWjets");
    TH1* MT2_Zjets = (TH1*) file->Get("MT2_Zjets");
    TH1* MT2_Top = (TH1*) file->Get("MT2_Top");
    TH1* MT2_QCD = (TH1*) file->Get("MT2_QCD");
  
    THStack* h_stack = new THStack("stack", "");
    TH1F* hnew_data = (TH1F*)MT2_data->Rebin(NumberOfBins,"hnew_data",xbin);
    TH1F* hnew_susy = (TH1F*)MT2_susy->Rebin(NumberOfBins,"hnew_susy",xbin);

    TH1F* hnew_MC = (TH1F*)MT2_MC->Rebin(NumberOfBins,"hnew_MC",xbin);
    hnew_MC->Scale(0);
    TH1F* hnew_QCD = (TH1F*)MT2_QCD->Rebin(NumberOfBins,"hnew_QCD",xbin);
    hnew_MC->Add(hnew_QCD);

    TH1F* hnew_Wjets = (TH1F*)MT2_Wjets->Rebin(NumberOfBins,"hnew_Wjets",xbin);
    hnew_MC->Add(hnew_Wjets);
    TH1F* hnew_WWjets = (TH1F*)MT2_WWjets->Rebin(NumberOfBins,"hnew_WWjets",xbin);
    hnew_MC->Add(hnew_WWjets);

    TH1F* hnew_Zjets = (TH1F*)MT2_Zjets->Rebin(NumberOfBins,"hnew_Zjets",xbin);
    hnew_MC->Add(hnew_Zjets);
    TH1F* hnew_Top = (TH1F*)MT2_Top->Rebin(NumberOfBins,"hnew_Top",xbin);
    hnew_MC->Add(hnew_Top);

    h_stack -> Add(hnew_QCD);
    h_stack -> Add(hnew_Wjets);
    h_stack -> Add(hnew_WWjets);
    h_stack -> Add(hnew_Zjets);
    h_stack -> Add(hnew_Top);

    TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
    Legend1->AddEntry(hnew_QCD, "QCD", "f");
    Legend1->AddEntry(hnew_Wjets, "W+jets", "f");
    Legend1->AddEntry(hnew_WWjets, "WW+jets", "f");
    Legend1->AddEntry(hnew_Zjets, "Z+jets", "f");
    Legend1->AddEntry(hnew_Top, "Top", "f");
    Legend1->AddEntry(hnew_susy, "SMS", "l");
    Legend1->AddEntry(hnew_data, "data", "l");


    cout << "------------------------------------" << endl
<< "QCD Integral: " << hnew_QCD->Integral() << endl
<< "W+jets Integral: " << hnew_Wjets->Integral() << endl
<< "WW+jets Integral: " << hnew_WWjets->Integral() << endl
<< "Z+jets Integral: " << hnew_Zjets->Integral() << endl
<< "Top Integral: " << hnew_Top->Integral() << endl
<< "TOTAL BG: " << hnew_QCD->Integral()+hnew_Wjets->Integral()+hnew_WWjets->Integral()+hnew_Zjets->Integral()+hnew_Top->Integral() <<endl
<< "SUSY: " << hnew_susy->Integral() << endl
<< "Data: " << hnew_data->Integral() << endl;
    
    //Same, but with errors
//     string OverallSamples[7] = {"QCD","W+jets","Z+Jets","Top","SUSY","Data","TOTAL BG"};
    
    TH1F* clone = (TH1F*) hnew_QCD->Clone();
    clone->Rebin(clone->GetNbinsX());
    cout << "QCD : " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;

    clone = (TH1F*) hnew_Wjets->Clone();
    clone->Rebin(clone->GetNbinsX());
    cout << "Wjets : " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;

    clone = (TH1F*) hnew_WWjets->Clone();
    clone->Rebin(clone->GetNbinsX());
    cout << "WWjets : " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;

    clone = (TH1F*) hnew_Zjets->Clone();
    clone->Rebin(clone->GetNbinsX());
    cout << "Zjets : " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;

    clone = (TH1F*) hnew_Top->Clone();
    clone->Rebin(clone->GetNbinsX());
    cout << "Top : " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
   
    clone = (TH1F*) hnew_MC->Clone();
    clone->Rebin(clone->GetNbinsX());
    cout << "Total BG : " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
    
    clone = (TH1F*) hnew_data->Clone();
    clone->Rebin(clone->GetNbinsX());
    cout << "data : " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
 
    clone = (TH1F*) hnew_susy->Clone();
    clone->Rebin(clone->GetNbinsX());
    cout << "SUSY : " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;

    printHisto(h_stack, hnew_data, hnew_MC, hnew_susy, Legend1 , "MTC", "hist", true, "MT2", "Events", 0, -10, 2, true);

    plotRatioStack(h_stack, hnew_MC, hnew_data, hnew_susy, true, false, "MT2_ratio", Legend1, "MT2", "Events", 0, -10, 2, true,"png");
  }else{
    THStack* h_stack = (THStack*) file->Get("MT2");
   
    TLegend* Legend1 = (TLegend*)file->Get("TPave");

    printHisto(h_stack, MT2_data, MT2_MC, MT2_susy, Legend1 , "MTC", "hist", true, "MT2", "Events", 0, -10, 2, true);

    plotRatioStack(h_stack, MT2_MC, MT2_data, MT2_susy, true, false, "MT2_ratio", Legend1, "MT2", "Events", 0, -10, 2, true);
  }


}



void MassPlotter::setLimitCounting(TString channels, TString cuts, TString hypothesis) {

  // === O U T - F I L E S ===                                                                                                                                                     
  TFile *fout = new TFile("counting.root", "RECREATE");
  fout->cd();

  // === B O O K I N G ===                                                                                                                                                         
  TH2D *h_PN_MLSP_MChi = new TH2D("h_PN_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
  h_PN_MLSP_MChi->Sumw2();

  TH1D *h_PN_Bkg = new TH1D("h_PN_Bkg", "", 1, 0, 1e6);
  h_PN_Bkg->Sumw2();


  // === F I L L I N G ===                                                                                                                                                         

  for (size_t i = 0; i < fSamples.size(); ++i) {

    cout << ">> " << fSamples[i].name << endl;

    //af: To remove inf                                                                                                                                                          
    if (fSamples[i].nevents == 0) fSamples[i].nevents = 1;

    TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
    TString ChannelSpecificSF = "1.00";
    if (fSamples[i].type != "data") {
      btagweight = "SFWeight.BTagCSV40eq0";

      if (channels == "muTau") {
	ChannelSpecificSF += " * muTau[0].muIdSF * muTau[0].muIsoSF * muTau[0].muTrgSF * muTau[0].tauTrgSF * muTau[0].tauEnergySF";
	if (fSamples[i].sname == "Wtolnu") ChannelSpecificSF += " * muTau[0].tauWjetsSF";
      } else if (channels == "eleTau") {
	ChannelSpecificSF += " * eleTau[0].tauTrgSF * eleTau[0].eleTrgSF * eleTau[0].eleIdIsoSF ";
	if (fSamples[i].sname == "Wtolnu") ChannelSpecificSF += " * eleTau[0].tauWjetsSF";
      } else if (channels == "eleMu") {
	ChannelSpecificSF += "* eleMu[0].MuIdIsoSF * eleMu[0].MuTrgSF * eleMu[0].EleTrgSF ";
      }
    }

    if (hypothesis == "sgn" && fSamples[i].type == "susy") {

      Double_t weight = fSamples[i].lumi;
      TString selection = TString::Format("(%.15f*%s*%s) * (%s)", weight, btagweight.Data(), ChannelSpecificSF.Data(), cuts.Data());
      TString variable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_PN_MLSP_MChi->GetName());

      fSamples[i].tree->Draw(variable, selection, "goff");

      h_SMSEvents = (TH2F*) fSamples[i].file->Get("h_SMSEvents");
      h_SMSEvents->Rebin2D(4, 4);
      h_PN_MLSP_MChi->Divide(h_SMSEvents);
      TH2* hXsec = (TH2*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");
      h_PN_MLSP_MChi->Multiply(hXsec);

    } else if (hypothesis == "bkg" && fSamples[i].type == "mc") { // WARNING: checking statistical fluctuation in this pt bin!!!                                                       
      Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents * fSamples[i].PU_avg_weight);
      TString selection = TString::Format("(%.15f*pileUp.Weight*%s*%s) * (%s)", weight, btagweight.Data(), ChannelSpecificSF.Data(), cuts.Data());
      TString variable = TString::Format("misc.MET>>+%s", h_PN_Bkg->GetName());

      fSamples[i].tree->Draw(variable, selection, "goff");


    }
  }//i fsamples.size()   

  fout->Write();
  fout->Close();
}


void MassPlotter::makeCard(double N, double S, double dS, double B, double dB, string sOut) {

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 1  number of backgrounds" << std::endl;
    fOut << "kmax 2  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
    fOut << "observation " << N << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1" << std::endl;
    fOut << "process         SMS    All" << std::endl;
    fOut << "process          0     1  " << std::endl;
    fOut << "rate           " << S << "\t" << B << std::endl;
    fOut << "---" << std::endl;
    fOut << "dS  lnN    " << 1 + dS << "\t-" << std::endl;
    fOut << "dB  lnN    - \t " << 1 + dB << std::endl;
    fOut.close();

}

TGraphAsymmErrors* MassPlotter::plotSig(TH1 *hSgn, TH1 *hBkg, TString xtitle, TString cutType, int type, double sys) {
    int nbins = hSgn->GetXaxis()->GetNbins();
    float *x = new float[nbins + 1];
    float *ex = new float[nbins + 1];
    float *y = new float[nbins + 1];
    float *ey = new float[nbins + 1];
    float *eyp = new float[nbins + 1];
    float *eym = new float[nbins + 1];

    for (int i = 1; i <= nbins+1; i++) {
 
        x[i - 1] = hSgn->GetBinLowEdge(i);
        ex[i - 1] = hSgn->GetBinWidth(i)/2.0;

        double s = (cutType == "Lower Cut") ? hSgn->Integral(i, nbins + 1) : hSgn->Integral(0, i-1);
        double ds = sqrt(s) + s * sys;
        double b = (cutType == "Lower Cut") ? hBkg->Integral(i, nbins + 1) : hBkg->Integral(0, i-1);
        double db = sqrt(b) + b * sys;


        if (b == 0 || s == 0) {
            y[i - 1] = .0;
            ey[i - 1] = .0;
        } else {
            if (type == 0) {
                y[i - 1] = s / sqrt(b);
                ey[i - 1] = y[i - 1] * (ds / s + db / (2 * b));
            }
            if (type == 1) {
                y[i - 1] = s / sqrt(s + b);
                ey[i - 1] = y[i - 1] * (ds / s + (db + ds) / (2 * (b + s)));
            }
            if (type == 2) {
                y[i - 1] = s / b;
                ey[i - 1] = y[i - 1] * (ds / s + db / b);
            }
            if (type == 3) {
	      cout<<cutType<<" bin "<<i<<endl;
                makeCard(b, s, sys, b, sys, "datacard");
                if (!(std::ifstream("datacard")).good()) continue;
                system("combine -M Asymptotic datacard");
                TTree* tree;

		TFile * flimit = new TFile("higgsCombineTest.Asymptotic.mH120.root");

                flimit->GetObject("limit", tree);

                Double_t limit;
                TBranch *b_limit; //!
                tree->SetBranchAddress("limit", &limit, &b_limit);

                Float_t quantileExpected;
                TBranch *b_quantileExpected; //!
                tree->SetBranchAddress("quantileExpected", &quantileExpected, &b_quantileExpected);

                std::vector<double> vLimit;
                Long64_t nEntrs = tree->GetEntriesFast();
                for (Long64_t iEntr = 0; iEntr < nEntrs; iEntr++) {
                    tree->GetEntry(iEntr);
                    cout << ">> quantileExpected: " << quantileExpected << "\tlimit: " << limit << endl;
                    vLimit.push_back(limit);
                }

                double SgmP2(vLimit[0]), SgmP1(vLimit[1]), Mdn(vLimit[2]), SgmM1(vLimit[3]), SgmM2(vLimit[4]), Obs(vLimit[5]);

                y[i - 1] = Mdn;
                eyp[i - 1] = SgmM1 - y[i - 1];
                eym[i - 1] = y[i - 1] - SgmP1;

		system("rm -f higgsCombineTest.Asymptotic.mH120.root");
		system("rm -f datacard");
                system("rm -f roostat*");

            }
        }
    }

    TGraphAsymmErrors *sig = new TGraphAsymmErrors(nbins +1 , x, y, ex, ey);
    if (type == 3) sig = new TGraphAsymmErrors(nbins+1 , x, y, ex, ex, eym, eyp);

    sig->SetTitle("");
    sig->GetXaxis()->SetTitle(xtitle + "_" + cutType);
    sig->SetMarkerStyle(20);
    sig->SetFillColor(kBlue-7);
    sig->SetFillStyle(3005);
    if (type == 0) sig->GetYaxis()->SetTitle("S/#sqrt{B}");
    if (type == 1) sig->GetYaxis()->SetTitle("S/#sqrt{S+B}");
    if (type == 2) sig->GetYaxis()->SetTitle("S/B");
    if (type == 3) sig->GetYaxis()->SetTitle("signal strength (r)");

    return sig;
}

void MassPlotter::LeptonEfficiency(TString cuts, unsigned int nevents){
  TH1::SetDefaultSumw2();
 
  float xbin1[4] = {0.0, 50.0, 100.0, 1000.0}; 

  TH1F *hPtTauAllBarrel = new TH1F("PtTauAllBarrel","PtTauAllBarrel", 3, xbin1);
  TH1F *hPtMuAllBarrel  = new TH1F("PtMuAllBarrel", "PtMuAllBarrel",  3, xbin1);
  TH1F *hPtEleAllBarrel = new TH1F("PtEleAllBarrel","PtEleAllBarrel", 3, xbin1);

  TH1F *hPtTauAllEndcap = new TH1F("PtTauAllEndcap","PtTauAllEndcap", 3, xbin1);
  TH1F *hPtMuAllEndcap  = new TH1F("PtMuAllEndcap", "PtMuAllEndcap",  3, xbin1);
  TH1F *hPtEleAllEndcap = new TH1F("PtEleAllEndcap","PtEleAllEndcap", 3, xbin1);

  TH1F *hPtTauPassBarrel = new TH1F("PtTauPassBarrel","PtTauPassBarrel", 3, xbin1);
  TH1F *hPtMuPassBarrel  = new TH1F("PtMuPassBarrel", "PtMuPassBarrel",  3, xbin1);
  TH1F *hPtElePassBarrel = new TH1F("PtElePassBarrel","PtElePassBarrel", 3, xbin1);

  TH1F *hPtTauPassEndcap = new TH1F("PtTauPassEndcap","PtTauPassEndcap", 3, xbin1);
  TH1F *hPtMuPassEndcap  = new TH1F("PtMuPassEndcap", "PtMuPassEndcap",  3, xbin1);
  TH1F *hPtElePassEndcap = new TH1F("PtElePassEndcap","PtElePassEndcap", 3, xbin1);

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
   TString myCuts = cuts;
 

   sample Sample = fSamples[ii];
    
   if(Sample.shapename != "ZJetsToLL" && Sample.name != "TTbar" && Sample.sname != "SUSY")
     continue;

      fMT2tree = new MT2tree();
   Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

   float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

   std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
   if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){

       Weight = Sample.lumi;
       if(fPUReweight) Weight /= Sample.PU_avg_weight;
       
   }

   Sample.tree->Draw(">>selList", myCuts);
   TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
   
   unsigned int nentries =  myEvtList->GetN();//Sample.tree->GetEntries();
   
   for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
     
     Sample.tree->GetEntry(myEvtList->GetEntry(jentry));
     
     if ( fVerbose>2 && jentry % 100000 == 0 ){  
       fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
       fflush(stdout);
     }
     
     float weight = Weight * fMT2tree->pileUp.Weight;
     
     for(int j = 0; j < fMT2tree->NGenLepts; j++){
       if(abs(fMT2tree->genlept[j].ID) != 15 && abs(fMT2tree->genlept[j].ID) != 13 && abs(fMT2tree->genlept[j].ID) != 11)
	 continue;
		
	int HadTau = 1;
	
	if(abs(fMT2tree->genlept[j].ID) == 15){   
	  for(int i = 0; i < fMT2tree->NGenLepts; i++){
	    if((abs(fMT2tree->genlept[i].ID) == 11 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID) || 
	       (abs(fMT2tree->genlept[i].ID) == 13 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID)){
	      HadTau = 0;
	    }
	  }
	}

 	if(HadTau == 0)
	  continue;
	
	TLorentzVector genLepLV;

	genLepLV = fMT2tree->genlept[j].lv;
	

	//To find the recLep matched to this genLep
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();

	float genleptPhi = fMT2tree->genlept[j].lv.Phi();
	
	//To find the rec/id tau matched to this gen had tau
	
	if(abs(fMT2tree->genlept[j].ID) == 15){ 
	  for(int k = 0; k < fMT2tree->NTaus; k++){
	    float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->tau[k].lv.Eta(), genleptPhi, fMT2tree->tau[k].lv.Phi());
	    if(fMT2tree->tau[k].PassTau_MuTau == 1)
	      if(deltaR < minDR){
		minDR = deltaR;
	      }
	  }
	  if(fabs(genleptEta) < 1.479)
	    hPtTauAllBarrel->Fill(genLepLV.Pt(), weight);
	  else
	    hPtTauAllEndcap->Fill(genLepLV.Pt(), weight);  
	}

	if(abs(fMT2tree->genlept[j].ID) == 13){ 
	  for(int k = 0; k < fMT2tree->NMuons; k++){
	    float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->muo[k].lv.Eta(), genleptPhi, fMT2tree->muo[k].lv.Phi());
	    if(fMT2tree->muo[k].PassMu0_TauMu == 1)
	      if(deltaR < minDR){
		minDR = deltaR;
	      }
	  }
	  if(fabs(genleptEta) < 1.479)
	    hPtMuAllBarrel->Fill(genLepLV.Pt(), weight);
	  else
	    hPtMuAllEndcap->Fill(genLepLV.Pt(), weight);  
	}

	if(abs(fMT2tree->genlept[j].ID) == 11){ 
	  for(int k = 0; k < fMT2tree->NEles; k++){
	    float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->ele[k].lv.Eta(), genleptPhi, fMT2tree->ele[k].lv.Phi());
	    if(fMT2tree->ele[k].IDSelETau == 1)
	      if(deltaR < minDR){
		minDR = deltaR;
	      }
	  }
	  if(fabs(genleptEta) < 1.479)
	    hPtEleAllBarrel->Fill(genLepLV.Pt(), weight);
	  else
	    hPtEleAllEndcap->Fill(genLepLV.Pt(), weight);  
	}

	if(minDR < 0.1){
	  if(abs(fMT2tree->genlept[j].ID) == 15){
	    if(fabs(genleptEta) < 1.479)
	      hPtTauPassBarrel->Fill(genLepLV.Pt(), weight);
	    else
	      hPtTauPassEndcap->Fill(genLepLV.Pt(), weight);}

	  if(abs(fMT2tree->genlept[j].ID) == 13){
	    if(fabs(genleptEta) < 1.479)
	      hPtMuPassBarrel->Fill(genLepLV.Pt(), weight);
	    else
	      hPtMuPassEndcap->Fill(genLepLV.Pt(), weight);}
	 
	  if(abs(fMT2tree->genlept[j].ID) == 11){
	    if(fabs(genleptEta) < 1.479)
	      hPtElePassBarrel->Fill(genLepLV.Pt(), weight);
	    else
	      hPtElePassEndcap->Fill(genLepLV.Pt(), weight);}
	}
     }
    }
  }

  hPtTauPassBarrel->Divide(hPtTauAllBarrel);
  hPtMuPassBarrel ->Divide(hPtMuAllBarrel);
  hPtElePassBarrel->Divide(hPtEleAllBarrel);

  hPtTauPassEndcap->Divide(hPtTauAllEndcap);
  hPtMuPassEndcap ->Divide(hPtMuAllEndcap);
  hPtElePassEndcap->Divide(hPtEleAllEndcap);
  
  TCanvas *MyC = new TCanvas("MyC", "MyC");
  MyC->Divide(3,4);
  MyC->cd(1);
  hPtEleAllBarrel->Draw();
  MyC->cd(2);
  hPtMuAllBarrel->Draw();
  MyC->cd(3);
  hPtTauAllBarrel->Draw();
  MyC->cd(4);
  hPtElePassBarrel->Draw();
  MyC->cd(5);
  hPtMuPassBarrel->Draw();
  MyC->cd(6);
  hPtTauPassBarrel->Draw();
  MyC->cd(7);
  hPtEleAllEndcap->Draw();
  MyC->cd(8);
  hPtMuAllEndcap->Draw();
  MyC->cd(9);
  hPtTauAllEndcap->Draw();
  MyC->cd(10);
  hPtElePassEndcap->Draw();
  MyC->cd(11);
  hPtMuPassEndcap->Draw();
  MyC->cd(12);
  hPtTauPassEndcap->Draw();
}


void MassPlotter::eeFakeRateRatio(TString cuts, TString trigger, unsigned int nevents, TString myfileName){

 TH1::SetDefaultSumw2();
 
 TH2F *hPtEtaAll  = new TH2F("hPtEtaAll", "hPtEtaAll",  60, -3.0, 3.0, 1000, 0, 1000);
 TH2F *hPtEtaPass = new TH2F("hPtEtaPass","hPtEtaPass", 60, -3.0, 3.0, 1000, 0, 1000);

 float xBins[4] = {0, 50, 100, 300};

 TH1F *hMETAll   = new TH1F("hMETAll"  ,"hMETAll"  ,3, xBins);
 TH1F *hMETPass  = new TH1F("hMETPass" ,"hMETPass" ,3, xBins);
 TH1F *hElePtAll  = new TH1F("hElePtAll" ,"hElePtAll" ,3, xBins);
 TH1F *hElePtPass = new TH1F("hElePtPass","hElePtPass",3, xBins);
 TH1F *hMT2All   = new TH1F("hMT2All"  ,"hMT2All"  ,3, xBins);
 TH1F *hMT2Pass  = new TH1F("hMT2Pass" ,"hMT2Pass" ,3, xBins);

 for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
   TString myCuts = cuts;
 
   int data = 0;
   sample Sample = fSamples[ii];
    
   if(Sample.type == "data"){
     //     if(Sample.sname == "Wtolnu"){
  data = 1;
     myCuts += " && " + trigger;
 }
 else 
 continue;

//           data = 1;                                                                                                                                  
   

   fMT2tree = new MT2tree();
   Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

   float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

   std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry); 
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

//          int genEEstatus = fMT2tree->GenLeptonAnalysis(2,0,0);
//              if(genEEstatus == 113 || genEEstatus == 123)
// 	       {data=0;}
//       	     else 
//                        continue;

     float weight = Weight;

      if(data == 1)
 	weight = 1.0;
     else{
	  weight *= (fMT2tree->pileUp.Weight * fMT2tree->SFWeight.BTagCSV40eq0/Sample.PU_avg_weight);
      }      

      //      if ((fMT2tree->NEles < 1 && fMT2tree->NJets == 0)/*|| fMT2tree->ele[0].lv.Pt() < 30*/ || fMT2tree->ele[0].PassE0_EE != 1)
	//      	continue;
      //if(!((fMT2tree->NEles == 2) && (fMT2tree->ele[0].lv.Pt() > 30) && (fMT2tree->ele[0].PassE0_EE == 1) && (fMT2tree->ele[0].Charge+fMT2tree->ele[1].Charge == 2 || fMT2tree->ele[0].Charge+fMT2tree->ele[1].Charge == -2)))

      //   continue;

      //      bool ok =false;
      //for (int kk=0; kk < fMT2tree->NEles; kk++)
	
	  //    for (int jj=0; jj < fMT2tree->NJets; jj++)
      //  {  float deltaR = Util::GetDeltaR(fMT2tree->ele[kk].lv.Eta(), fMT2tree->jet[0].lv.Eta(), fMT2tree->ele[kk].lv.Phi(), fMT2tree->jet[0].lv.Phi());
	  	
     //     cout << "deltaR..." << deltaR << endl;
    // if(deltaR > 0.5 && deltaR < 6 /*&& fMT2tree->ele[kk].PassE0_EE == 1 /*&& (fMT2tree->ele[kk].MT < 70 || fMT2tree->ele[kk].MT > 90)*/){

	       

	//  if((fMT2tree->ele[0].Charge)+(fMT2tree->ele[kk].Charge) != 2 || (fMT2tree->ele[0].Charge)+(fMT2tree->ele[kk].Charge) != -2)

	  //	if((fMT2tree->ele[0].MT) < 60 || (fMT2tree->ele[0].MT) > 90); 
	  //        {continue;}
		
		//		if (kk==0)
		//		  {myMT2 =  fMT2tree->CalcMT2(0, 0, fMT2tree->ele[0].lv, fMT2tree->ele[1].lv, fMT2tree->pfmet[0]);}//fMT2tree->doubleEle[0].MT2;
		//		else
//    float myMT2=0;
//    if(kk==0)
//        myMT2 = fMT2tree->CalcMT2(0, 0, fMT2tree->ele[0].lv, fMT2tree->ele[1].lv, fMT2tree->pfmet[0]);
//    else
//        myMT2 = fMT2tree->CalcMT2(0, 0, fMT2tree->ele[0].lv, fMT2tree->ele[kk].lv, fMT2tree->pfmet[0]);
 //myChannelCuts.push_back("ele[1].QCDSyst03E0_EE == 1");                                                                                            
      //if(deltaR < 0.5 && deltaR > 6) /*&& (fMT2tree->ele[kk].MT < 70 || fMT2tree->ele[kk].MT > 90)*/
      //  continue;
  
// if(fMT2tree->ele[kk].PassE0_EE == 1)
//   {
//      for(int jj=0; jj < fMT2tree->NEles; jj++)
//        {
//          if(jj != kk && fMT2tree->ele[jj].PassE0_EE == 1) 
//          ok = true;
//          break;
//        }
//   }
// if(ok)
//   break; 



//	  if (GenLeptFromW(11, 10, 2.4, true) == true)



 //        if(fMT2tree->ele[0].QCDSyst03E0_EE == 1)
        float myMT2 = fMT2tree->CalcMT2(0, 0, fMT2tree->ele[0].lv, fMT2tree->ele[1].lv, fMT2tree->pfmet[0]);

      	hPtEtaAll->Fill(fMT2tree->ele[1].lv.Eta(), fMT2tree->ele[1].lv.Pt()); 
      	hMETAll->Fill(fMT2tree->misc.MET);
      	hMT2All->Fill(myMT2);
	hElePtAll->Fill(fMT2tree->ele[1].lv.Pt());
      
	if(fMT2tree->ele[1].PassE1_EE == 1){

        hPtEtaPass->Fill(fMT2tree->ele[1].lv.Eta(), fMT2tree->ele[1].lv.Pt());
	hMETPass->Fill(fMT2tree->misc.MET);
	hMT2Pass->Fill(myMT2);
	hElePtPass->Fill(fMT2tree->ele[1].lv.Pt());      
	}

	//      }
	//	}
  
     
 }
 }
 TCanvas *MyC = new TCanvas("Fake","Fake");
 MyC->Divide(4,2);

 MyC->cd(1);
 hElePtAll->Draw();

 MyC->cd(2);
 hMETAll->Draw();

 MyC->cd(3);
 hMT2All->Draw();

 MyC->cd(4);
 hPtEtaAll->Draw();

 MyC->cd(5);
 hElePtPass->Divide(hElePtAll);
 hElePtPass->Draw();

 MyC->cd(6);
 hMETPass->Divide(hMETAll);
 hMETPass->Draw();

 MyC->cd(7);
 hMT2Pass->Divide(hMT2All);
 hMT2Pass->Draw();

 MyC->cd(8);
 hPtEtaPass->Divide(hPtEtaAll);
 hPtEtaPass->Draw();

 TString fileName = fOutputDir;
 if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName  + myfileName +"_FRHistos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();

  hElePtAll->Write();
  hMETAll->Write();
  hMT2All->Write();
  hPtEtaAll->Write();

  hElePtPass->Write();
  hMETPass->Write();
  hMT2Pass->Write();
  hPtEtaPass->Write();
  MyC->Write();

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;
}

void MassPlotter::eeWJetsEstimation(TString cuts, TString trigger, TString myfileName){

  TH1::SetDefaultSumw2();
  //  setFlags(10);
  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  fileName = fileName + myfileName;
  TFile *file = new TFile(fileName.Data(), "READ");

  TH2* hPtEtaAll  = (TH2*) file->Get("hPtEtaAll");
  TH2* hPtEtaPass = (TH2*) file->Get("hPtEtaPass");

  float xbin1[4] = {0.0, 50.0, 100.0, 1000.0}; //Barrel
  TH1F* hPt1All  = new TH1F("hPt1All",  "Barrel",  3,xbin1);
  TH1F* hPt1Pass = new TH1F("hPt1Pass", "Barrel", 3,xbin1);

  float xbin2[4] = {0.0, 50.0, 100.0, 1000.0}; //Endcap
  TH1F* hPt2All  = new TH1F("hPt2All",  "Endcap",  3,xbin2);
  TH1F* hPt2Pass = new TH1F("hPt2Pass", "Endcap", 3,xbin2);
 

  for(int i = 1; i < 61; i++){
    for(int j = 1; j < 1001; j++){

      double binIJAll  = hPtEtaAll ->GetBinContent(i,j);
      double binIJPass = hPtEtaPass->GetBinContent(i,j);
      
      if(i < 17 || i > 44){//Endcap
	hPt2All ->Fill(j - 0.5, binIJAll);
	hPt2Pass->Fill(j - 0.5, binIJPass);
      }else{
	hPt1All ->Fill(j - 0.5, binIJAll);
	hPt1Pass->Fill(j - 0.5, binIJPass);
      }
    }
  }

  hPt1Pass->Divide(hPt1All);
  
  hPt2Pass->Divide(hPt2All);



      static const int nbins = 6;//11;
//   //  double bins[nbins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,125.0,150.0,175.0,200.0,250.0,300.0,400.0};      //MT2
//   //  double bins[nbins+1] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0,120.0,200.0};      //MT2
   double bins[nbins+1] = {50.0,60.0,70.0,80.0,100.0,120.0,200.0};      //MT2

//   fileName = fOutputDir + "/FakePromptRatio-Single_FRHistos.root";//MT2.root"; "/MT2_MuTau_Over_QCDMuTau_SignalSelectionNoZVeto_PRHistos.root";
  
//   TFile *file2 = new TFile(fileName.Data(), "READ");

//   TH1* hOldMT2All  = (TH1*) file2->Get("hPtEtaAll");
//   TH1* hOldMT2Pass = (TH1*) file2->Get("hPtEtaPass");

// //   TH1* hMT2All   = hOldMT2All ->Rebin(nbins, "hMT2All" , bins);
// //   TH1* hMT2Pass  = hOldMT2Pass->Rebin(nbins, "hMT2Pass", bins);
//   TH1* hMT2All   = hOldMT2All->Rebin(10, "hMT2All");
//   TH1* hMT2Pass  = hOldMT2Pass->Rebin(10, "hMT2Pass");

 
 
//   for(int j = 1; j < hMT2All->GetNbinsX(); j++){
//     cout<<" bin "<<j<<" hMT2All  "<<hMT2All ->GetBinContent(j)<<" +- "<<hMT2All ->GetBinError(j)<<endl;
//     cout<<" bin "<<j<<" hMT2Pass "<<hMT2Pass->GetBinContent(j)<<" +- "<<hMT2Pass->GetBinError(j)<<endl;
//   }

//   TH1F* hMT2All  = new TH1F("hMT2All",  "All",  nbins, bins);
//   TH1F* hMT2Pass = new TH1F("hMT2Pass", "Pass", nbins, bins);

//   for(int j = 1; j < 201; j++){
    
//     double binJAll  = hOldMT2All ->GetBinContent(j);
//     double binJPass = hOldMT2Pass->GetBinContent(j);
    
//     hMT2All ->Fill(j - 0.5, binJAll);
//     hMT2Pass->Fill(j - 0.5, binJPass);
             
//   }
  
// //  hMT2Pass->Divide(hMT2All);
//   for(int j = 1; j < hMT2All->GetNbinsX(); j++){
//     cout<<" bin "<<j<<" hMT2Pass "<<hMT2Pass->GetBinContent(j)<<" +- "<<hMT2Pass->GetBinError(j)<<endl;
//   }


  TH1F *myWeights       = new TH1F("myWeights",       "myWeights",       30, -0.1, 0.2);
  TH2F *WeightsFakeRate = new TH2F("WeightsFakeRate", "WeightsFakeRate", 30, 0, 0.025, 30, -0.1, 0.2);
 
  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "WWjets", "MC", "susy","data"};
  int      ccolor[NumberOfSamples] = { 401,       417,    419,   855,       603,  603,      1, 632};
  TString varname = "MT2";
  for (int i=0; i<(NumberOfSamples); i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i], "", nbins, bins);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[7] -> SetMarkerStyle(20);
  MT2[7] -> SetMarkerColor(kBlack);
  MT2[7] -> SetLineColor(kBlack);
  
  MT2[5] -> SetFillStyle(3004);
  MT2[5] -> SetFillColor(kBlack);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  // vector of all histos
  vector<TH1D*> h_samples;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){


    TString myCuts = cuts;
    
    sample Sample = fSamples[ii];
    
 
    h_samples.push_back(new TH1D(varname+"_"+Sample.name, "", nbins, bins));
    h_samples[ii] -> Sumw2();
    h_samples[ii] -> SetLineColor(Sample.color);
   
    h_samples[ii] -> SetMarkerColor(Sample.color);
    h_samples[ii] -> SetStats(false);
    if(Sample.type == "susy" ){
      h_samples[ii] -> SetLineColor(kBlack);
      h_samples[ii] -> SetLineStyle(kDotted);
    }

    // myChannelCuts.push_back(std::string(std::string(myChan) + ".Ele0Ind>=0"));                                                                         
    // myChannelCuts.push_back(std::string(std::string(myChan) + ".Ele1Ind>=0"));                                                                         
    // myChannelCuts.push_back(std::string(std::string(myChan) + ".Isolated==1"));                                                                        
    // myChannelCuts.push_back("((ele[doubleEle[0].Ele0Ind].Charge + ele[doubleEle[0].Ele1Ind].Charge) == 0)");                                           
    // myChannelCuts.push_back("((doubleEle[0].lv.M() > 15 && doubleEle[0].lv.M() < 76) || (doubleEle[0].lv.M() > 106))");  

    if(Sample.sname != "Wtolnu")
      continue;

    if(Sample.type == "data"){
      myCuts += "&&" + trigger;
    }
    //else
    //      myCuts += "&& (doubleEle[0].Isolated == 1)";// && (tau[muTau[0].tau0Ind].Isolation3Hits <= 3.)";

  
    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

//     if (Sample.name == "QCD-Pt-15-20-MuEnriched"){
//       fPUReweight = false;}
//     else {fPUReweight = true;}

    float Weight =0;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    //    if(Sample.type == "data"){
    if(Sample.sname == "Wtolnu"){

      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
      
      unsigned int nentries = myEvtList->GetN();
      
      for (unsigned int jentry=0; jentry<nentries; jentry++) {
	
	Sample.tree->GetEntry(myEvtList->GetEntry(jentry));
	
	if ( fVerbose>2 && jentry % 100000 == 0 ){
	  fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	  fflush(stdout);
	}
	//	int muIndex = fMT2tree->muTau[0].GetMuIndex0();
	
	if (fMT2tree->NEles < 2)
	  continue;

	int Ele0Ind = -2;
	int Ele1Ind = -2;

	for (int i=0 ; i < fMT2tree->NEles ; i++){

	  if  (fMT2tree->ele[i].PassE0_EE)
	    {
	      Ele0Ind = i ;
	      break;
	    }
	}

	for (int j = Ele0Ind+1 ; j < fMT2tree->NEles ; j++)

	    {
              int PairCharge = fMT2tree->ele[Ele0Ind].Charge + fMT2tree->ele[j].Charge;
              if (PairCharge == 0)
		{
                 Ele1Ind = j ;
	         break;
		}
	    }

 if (Ele0Ind==-2 || Ele1Ind == -2)
     continue;

 TLorentzVector InvMasslv = fMT2tree->ele[Ele0Ind].lv+fMT2tree->ele[Ele1Ind].lv;
 float InvMass = InvMasslv.M();
 if(InvMass > 76 && InvMass < 106)
   continue;

	if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){

	  Weight = Sample.lumi;
	  if(fPUReweight) Weight /= Sample.PU_avg_weight;
	  
	}

	float weight = Weight;
		

	float Ele1Eta = fMT2tree->ele[Ele1Ind].lv.Eta();

	float Ele1Pt  = fMT2tree->ele[Ele1Ind].lv.Pt();

	float fakeRate    = 0.0;
// 	float fakeRateErr = 0.0;

	if(fabs(Ele1Eta) < 1.479){
	  int binNumber = hPt1Pass->FindBin(Ele1Pt);
	  fakeRate      = hPt1Pass->GetBinContent(binNumber);
// 	  fakeRateErr   = hPt1Pass->GetBinError(binNumber);
	}else{

	  int binNumber = hPt2Pass->FindBin(Ele1Pt);
	  fakeRate      = hPt2Pass->GetBinContent(binNumber);
// 	  fakeRateErr   = hPt2Pass->GetBinError(binNumber);
	}
	// P +  F = Loose
	//pP + fF = Loose - LooseNonTight
	//FakeContribution = f * F 
	//F * (f - p) = (1 - p) Loose - LooseNonTight

	float promptRate = 0.9;
	//	int binNumber = hMT2Pass->FindBin(myQuantity);
	//	promptRate    = hMT2Pass->GetBinContent(binNumber);
	
	weight = fakeRate * (1 - promptRate);

	//	bool IsoNonTightTrue = false;
	//        float IsoNonTight = fMT2tree->ele[Ele1Ind].Iso04;
	//	if (Ele1Eta < 1.479){
	//	  if(IsoNonTight  >= 0.15 )
	//            IsoNonTightTrue = true;
	//	  	}
	//        else{
	//	  if(IsoNonTight >= 0.10)
	//	IsoNonTightTrue = true;
	//	}

	if (!(fMT2tree->ele[Ele1Ind].PassE1_EE))
	  {weight -= fakeRate;}

	weight /= (fakeRate - promptRate);


	myWeights->Fill(weight);
	
	WeightsFakeRate->Fill(fakeRate, weight);

	float myQuantity =  fMT2tree->ele[Ele1Ind].lv.Pt();//CalcMT2(0, 0, fMT2tree->ele[Ele0Ind].lv, fMT2tree->ele[Ele1Ind].lv, fMT2tree->pfmet[0]);

	MT2[7]->Fill(myQuantity, weight);//data
      }
    }
    //else{
        if(Sample.sname == "Wtolnu"){
    
      //       if(Sample.sname != "Wtolnu" && Sample.sname != "QCD" && Sample.sname != "SUSY") 
      //	 continue;

      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

      unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();

      for (unsigned int jentry=0; jentry<nentries;jentry++) {
    
	Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

	if ( fVerbose>2 && jentry % 100000 == 0 ){
	  fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	  fflush(stdout);
	}

	if(fMT2tree->NEles < 2)
	  continue;

	int Ele0IndMC = -2;
	int Ele1IndMC = -2;

	for (int i=0 ; i < fMT2tree->NEles ; i++){

	  if  (fMT2tree->ele[i].PassE0_EE)
	    {
	      Ele0IndMC = i ;
	      break;
	    }
	}
	for (int jj = Ele0IndMC+1 ; jj < fMT2tree->NEles ; jj++){

              int PairChargeMC = fMT2tree->ele[Ele0IndMC].Charge + fMT2tree->ele[jj].Charge;
              if (fMT2tree->ele[jj].PassE1_EE && PairChargeMC == 0)
		{
                 Ele1IndMC = jj ;
	         break;
		}
	}

 if (Ele0IndMC==-2 || Ele1IndMC == -2)
     continue;
 //int PairChargeMC = fMT2tree->ele[Ele0IndMC].Charge + fMT2tree->ele[Ele1IndMC].Charge;
 // if (PairChargeMC != 0)
 //  continue;
 TLorentzVector InvMasslvMC = fMT2tree->ele[Ele0IndMC].lv+fMT2tree->ele[Ele1IndMC].lv;
 float InvMassMC = InvMasslvMC.M();
 if(InvMassMC > 76 && InvMassMC < 106)
   continue;

// 	int genEEStatus = fMT2tree->GenLeptonAnalysis(2,0,0);
// 	if(genEEStatus != 113 && genEEStatus != 123 && genEEStatus != 133 && genEEStatus != -1)
// 	  continue;       


	if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){

	  Weight = Sample.lumi;
	  if(fPUReweight) Weight /= Sample.PU_avg_weight;
	  
	}
	

	float weight = Weight;

        //////////////////////////////////
		if(fPUReweight)
	   weight *= fMT2tree->pileUp.Weight;
	/////////////////////////////////
	
	//	cout << "before.. "<< weight << endl; 

	//float DoubleEleTrgSF= ((fMT2tree->ele[Ele0IndMC].GetEleTrigSFleg17() * fMT2tree->ele[Ele1IndMC].GetEleTrigSFleg8()) + (fMT2tree->ele[Ele0IndMC].GetEleTrigSFleg8()  * fMT2tree->ele[Ele1IndMC].GetEleTrigSFleg17()) - (fMT2tree->ele[Ele0IndMC].GetEleTrigSFleg17() * fMT2tree->ele[Ele1IndMC].GetEleTrigSFleg17())) ;

	weight *=fMT2tree->SFWeight.BTagCSV40eq0;// * fMT2tree->ele[Ele0IndMC].EleIDIsoSFEleEle * fMT2tree->ele[Ele1IndMC].EleIDIsoSFEleEle * DoubleEleTrgSF;

// cout <<"after.. "<< weight <<endl;

	float myQuantityMC =  fMT2tree->ele[Ele1IndMC].lv.Pt();//CalcMT2(0, 0, fMT2tree->ele[Ele0IndMC].lv, fMT2tree->ele[Ele1IndMC].lv, fMT2tree->pfmet[0]);
	h_samples[ii]->Fill(myQuantityMC, weight);
      
      }

      AddOverAndUnderFlow(h_samples[ii]);

      if(Sample.type != "susy")
	MT2[5]->Add(h_samples[ii]);

      if (Sample.sname == "QCD") {
	MT2[0]->Add(h_samples[ii]);
      }
      else if(Sample.sname == "Top" ){
	MT2[3]->Add(h_samples[ii]);
      }
      else if(Sample.sname=="Wtolnu") {
	MT2[1]->Add(h_samples[ii]);
      }
      else if(Sample.sname=="DY") {
	MT2[2]->Add(h_samples[ii]);
      }
      else if(Sample.sname == "VV") {
	MT2[4]->Add(h_samples[ii]);
      }
      else if(Sample.type == "susy"){
	MT2[6]->Add(h_samples[ii]);

      }

    }
  }
 
 for(int j = 0; j < (NumberOfSamples); j++){
    AddOverAndUnderFlow(MT2[j], true, true);
  }


  printYield();

  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j < (NumberOfSamples); j++){
    // MT2[j]->Rebin(3);
    if(j < (NumberOfSamples - 3))
      h_stack -> Add(MT2[j]);
  }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW+jets", "f");
  Legend1->AddEntry(MT2[6], "SMS", "l");
  Legend1->AddEntry(MT2[7], "data", "l");
 
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;
  
  TCanvas *MyC = new TCanvas("MyC","MyC");
  MyC->Divide(2,1);
  MyC->cd(1);
  hPt1Pass->Draw();
  MyC->cd(2);
  hPt2Pass->Draw();
  //MyC->cd(3);
  //hMT2Pass->Draw();
  //  MyC->cd(4);
  //AddOverAndUnderFlow(WeightsFakeRate);
  //WeightsFakeRate->Draw();
  //AddOverAndUnderFlow(myWeights);
  //myWeights->Draw();


  myfileName =  fileName.ReplaceAll("_FRHistos.root", "");

  printHisto(h_stack, MT2[7], MT2[5], MT2[6], Legend1 , myfileName, "hist", true, myfileName , "Events", 0, -10, 2, true);

  plotRatioStack(h_stack, MT2[5], MT2[7], MT2[6], true, false, myfileName, Legend1, myfileName , "Events", 0, -10, 2, true);

}

void MassPlotter::eeFakePromptCategory(TString cuts, TString trigger, unsigned int nevents, TString myfileName){

  TH1::SetDefaultSumw2();

  TString  cnames[NumberOfSamples+1] = {"QCD", "Wjets", "Zjets", "Top", "WWjets", "MC", "susy","data"};
  int      ccolor[NumberOfSamples+1] = { 401,       417,    419,   855,       603,  603,      1, 632};
  TString varname = "MT2";
  for (int i=0; i<(NumberOfSamples); i++){
    //    MT2[i] = new TH1D(varname+"_"+cnames[i], "", 50, 0, 250);
    //following to check the promptness and faking of the events
    MT2[i] = new TH1D(varname+"_"+cnames[i], "", 8, 0, 8);
    
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);

    //following to check the promptness and faking of the events: 
    TString  genStatus[8] = {"pp", "pf", "ff", "pt", "tf", "tt" ,"nothing", "wrong"};
    for(int k = 0; k < MT2[i]->GetNbinsX(); k++)
      MT2[i] ->GetXaxis()->SetBinLabel(k+1,genStatus[k]);
      
}


  MT2[7] -> SetMarkerStyle(20);
  MT2[7] -> SetMarkerColor(kBlack);
  MT2[7] -> SetLineColor(kBlack);
  
  MT2[5] -> SetFillStyle(3004);
  MT2[5] -> SetFillColor(kBlack);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data"){
      data = 1;
      myCuts += " && " + trigger;
    }

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    if (Sample.name == "QCD-Pt-15-20-MuEnriched"){
      fPUReweight = false;}
    else {fPUReweight = true;}

    
    float Weight;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry);
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

      if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){

	Weight = Sample.lumi;
	if(fPUReweight) Weight /= Sample.PU_avg_weight;
      
      }
 
      float weight = Weight;

      if(data == 1)
	weight = 1.0;
      else{
	weight *= fMT2tree->SFWeight.BTagCSV40eq0 ;//* fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;
	
	if(fPUReweight)
	  weight *= fMT2tree->pileUp.Weight;
	
      }
 
//    bool ok =false;

// for (int kk=0; kk < fMT2tree->NEles; kk++)
// 	{
//     for (int jj=0; jj < fMT2tree->NJets; jj++)
//       {  float deltaR = Util::GetDeltaR(fMT2tree->ele[kk].lv.Eta(), fMT2tree->jet[jj].lv.Eta(), fMT2tree->ele[kk].lv.Phi(), fMT2tree->jet[jj].lv.Phi());
	  	
// cout << "deltaR..." << deltaR << endl;
//  if(deltaR > 0.5 && deltaR < 5)
//                  ok =true;
// 	       // continue;		 
	              
//               else
// 		ok =false;
// 	       break;      
// }

//                 if (!ok)
// 		  continue;
//                  else		  

// {
//   if(fMT2tree->ele[kk].PassE0_EE == 1){


//following to check the promptness and fakeness of the events:


// bool ok =false;      

// for (int ll=0; ll < fMT2tree->NMuons; ll++)
// {
//   //float deltaR1 = Util::GetDeltaR(fMT2tree->muo[ll].lv.Eta(), fMT2tree->jet[0].lv.Eta(), fMT2tree->muo[ll].lv.Phi(), fMT2tree->jet[0].lv.Phi());
// if(fMT2tree->muo[ll].PassMu0_EleMu ==1)
//     ok =true;
//     break;	
// }
// if(ok)
//     continue;
 

// for(int qq=1; qq < fMT2tree->NEles; qq++)
//   {
//     //float deltaR = Util::GetDeltaR(fMT2tree->ele[kk].lv.Eta(), fMT2tree->jet[0].lv.Eta(), fMT2tree->ele[kk].lv.Phi(), fMT2tree->jet[0].lv.Phi());

// if(fMT2tree->ele[qq].PassE0_EE == 1) //&& (fMT2tree->ele[kk].MT < 70) 
//   //  continue;

//   {
//        for(int jj=1; jj < fMT2tree->NEles; jj++)
//        {
//          if(jj != qq && fMT2tree->ele[jj].PassE0_EE == 1) 
//          ok = true;
//          break;
//        }
//   }

// if(ok)
//   break; 
//   }

//  if(ok)
//    continue;


// for (int kk=1; kk<fMT2tree->NEles ; kk++){

//   //  float deltaR = Util::GetDeltaR(fMT2tree->ele[kk].lv.Eta(), fMT2tree->jet[0].lv.Eta(), fMT2tree->ele[kk].lv.Phi(), fMT2tree->jet[0].lv.Phi());
//   //  if(deltaR > 0.5 && deltaR < 6 && /*fMT2tree->ele[kk].PassE0_EE == 1*/ fMT2tree->ele[kk].IDSelQCD == 1 && fMT2tree->ele[kk].MT < 40)
//   if((((fMT2tree->ele[0].Charge)+(fMT2tree->ele[kk].Charge) == 2) || ((fMT2tree->ele[0].Charge)+(fMT2tree->ele[kk].Charge) == -2)) /*&&  fMT2tree->ele[0].IDSelQCD == 1*/ &&  fMT2tree->ele[0].Iso04 > 0.15 && fMT2tree->ele[kk].MT < 30 && fMT2tree->ele[0].MT < 30 && fMT2tree->ele[kk].PassE0_EE == 1)
//    { 

       TString myQuantity = fMT2tree->GenLeptonAnalysisInterpretation(2,0,0, false);

      if(data == 1){
      
	MT2[7]->Fill(myQuantity, weight);//data
      
      }else{
	if(Sample.sname == "SUSY")
	  MT2[6]->Fill(myQuantity, weight);
	else
	  MT2[5]->Fill(myQuantity, weight);
      
	if(Sample.sname == "Top")
	  MT2[3]->Fill(myQuantity, weight);
	else
	  if(Sample.sname == "DY")	
	    MT2[2]->Fill(myQuantity, weight);
	  else
	    if(Sample.sname == "Wtolnu"){
// 	      float pt = fMT2tree->tau[tauIndex].lv.Pt();
// 	      weight *= 1.157 - 7.361E-3 * pt + 4.370E-5 * pt * pt - 1.188E-7*pt * pt * pt;
	      MT2[1]->Fill(myQuantity, weight);}
	    else
	      if(Sample.sname == "QCD")
		MT2[0]->Fill(myQuantity, weight);
	      else
		if(Sample.sname == "VV")
		  MT2[4]->Fill(myQuantity, weight);}

      // }
      // }//for electrons
 }//for events
 }//for samples

  for(int j = 0; j < (NumberOfSamples); j++){
    AddOverAndUnderFlow(MT2[j], true, true);
  }
   printYield();

  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j < (NumberOfSamples); j++){
    // MT2[j]->Rebin(3);
    if(j < (NumberOfSamples - 3))
      h_stack -> Add(MT2[j]);
  }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW+jets", "f");
  Legend1->AddEntry(MT2[6], "SMS", "l");
  Legend1->AddEntry(MT2[7], "data", "l");


  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + myfileName +"_Histos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  h_stack->Write();
  MT2[0]->Write();
  MT2[1]->Write();
  MT2[2]->Write();
  MT2[3]->Write();
  MT2[4]->Write();
  MT2[5]->Write();
  MT2[6]->Write();
  MT2[7]->Write();
  Legend1->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  printHisto(h_stack, MT2[7], MT2[5], MT2[6], Legend1 , myfileName, "hist", true, myfileName, "Events", 0, -10, 2, true);

  plotRatioStack(h_stack, MT2[5], MT2[7], MT2[6], true, false, myfileName, Legend1, myfileName, "Events", 0, -10, 2, true);


}

void MassPlotter::eeAnalysis(TString cuts, TString trigger, unsigned int nevents, TString action, TString variable, TString myfileName){

  TH1::SetDefaultSumw2();

  //  TH2D *h_PN_MLSP_MChi = new TH2D("h_PN_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
  //  h_PN_MLSP_MChi->Sumw2();

  //  TH2D *h_N_MLSP_MChi = new TH2D("h_N_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
  //  h_N_MLSP_MChi->Sumw2();


  TString  cnames[NumberOfSamples+1] = {"QCD", "W", "ZX", "Top", "WW", "Higgs", "MC", "SUSY", "data" };
  int      ccolor[NumberOfSamples+1] = { 401 , 417,  419,   855,  603,    kRed,  603,      1,    632 };

  static const int nbins = 5;
  double bins[nbins+1] = {40, 60, 90, 120, 150, 200};//MT2
 //  static const int nbins = 18;
 //  double bins[nbins+1] = {-500, -200, -150, -120, -100, -80, -65, -50, -40, 0, 40, 50, 65, 80, 100, 120, 150, 200, 500};//JZB
  //  static const int nbins = 8;
  //  double bins[nbins+1] = {0.0, 20.0, 40.0, 60.0, 90.0, 120.0, 150.0, 200.0, 500.0};

  TString varname = "MT2";
  for (int i=0; i <= (NumberOfSamples) ; i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i],"",nbins, bins);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[NumberOfSamples] -> SetMarkerStyle(20);
  MT2[NumberOfSamples] -> SetMarkerColor(kBlack);
  MT2[NumberOfSamples] -> SetLineColor(kBlack);
  
  MT2[NumberOfSamples-2] -> SetFillStyle(3004);
  MT2[NumberOfSamples-2] -> SetFillColor(kBlack);
  
  MT2[NumberOfSamples-1] -> SetMarkerStyle(20);
  MT2[NumberOfSamples-1] -> SetMarkerColor(kRed-4);
  MT2[NumberOfSamples-1] -> SetLineColor(kRed-4);
  MT2[NumberOfSamples-1] -> SetLineWidth(3);


  for(unsigned int ii = 0; ii < fSamples.size(); ii++)
    {
      TString myCuts = cuts;
      
      int data = 0;
      sample Sample = fSamples[ii];

      if(Sample.type == "data")
	{
	  data = 1;
	  myCuts += " && " + trigger;
	}


  cout<<" trigger "<< trigger<<endl;
  cout<<" cuts "<< myCuts<<endl;

      fMT2tree = new MT2tree();
      Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      std::cout << setfill('=') << std::setw(70) << "" << std::endl;
      cout << "looping over :     " <<endl;	
      cout << "   Name:           " << Sample.name << endl;
      cout << "   File:           " << (Sample.file)->GetName() << endl;
      cout << "   Events:         " << Sample.nevents  << endl;
      cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
      cout << "   Xsection:       " << Sample.xsection << endl;
      cout << "   kfactor:        " << Sample.kfact << endl;
      cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
      cout << "   Weight:         " << Weight <<endl;
      std::cout << setfill('-') << std::setw(70) << "" << std::endl;
      
      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

      unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();

      float Weight=0.0;



      
      if(fPUReweight) 
	Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
       else
	 Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

      for ( unsigned int jentry=0 ; jentry<min(nentries , nevents);jentry++ )
	{

// 	  bool RejMu_ee = false;
// 	  for(int i=0; i<fMT2tree->NMuons; i++){ 
// 	    if(fMT2tree->muo[i].RejMu_EE)
// 	      {
// 		RejMu_ee = true;
// 		break;
// 	      }
// 	  }
// 	  //	  cout << "RejMu_ee" << RejMu_ee << endl; 


//  	  bool RejE_ee = false;
//  	  for (int i=0 ; i < fMT2tree->NEles ; i++){
//  	    if (i != fMT2tree->doubleEle[0].Ele0Ind && i != fMT2tree->doubleEle[0].Ele1Ind)
//  	      {
// 		if(fMT2tree->ele[i].RejE2_EE)// ( fMT2tree->ele[i].Iso04 < 0.5 && fMT2tree->ele[i].lv.Eta() < 1.479 )
//  		{		
//  		  RejE_ee = true;
//  		  break;
//  		}
//  	      }
//  	  }
//  	  //	  cout << "Reje_ee" << Reje_ee << endl; 

// 	  // 	  if (RejMu_ee == true)
// 	  // 	    continue;
//  	  if (RejE_ee == true)
//  	    continue;

//	  if (fMT2tree->NMuons > 0)
//          continue;                                                                                                                                 
//	  if (fMT2tree->NEles > 2)
//          continue;                                                                                                                                 


	  //Sample.tree->GetEntry(jentry);
	  Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

	  if (fVerbose>2 && jentry % 100000 == 0 )
	    {
	      fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	      fflush(stdout);
	    }

	  if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")))
	    {
	      Weight = Sample.lumi;
	      if(fPUReweight) 
		Weight /= Sample.PU_avg_weight;
	    }

	  float weight = Weight;
	  
	  if(data == 1)
	    weight = 1.0;
	  else
	    {

	      //	 if(fMT2tree->NBJetsCSVL == 0)// || fMT2tree->NBJetsCSVM == 0 || fMT2tree->NBJetsCSVT == 0)
	   weight *=fMT2tree->SFWeight.BTagCSV40eq0;
	   //	 if(fMT2tree->NBJetsCSVL == 1)// || fMT2tree->NBJetsCSVM == 1 || fMT2tree->NBJetsCSVT == 1)
	   //	   weight *=fMT2tree->SFWeight.BTagCSV40eq1;
	   //	 if(fMT2tree->NBJetsCSVL == 2)// || fMT2tree->NBJetsCSVM == 2 || fMT2tree->NBJetsCSVT == 2)
	   //	   weight *=fMT2tree->SFWeight.BTagCSV40eq2;
	   //	 if(fMT2tree->NBJetsCSVL >= 3)// || fMT2tree->NBJetsCSVM >= 3 || fMT2tree->NBJetsCSVT >= 3)
	   //	   weight *=fMT2tree->SFWeight.BTagCSV40ge3;

	
	 //	      weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	      weight *= fMT2tree->pileUp.Weight;
	      weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      

	    }




	  //    double myQuantity = fMT2tree->CalcMT2(0, false, fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv, tau_mutau_down, MET_down); 
	  //      double myQuantity = fMT2tree->muTau[0].GetMT2();
	  //         float myQuantity = fMT2tree->doubleEle[0].MT2;
	  //    double myQuantity = fMT2tree->muTau[0].GetMT2();
	  //      myChannelCuts.push_back("ele[0].lv.Pt() > 30");
	  //      myChannelCuts.push_back("NEles == 2");
	  //      myChannelCuts.push_back("ele[0].PassE0_EE == 1");
	  //      myChannelCuts.push_back("ele[1].PassE0_EE == 1 && (ele[0].Charge + ele[1].Charge == 2 || ele[0].Charge + ele[1].Charge == -2)");
	  //      if(!((fMT2tree->NEles == 2) && /*(fMT2tree->ele[0].lv.Pt() > 30) &&*/ (fMT2tree->ele[0].PassE0_EE == 1) &&/* (fMT2tree->ele[0].Charge+fMT2tree->ele[1].Charge == 2 || fMT2tree->ele[0].Charge+fMT2tree->ele[1].Charge == -2) && (fMT2tree->ele[1].PassE0_EE == 1*/(fMT2tree->ele[1].QCDSyst03E0_EE==1)))
	  //   continue;
	  //   int Neles= -10;
	  //   int EleTightIndex= -10;
	  //   int MuTightIndex = -10;
	  //---------------------jet-----------------
	  // bool ok =false;      
	  //   //float deltaR1 = Util::GetDeltaR(fMT2tree->muo[ll].lv.Eta(), fMT2tree->jet[0].lv.Eta(), fMT2tree->muo[ll].lv.Phi(), fMT2tree->jet[0].lv.Phi());
	  // if(fMT2tree->muo[ll].PassMu0_EleMu ==1)
	  //     ok =true;
	  //     break;	
	  // }
	  // if(ok)
	  //     continue;
	  
	  // for(int qq=1; qq < fMT2tree->NEles; qq++)
	  //   {
	  //     //float deltaR = Util::GetDeltaR(fMT2tree->ele[kk].lv.Eta(), fMT2tree->jet[0].lv.Eta(), fMT2tree->ele[kk].lv.Phi(), fMT2tree->jet[0].lv.Phi());
	  
	  // if(fMT2tree->ele[qq].PassE0_EE == 1) //&& (fMT2tree->ele[kk].MT < 70) 
	  //   //  continue;
	  
	  //   {
	  //        for(int jj=1; jj < fMT2tree->NEles; jj++)
	  //        {
	  //          if(jj != qq && fMT2tree->ele[jj].PassE0_EE == 1) 
	  //          ok = true;
	  //          break;
	  //        }
	  //   }
	  
	  // if(ok)
	  //   break; 
	  //   }
	  
	  //  if(ok)
	  //    continue;
 
	  
	  //for (int kk=1; kk<fMT2tree->NEles ; kk++){
	  
	  //  float deltaR = Util::GetDeltaR(fMT2tree->ele[kk].lv.Eta(), fMT2tree->jet[0].lv.Eta(), fMT2tree->ele[kk].lv.Phi(), fMT2tree->jet[0].lv.Phi());
	  //  if(deltaR > 0.5 && deltaR < 6 && /*fMT2tree->ele[kk].PassE0_EE == 1*/ fMT2tree->ele[kk].IDSelQCD == 1 && fMT2tree->ele[kk].MT < 40)
	  
	  //  if((((fMT2tree->ele[0].Charge)+(fMT2tree->ele[kk].Charge) == 2) || ((fMT2tree->ele[0].Charge)+(fMT2tree->ele[kk].Charge) == -2)) &&  fMT2tree->ele[0].IDSelQCD == 1 &&  fMT2tree->ele[0].Iso04 > 0.15 && fMT2tree->ele[kk].MT < 30 && fMT2tree->ele[0].MT < 50)
	  //  cout << "IDSelQD...." << kk <<"..." <<fMT2tree->ele[kk].IDSelQCD << endl;
	  //  if(fMT2tree->ele[kk].IDSelQCD == 1)   { 
	  //    {
	  //---------------------singleEle------------------------
	  //   	    {ok = true; break;}
	  //   	  else 
	  //   	    {ok=false; continue;}
	  //  	}

	  // 	//  	cout << "ok:-----" << ok << endl;
	  //   	if (ok){
	  //for (int kk=1; kk < fMT2tree->NEles; kk++)
	  //{
	  
	  //  if(fMT2tree->ele[kk].PassE0_EE == 1){
	  //  int PairCharge = (fMT2tree->ele[0].Charge)+(fMT2tree->ele[kk].Charge) ;
	  //cout << "PairCharge Before:-----" << PairCharge << endl;
	  
	  //   if((PairCharge  == -2 || PairCharge == 2) && (fMT2tree->ele[kk].PassE0_EE == 1))
	  //{
	  //cout << "PairCharge After:-----" << PairCharge << endl;
	  //  float deltaR = Util::GetDeltaR(fMT2tree->ele[kk].lv.Eta(), fMT2tree->jet[0].lv.Eta(), fMT2tree->ele[kk].lv.Phi(), fMT2tree->jet[0].lv.Phi());
	  
	  //	if(deltaR < 0.5)
	  //	  continue;


	  //	if((fMT2tree->ele[0].MT) < 60 || (fMT2tree->ele[0].MT) > 90); 
	  //        {continue;}
	  
	  //        if (fMT2tree->NEles < 2)
	  // 	 continue; 
	  
	  //       int genLepIndexFromW = -2;;
	  
	  //        for(int i=0; i<fMT2tree->NGenLepts; ++i){
	  // 	  if(abs(fMT2tree->genlept[i].ID) !=11 &&  abs(fMT2tree->genlept[i].MID) !=24 ) 
	  //             continue;
	  // 	  //	  if( (!includeTau) && abs(genlept[i].MID) !=24       ) continue;
	  // 	  //	  if(   includeTau  && !((abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 ) || abs(genlept[i].MID)==24)) continue;
	  // 	  if(fMT2tree->genlept[i].lv.Pt() < 10 ) 
	  //             continue;
	  //  	  if(fabs(fMT2tree->genlept[i].lv.Eta())>2.5)
	  //             continue;
	  
	  // 	  //    deltaR = Util::GetDeltaR(fMT2tree->ele[0].lv.Eta(), fMT2tree->genlept[i].lv.Eta(), fMT2tree->ele[0].lv.Phi(), fMT2tree->genlept[i].lv.Phi());
	  
	  
	  // 	    genLepIndexFromW = i;
	  //             break;
	  //       }
	  
	  //        if (genLepIndexFromW == -2)
	  // 	 continue;
	  
	  //        //       float deltaR =0;
	  //        //       if (fMT2tree->ele[0].lv.Eta() < 5 && fMT2tree->genlept[genLepIndexFromW].lv.Eta() < 5)
	  //        //       if (fMT2tree->ele[0].lv.Eta() > 2.5)
	  //        //	 continue;
	  //        float deltaR =0; 
	  // deltaR = Util::GetDeltaR(fMT2tree->ele[1].lv.Eta(), fMT2tree->genlept[genLepIndexFromW].lv.Eta(), fMT2tree->ele[1].lv.Phi(), fMT2tree->genlept[genLepIndexFromW].lv.Phi());
	  
	  //	     cout << "deltaR.. " << deltaR << endl;
	  //       if (deltaR < 6)
	  //	 deltaR =10;
	  //	 {
	  //      continue;
	  
	  //   if (genLepIndeFromW == -2)
	  //       cout << "genLepIndexFromW......" << genLepIndexFromW << endl;
	  //jjjj
	  //       else 
	  //           continue;
	  //if (deltaR < 0.001)
	  //    deltaR = 0;
	  
	   //if (deltaR > 6)
	   //  continue; 
	  //      float myQuantity;
	  //      if (fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Pt() > fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Pt())
	  
	  //	myQuantity = fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].MT;
	  //      else 
	  
	  
	  //          double myQuantity =fMT2tree->misc.MinMetJetDPhiPt40;

	  double myQuantity=0; 

	  if(variable == "MT2")
	   myQuantity = fMT2tree->doubleEle[0].MT2;

	  if(variable == "MET")
	    myQuantity = fMT2tree->misc.MET;

	  if(variable == "JZB")
	    myQuantity = fMT2tree->eeJZBInDirect();

	  if(variable == "minMETJet")
	    myQuantity = fMT2tree->misc.MinMetJetDPhiPt40;

	  if(variable == "sumMT")
	    myQuantity = fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].MT + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].MT;



		 
	  if(data == 1){
      
	    MT2[NumberOfSamples]->Fill(myQuantity, weight);//data
	    
	  }else{
 	    if(Sample.sname == "SUSY")
 	      MT2[NumberOfSamples-1]->Fill(myQuantity, weight);
 	    else
 	      MT2[NumberOfSamples-2]->Fill(myQuantity, weight);
	    
	    if(Sample.sname == "Top")
	      MT2[3]->Fill(myQuantity, weight);
	    else
	      if(Sample.sname == "DY")	
		MT2[2]->Fill(myQuantity, weight);
	      else
		if(Sample.sname == "Wtolnu"){
		  MT2[1]->Fill(myQuantity, weight);}
		else
		  if(Sample.sname == "QCD")
		    MT2[0]->Fill(myQuantity, weight);
		  else
		    if(Sample.sname == "VV")
		      MT2[4]->Fill(myQuantity, weight);
		    else
		      if(Sample.sname == "Higgs")
			MT2[5]->Fill(myQuantity, weight);

	  }
	   
	    



	    



	  AddOverAndUnderFlow(MT2[0] , true , true);
	  AddOverAndUnderFlow(MT2[1] , true , true);
	  AddOverAndUnderFlow(MT2[2] , true , true);
	  AddOverAndUnderFlow(MT2[3] , true , true);
	  AddOverAndUnderFlow(MT2[4] , true , true);
	  AddOverAndUnderFlow(MT2[5] , true , true);
	  AddOverAndUnderFlow(MT2[6] , true , true);
	  AddOverAndUnderFlow(MT2[7] , true , true);
	  AddOverAndUnderFlow(MT2[8] , true , true);
	  //	  AddOverAndUnderFlow(MT2[9] , true , true);
// 	  if(Sample.type=="susy"){
	    
// 	    TString SUSYselection = TString::Format("(%s)", weight, myCuts.Data());
// 	    TString SUSYvariable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_PN_MLSP_MChi->GetName());
// 	    Sample.tree->Draw(SUSYvariable, SUSYselection, "goff");
	    
// 	    SUSYselection = TString::Format(" (%s)", myCuts.Data()); 
// 	    SUSYvariable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_N_MLSP_MChi->GetName());
// 	    Sample.tree->Draw(SUSYvariable, SUSYselection, "goff");
	    
// 	  }
	  


	  //       if(data == 1)
	  // 	{
	  // 	  MT2[7]->Fill(myQuantity, weight);//data
	  // 	}
	  //       else
	  // 	{
	  //           if(Sample.sname == "SUSY")
	  // 	    MT2[6]->Fill(myQuantity, weight);
	  // 	  else
	  // 	    MT2[5]->Fill(myQuantity, weight);//MC
	  
	  //              if(Sample.sname == "Top")
	  // 	       MT2[3]->Fill(myQuantity, weight);
	  // 	     else
	  // 	     if(Sample.sname == "DY")	
	  // 	       MT2[2]->Fill(myQuantity, weight);
	  //              else
	  // 	     if(Sample.sname == "Wtolnu")
	  // 	       MT2[1]->Fill(myQuantity, weight);
	  // 	     else
	  // 	     if(Sample.sname == "QCD")
	  // 	       MT2[0]->Fill(myQuantity, weight);
	  // 	     else
	  // 	     if(Sample.sname == "VV")
	  // 	       MT2[4]->Fill(myQuantity, weight);
	  //       }
	  
	  
	  //	    break; 
	  //        }
	  //    else 
	  //     continue;
	  
	  //    }
	  //}//for(electrons)

	  
	}//for(events)
      
      
      double e1 =0;
      for(int m = 0; m <= (NumberOfSamples); m++){
	
	e1 = MT2[m]->Integral();
	std::cout << setfill('#') << std::setw(50) << "" << std::endl;
	cout << "sample:" << cnames[m] << endl;
	cout << "all_content:..."<< e1 << endl;   
	
      } 
      
    }  //for(samples)

	  if (action == "topValid")
	    {

	      MT2[NumberOfSamples]->Add(MT2[NumberOfSamples-2], -1);
	      MT2[NumberOfSamples]->Add(MT2[3], +1);
	      MT2[2]->Scale(0);
	      MT2[1]->Scale(0);
	      MT2[0]->Scale(0);
	      MT2[4]->Scale(0);
	      MT2[5]->Scale(0);


	    }

  //  for(int j = 0; j <= 8; j++){
  //    AddOverAndUnderFlow(MT2[j] , true , true);
    //    AddOverAndUnderFlow(MT2[j]);
    //AddOverAndUnderFlow(MT2[j], true, true);
  //  }
   printYield();

  
  double e ;
  double f ;
  
  for(int kk = 0; kk <= NumberOfSamples; kk++){
    
    e = 0;
    f = 0;
    
    for(int mm = 0; mm < nbins+1; mm++){
      
      
      e += MT2[kk]->GetBinContent(mm);
      f += MT2[kk]->GetBinError(mm)*MT2[kk]->GetBinError(mm);
      
    } 
    double g = sqrt(f);
    
    std::cout << setfill('-') << std::setw(50) << "" << std::endl;
    cout << "sample:" << cnames[kk] << endl;
    cout << "all_bin_content:..."<< e << endl;   
    cout << "all_bin_error:..."<< g << endl;   
  }

    std::cout << setfill('#') << std::setw(50) << "" << std::endl;

  for(int ll = 0; ll <= NumberOfSamples ; ll++){
    
    double err[NumberOfSamples+1];
    
    std::cout << setfill('-') << std::setw(50) << "" << std::endl;
    cout << "sample:" << cnames[ll] << endl;
    cout << "all_bin_content:..."<< MT2[ll]-> IntegralAndError(0, MT2[ll]->GetNbinsX()+1, err[ll]) << endl;   
    cout << "all_bin_error:..."<< err[ll] << endl;   
  }


  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j <= (NumberOfSamples); j++){
    // MT2[j]->Rebin(3);
    if(j <= (NumberOfSamples - 3))
      h_stack -> Add(MT2[j]);
  }
  
  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W", "f");
  Legend1->AddEntry(MT2[2], "ZX", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW", "f");
  Legend1->AddEntry(MT2[5], "Higgs", "f"); 
  Legend1->AddEntry(MT2[NumberOfSamples-1], "SMS", "l");
  Legend1->AddEntry(MT2[NumberOfSamples], "data", "l");

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + myfileName +"_Histos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  h_stack->Write();
  
  for(int j = 0; j <= (NumberOfSamples); j++)
    MT2[j]->Write();
  Legend1->Write();
  
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;
  
  
//   TH1F *h_PN_Bkg  = (TH1F*)MT2[NumberOfSamples-2]->Clone();
//   h_PN_Bkg->SetName("h_PN_Bkg");
//   h_PN_Bkg->Rebin(h_PN_Bkg->GetNbinsX());
  
//   TH1F *h_PN_Data  = (TH1F*)MT2[NumberOfSamples]->Clone();
//   h_PN_Data->SetName("h_PN_Data");
//   h_PN_Data->Rebin(h_PN_Data->GetNbinsX());

//   h_SMSEvents->Rebin2D(4, 4);
//   h_PN_MLSP_MChi->Divide(h_SMSEvents);
  
//   TH2* hXsec = (TH2*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");
//   h_PN_MLSP_MChi->Multiply(hXsec);
  
//   TString fileName1 = fOutputDir;
//   if(!fileName1.EndsWith("/")) fileName1 += "/";
//   Util::MakeOutputDir(fileName1);
//   fileName1 = fileName1 + "countingForExclusion_" + myfileName +"_Histos.root";
//   TFile *savefile1 = new TFile(fileName1.Data(), "RECREATE");
//   savefile1 ->cd();
  
//   h_stack->SetName("h_stack");
  
//   h_stack->Write();
  //  h_PN_Bkg->Write();
  //  h_PN_Data->Write();
  //  h_PN_MLSP_MChi->Write();
  //  h_N_MLSP_MChi->Write();
  
  //  savefile1->Close();	
  //  std::cout << "Saved histograms in " << savefile1->GetName() << std::endl;
	  
  printHisto(h_stack, MT2[NumberOfSamples], MT2[NumberOfSamples-2], MT2[NumberOfSamples-1], Legend1 , myfileName, "hist", true, myfileName, "Events", -10, 0, -10, true);
  
  plotRatioStack(h_stack, MT2[NumberOfSamples-2], MT2[NumberOfSamples], MT2[NumberOfSamples-1], true, false, myfileName, Legend1, myfileName, "Events", -10, 0, -10, true, "C");
  


  //   printHisto(h_stack, MT2[7], MT2[5], MT2[6], Legend1 , myfileName, "hist", true, myfileName, "Events", -10, 0, -10, true);                        
  //   plotRatioStack(h_stack, MT2[5], MT2[7], MT2[6], true, false, myfileName, Legend1, myfileName, "Events", -10, 0, -10, true);    


   TCanvas *MyC = new TCanvas("MyC", "MyC");
    MyC->Divide(2,1);
   MyC->cd(1);
   TGraphAsymmErrors* sig1 = plotSig(MT2[NumberOfSamples-1], MT2[NumberOfSamples-2], "MT2", "Lower Cut", /*type = */ 1, /*sys = */ 0.1);
   sig1->Draw("ACP");

   MyC->cd(2);
   TGraphAsymmErrors* sig2 = plotSig(MT2[NumberOfSamples-1], MT2[NumberOfSamples-2], "MT2", "Upper Cut", /*type = */ 1, /*sys = */ 0.1);
   sig2->Draw("ACP");
//-----------------added---------------------^

  //  for(int j = 0; j < NumberOfSamples; j++){
  //     AddOverAndUnderFlow(MT2[j]);
//   }
//   printYield();

// THStack* h_stack     = new THStack(varname, "");
//  double xbin[7] = {0.0, 35.0, 70.0, 110.0, 150.0, 250.0, 2000};
//  TH1F *hnew[NumberOfSamples];
//   for(int j = 0; j < NumberOfSamples; j++){
//     hnew[j] = (TH1F*)MT2[j]->Rebin(6,"hnew",xbin);
//     TH1F* mt2 = (TH1F*)MT2[j]->Clone();
//     mt2->SetName("mt2");
//     if(j < (NumberOfSamples - 2))
//       //h_stack  -> Add(MT2[j]);
//       h_stack  -> Add(hnew[j]);
//     delete mt2;
//   }  

//    for(int j = 0; j < NumberOfSamples; j++){
//      //          AddOverAndUnderFlow(MT2[j], true, true);
//           AddOverAndUnderFlow(MT2[j]);
// //     //    MT2[j]->Multiply(pileup_data_up_histo);
//    }
//    //    printYield();

//    THStack* h_stack = new THStack(varname, "");
//    for(int j = 0; j < NumberOfSamples; j++){
//      // MT2[j]->Rebin(3);
//       if(j < (NumberOfSamples - 3))
//        h_stack -> Add(MT2[j]);
//    }

// //      double a ;
// //      double b ;


// //     for(int jj = 0; jj < NumberOfSamples; jj++){
// //     std::cout << setfill('#') << std::setw(70) << "" << std::endl;
// //      cout << "sample:" << cnames[jj] << endl;

// //       a = 0;
// //       b = 0;


// //    for(int k = 0; k < 35; k++){
// //     std::cout << setfill('=') << std::setw(70) << "" << std::endl;
// //     //    cout << "bin" << k <<":" << bins[k]<<"-"<<bins[k+1] <<"..."<<"its content:" <<MT2[jj]->GetBinContent(k)<< "its error: " << MT2[jj]->GetBinError(k)<<endl;


// //      a += MT2[jj]->GetBinContent(k);

// //      b += MT2[jj]->GetBinError(k)*MT2[jj]->GetBinError(k);

// //  } 
// //    double c =0;
// //    c = sqrt(b);
// //    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
// //    cout << "sample:" << cnames[jj] << endl;
// //    cout << "all_bin_content:..."<< a << endl;   
// //    cout << "all_bin_error:..."<< c << endl;   
// //    }


//       double e ;
//       double f ;


//     for(int kk = 0; kk < NumberOfSamples; kk++){

//       e = 0;
//       f = 0;

//    for(int mm = 0; mm < nbins; mm++){


//      e += MT2[kk]->GetBinContent(mm);

//      f += MT2[kk]->GetBinError(mm)*MT2[kk]->GetBinError(mm);

//  } 
//    double g = sqrt(f);

//    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
//    cout << "sample:" << cnames[kk] << endl;
//    cout << "all_bin_content:..."<< e << endl;   
//    cout << "all_bin_error:..."<< g << endl;   
//    }



//    //   static const int nbins = 8;//11;                                                                                                            
//    //   double bins[nbins+1] = {0.0,50.0,60.0,70.0,80.0,100.0,120.0,200.0,1000.0};      //MT2 

//   TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
//   Legend1->AddEntry(MT2[0], "QCD", "f");
//   Legend1->AddEntry(MT2[1], "W+jets", "f");
//   Legend1->AddEntry(MT2[2], "Z+jets", "f");
//   Legend1->AddEntry(MT2[3], "Top", "f");
//   Legend1->AddEntry(MT2[4], "WWto2L2Nu", "f");
//   Legend1->AddEntry(MT2[6], "Susy", "l");
//   Legend1->AddEntry(MT2[7], "data", "l");


//   TString fileName2 = fOutputDir;
//   if(!fileName2.EndsWith("/")) fileName2 += "/";
//   Util::MakeOutputDir(fileName2);
//   fileName2 = fileName2 + myfileName +"_eeA"+".root";
//   TFile *savefile2 = new TFile(fileName2.Data(), "RECREATE");
//   savefile2 ->cd();
//   h_stack->Write();
//   MT2[0]->Write();
//   MT2[1]->Write();
//   MT2[2]->Write();
//   MT2[3]->Write();
//   MT2[4]->Write();
//   MT2[5]->Write();
//   MT2[6]->Write();
//   MT2[7]->Write();
//   Legend1->Write();
//   savefile2->Close();
//   std::cout << "Saved histograms in " << savefile2->GetName() << std::endl;

//   cout<<" trigger "<<trigger<<endl;
//   cout<<" cuts "<<cuts<<endl;

//   printHisto(h_stack, MT2[7], MT2[5], MT2[6], Legend1 , myfileName, "hist", true, myfileName, "Events", -10, 0, -10, true);

//   plotRatioStack(h_stack, MT2[5], MT2[7], MT2[6], true, false, myfileName, Legend1, myfileName, "Events", -10, 0, -10, true);


}


void MassPlotter::eeAnalysisTESpUsys(TString cuts, TString trigger, unsigned int nevents, TString myfileName, TString giveStatus, TString sampleName ){


  TH1::SetDefaultSumw2();


  TH2D *h_PN_MLSP_MChi = new TH2D("h_PN_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
  h_PN_MLSP_MChi->Sumw2();

  TH2D *h_N_MLSP_MChi = new TH2D("h_N_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
  h_N_MLSP_MChi->Sumw2();


  TFile* pileup_data = new TFile("pileupData.root","READ");
  
  TH1F* pileup_data_histo = (TH1F*) pileup_data->Get("pileup");

  TFile* pileup_data_up = new TFile("pileupDataUp.root","READ");
  
  TH1F* pileup_data_up_histo = (TH1F*) pileup_data_up->Get("pileup");
 
  TFile* pileup_data_down = new TFile("pileupDataDown.root","READ");
  
  TH1F* pileup_data_down_histo = (TH1F*) pileup_data_down->Get("pileup");

  
  pileup_data_histo->Scale(1.0/pileup_data_histo->Integral());
  pileup_data_up_histo->Scale(1.0/pileup_data_up_histo->Integral());
  pileup_data_down_histo->Scale(1.0/pileup_data_down_histo->Integral());

  pileup_data_up_histo->Divide(pileup_data_histo);
  pileup_data_down_histo->Divide(pileup_data_histo);

//   for(int i = 0; i < pileup_data_up_histo->GetNbinsX();i++){
//     pileup_data_up_histo->SetBinContent(i+1,  pileup_data_up_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
//   }

//   for(int i = 0; i < pileup_data_down_histo->GetNbinsX();i++){
//     pileup_data_down_histo->SetBinContent(i+1,  pileup_data_down_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
//   }


  //-------------------added-------------------
  TString  cnames[NumberOfSamples+1] = {"QCD",     "W",    "ZX",  "Top",       "WW",  "Higgs",   "MC", "susy","data"};
  int      ccolor[NumberOfSamples+1] = { 401,      417,     419,   855,         603,     kRed,    603,      1, 632};
  //-------------------added-------------------


  //  TString  cnames[NumberOfSamples+1] = { "QCD", "Wjets", "Zjets", "Top", "WWjets",   "MC", "susy","data" };
  //  int      ccolor[NumberOfSamples+1] = { 401  ,     417,     419,   855,      603,    603,      1,  632  };

//  TString  cnames[11] = {"DYToLL_M10To50", "DYToLL_M50", "DY1ToLL_M50", "DY2ToLL_M50", "DY3ToLL_M50", "DY4ToLL_M50", "ZZJetsTo2L2Nu","ZZJetsTo2L2Q", "ZZJetsTo4L", "WZJetsTo3LNu", "WZJetsTo2L2Q"};
//   int      ccolor[11] = { 401, 417, 419, 855, 603, 650, 670, 840, 230, 584, 632};


//  static const int nbins = NumberOfBins;
  //  double bins[nbins+1] = {0.0};
  //  bins[nbins+1] = xbin;
  //{-2000.0, 0.0, 20.0, 40.0, 60.0, 90.0, 120.0, 150.0, 200.0, 500.0 , 2000.0};

  static const int nbins = 10;
  double bins[nbins+1] = {-2000.0, 0.0, 20.0, 40.0, 60.0, 90.0, 120.0, 150.0, 200.0, 500.0 , 2000.0};
  //  static const int nbins = 12;
  //  double bins[nbins+1] = {0.0, 5, 10.0, 15, 20, 25, 30, 35, 40, 45 , 50, 55, 60};
  TString varname = "MT2";
  for (int i=0; i <= NumberOfSamples ; i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i],"",nbins, bins);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[8] -> SetMarkerStyle(20);
  MT2[8] -> SetMarkerColor(kBlack);
  MT2[8] -> SetLineColor(kBlack);
  
  MT2[6] -> SetFillStyle(3004);
  MT2[6] -> SetFillColor(kBlack);

  MT2[7] -> SetMarkerStyle(20);
  MT2[7] -> SetMarkerColor(kRed-4);
  MT2[7] -> SetLineColor(kRed-4);
  MT2[7] -> SetLineWidth(3);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++)
    {


    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];

    if(Sample.sname != "SUSY" || Sample.sname != "DY")
      continue;


    //    if (sampleName != "")

    //      {
    //	if(sampleName != Sample.sname)
    //          continue;
    //      }

    //    else
                                 {

    if(Sample.type == "data")
      {
      data = 1;
      myCuts += " && " + trigger;
      }

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();


    if (Sample.name == "QCD-Pt-15-20-MuEnriched")
      fPUReweight = false;
    else
      fPUReweight = true;

    float Weight=0.0;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    for ( unsigned int jentry=0 ; jentry<min(nentries , nevents);jentry++ )
      {

      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 )
	{
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	}

      if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")))
	{
	Weight = Sample.lumi;
	if(fPUReweight) 
        Weight /= Sample.PU_avg_weight;
	}


   TString status = giveStatus ;
   //   TString mutau_pu_up, mutau_pu_down, mutau_tes_up, mutau_tes_down, ditau_bin1_pu_up, ditau_bin1_pu_down, ditau_bin1_tes_up, ditau_bin1_tes_down, ditau_bin2_pu_up, ditau_bin2_pu_down ,ditau_bin2_tes_up, ditau_bin2_tes_down, mutau_nominal, ditau_bin1_nominal, ditau_bin2_nominal; 

 
   double myQuantity=0; 

 
   float weight = Weight;


  if(status == "ditau_bin1_tes_up" || status == "ditau_bin1_tes_down" || status == "ditau_bin2_tes_up" || status == "ditau_bin2_tes_down" || status == "mutau_tes_up" || status == "mutau_tes_down"  || status == "etau_tes_up" || status == "etau_tes_down")
       {
      double px_all_up = 0;
      double py_all_up = 0;
      double pz_all_up = 0;
      double E_all_up  = 0;

   for(int l=0; l <fMT2tree->NTaus; l++)
     {
       px_all_up += 0.03*fMT2tree->tau[l].lv.Px();
       py_all_up += 0.03*fMT2tree->tau[l].lv.Py();
       pz_all_up += 0.03*fMT2tree->tau[l].lv.Pz();
       E_all_up  += 0.03*fMT2tree->tau[l].lv.E();
     }

   TLorentzVector tau_delta_up    (px_all_up ,  py_all_up,  pz_all_up,  E_all_up);
   TLorentzVector tau_delta_down  (-px_all_up, -py_all_up, -pz_all_up, -E_all_up);


   TLorentzVector MET_up   = fMT2tree->pfmet[0] - tau_delta_up;
   TLorentzVector MET_down = fMT2tree->pfmet[0] - tau_delta_down;

   TLorentzVector tau_mutau_up  (1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Px(), 
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pz(),  
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_mutau_down (0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Px(),
                                  0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Py(),
	                          0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pz(),
                                  0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_etau_up  (1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Px(), 
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pz(),  
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_etau_down (0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Px(),
                                  0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Py(),
	                          0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pz(),
                                  0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.E());


   TLorentzVector tau0_ditau_up (1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Px(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pz(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau0_ditau_down (0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Px(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Py(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pz(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau1_ditau_up (1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Px(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pz(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.E());

   TLorentzVector tau1_ditau_down (0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Px(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Py(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pz(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.E());

   double tau_MT_mutau_up  = fMT2tree->GetMT( tau_mutau_up , 0., MET_up , 0.);
   double tau_MT_etau_up  = fMT2tree->GetMT( tau_etau_up , 0., MET_up , 0.);
   double tau0_MT_ditau_up = fMT2tree->GetMT( tau0_ditau_up , 0., MET_up , 0.); 
   double tau1_MT_ditau_up = fMT2tree->GetMT( tau1_ditau_up , 0., MET_up , 0.);                  

   double tau_MT_mutau_down  = fMT2tree->GetMT( tau_mutau_down , 0., MET_down , 0.);                            
   double tau_MT_etau_down  = fMT2tree->GetMT( tau_etau_down , 0., MET_down , 0.);                            
   double tau0_MT_ditau_down = fMT2tree->GetMT( tau0_ditau_down , 0., MET_down , 0.);
   double tau1_MT_ditau_down = fMT2tree->GetMT( tau1_ditau_down , 0., MET_down , 0.);                      
            
   double maxMT_up   = max(tau0_MT_ditau_up   , tau1_MT_ditau_up);
   double maxMT_down = max(tau0_MT_ditau_down , tau1_MT_ditau_down);
   double sumMT_up   = tau0_MT_ditau_up   + tau1_MT_ditau_up;
   double sumMT_down = tau0_MT_ditau_down + tau1_MT_ditau_down;

   double mt2_mutau_up = fMT2tree->CalcMT2(0, false, fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv, tau_mutau_up, MET_up);
   double mt2_mutau_down = fMT2tree->CalcMT2(0, false, fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv, tau_mutau_down, MET_down);
   double invMass_mutau_up  = (tau_mutau_up   + fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv).M(); 
   double invMass_mutau_down= (tau_mutau_down + fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv).M(); 

   double mt2_etau_up = fMT2tree->CalcMT2(0, false, fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv, tau_etau_up, MET_up);
   double mt2_etau_down = fMT2tree->CalcMT2(0, false, fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv, tau_etau_down, MET_down);
   double invMass_etau_up  = (tau_etau_up   + fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv).M(); 
   double invMass_etau_down= (tau_etau_down + fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv).M(); 

   double mt2_ditau_up   = fMT2tree->CalcMT2(0, false, tau0_ditau_up, tau1_ditau_up, MET_up);
   double mt2_ditau_down = fMT2tree->CalcMT2(0, false, tau0_ditau_down, tau1_ditau_down, MET_down);
   double invMass_ditau_up    = (tau0_ditau_up + tau1_ditau_up).M();
   double invMass_ditau_down  = (tau0_ditau_down + tau1_ditau_down).M();

   double  MinMetJetDPhi_up   = fMT2tree->MinMetJetDPhi(0,40,5.0,100); 
   double  MinMetJetDPhi_down = fMT2tree->MinMetJetDPhi(0,40,5.0,-100);  

   //-----------------ditau_bin1_tes_up--------------------------

   if(status == "ditau_bin1_tes_up")
     {      
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->pileUp.Weight;
	 }
       
       if (tau0_ditau_up.Pt() <= 45 || tau1_ditau_up.Pt() <= 45)
	 continue;
       if (invMass_ditau_up <= 15 || (invMass_ditau_up >= 55  && invMass_ditau_up <= 85  ))
       	 continue;
       if (MinMetJetDPhi_up <= 1)
	 continue;
       if (mt2_ditau_up <= 50)
       	 continue;

       myQuantity =mt2_ditau_up;
     }

 //-----------------ditau_bin1_tes_down--------------------------

   else if(status == "ditau_bin1_tes_down")
     {      
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->pileUp.Weight;
	 }
       
       if (tau0_ditau_down.Pt() <= 45 || tau1_ditau_down.Pt() <= 45)
	 continue;
       if (invMass_ditau_down <= 15 || (invMass_ditau_down >= 55  && invMass_ditau_down <= 85  ))
	 continue;
       if (MinMetJetDPhi_down <= 1)
	 continue;
       if (mt2_ditau_down <= 50)
	 continue;
             myQuantity =mt2_ditau_down;
      
     }

    

 //-----------------ditau_bin2_tes_up--------------------------

   else if(status == "ditau_bin2_tes_up")
     {      
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   weight *= fMT2tree->pileUp.Weight;
	 }
       
       if (tau0_ditau_up.Pt() <= 45 || tau1_ditau_up.Pt() <= 45)
	 continue;
       if (invMass_ditau_up <= 15 || (invMass_ditau_up >= 55  && invMass_ditau_up <= 85  ))
	 continue;
       if (MinMetJetDPhi_up <= 1)
	 continue;
       if(mt2_ditau_up >= 50 || mt2_ditau_up <= 40)
	 continue;
       if (sumMT_up <= 250)
	 continue;
       myQuantity =mt2_ditau_up;
     }

   //-----------------ditau_bin2_tes_down--------------------------

   else if(status == "ditau_bin2_tes_down")
     {      
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   weight *= fMT2tree->pileUp.Weight;
	 }

       if (tau0_ditau_down.Pt() <= 45 || tau1_ditau_down.Pt() <= 45)
	 continue;
       if (invMass_ditau_down <= 15 || (invMass_ditau_down >= 55  && invMass_ditau_down <= 85  ))
	 continue;
       if (MinMetJetDPhi_down <= 1)
	 continue;
       if(mt2_ditau_down >= 50 || mt2_ditau_down <= 40)
	 continue;
       if (sumMT_down <= 250)
	 continue;
       myQuantity =mt2_ditau_down;
     }



   //-----------------mutau_tes_up--------------------------

   
   else if(status == "mutau_tes_up")
     {
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   if(fPUReweight)
	     weight *= fMT2tree->pileUp.Weight;
	   weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

	   if(Sample.sname == "Wtolnu")
	     weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	 }
       
       if (MET_up.Pt() <= 30)
	 continue;
       if (tau_mutau_up.Pt() <= 20)
         continue;
       if (MinMetJetDPhi_up <= 1)
	 continue;
       if ( invMass_mutau_up <= 15 || (invMass_mutau_up >= 45 && invMass_mutau_up <= 75))
	 continue;
       if (mt2_mutau_up <= 90)
	 continue;
       if (tau_MT_mutau_up <= 200)
	 continue;
       myQuantity =mt2_mutau_up;
     }

   //-----------------mutau_tes_down--------------------------
   
   else if(status == "mutau_tes_down")
     {
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   
	   if(fPUReweight)
	     weight *= fMT2tree->pileUp.Weight;
	   weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

	   if(Sample.sname == "Wtolnu")
	     weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	}

       if (MET_down.Pt() <= 30)
	 continue;
       if (tau_mutau_down.Pt() <= 20)
	 continue;
       if (MinMetJetDPhi_down <= 1)
	 continue;
       if ( invMass_mutau_down <= 15 || (invMass_mutau_down >= 45 & invMass_mutau_down <= 75))
	 continue;
       if (mt2_mutau_down <= 90)
	 continue;
       if (tau_MT_mutau_down <= 200)
	 continue;
       myQuantity =mt2_mutau_down;
     }

   //-----------------etau_tes_up--------------------------

   else if(status == "etau_tes_up")
     {
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   
	   if(fPUReweight)
	     weight *= fMT2tree->pileUp.Weight;
	   weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	   
	   if(Sample.sname == "Wtolnu")
	     weight *= fMT2tree->eleTau[0].tauWjetsSF;
	 }
       
       if (MET_up.Pt() <= 30)
	 continue;
       if (MinMetJetDPhi_up <= 1)
	 continue;
       if ( invMass_etau_up <= 15 || (invMass_etau_up >= 45 && invMass_etau_up <= 75))
	 continue;
       if (tau_etau_up.Pt() <= 20)
         continue;
       if (mt2_etau_up <= 90)
         continue;
       if (tau_MT_etau_up <= 200)
	 continue;
       myQuantity =mt2_etau_up;
     }

   //-----------------etau_tes_down--------------------------

   else 
     //if(status == "etau_tes_down")
     {
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   
	   if(fPUReweight)
	     weight *= fMT2tree->pileUp.Weight;
	   weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	   
	   if(Sample.sname == "Wtolnu")
	     weight *= fMT2tree->eleTau[0].tauWjetsSF;
	 }
       
       if (MET_down.Pt() <= 30)
	 continue;
       if (MinMetJetDPhi_down <= 1)
	 continue;
       if ( invMass_etau_down <= 15 || (invMass_etau_down >= 45 & invMass_etau_down <= 75))
	 continue;
       if (tau_etau_down.Pt() <= 20)
	 continue;
       if (mt2_etau_down <= 90)
	 continue;
       if (tau_MT_etau_down <= 200)
	 continue;
       myQuantity =mt2_etau_up;

     }
      
       }

    //----------------------ditau_bin1_nominal----------------------------

    else if(status == "ditau_bin1_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    
	    weight *= fMT2tree->pileUp.Weight;
	  }
	

	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;
	if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->doubleTau[0].GetMT2() <= 50)
	  continue;

	//	if(fMT2tree->doubleTau[0].GetMT2() <= 90)
	//	  continue;

		myQuantity = fMT2tree->doubleTau[0].GetMT2();
	//       myQuantity = fMT2tree->pileUp.NVertices; 
     }

    //----------------------ditau_bin2_nominal----------------------------
    
    else if(status == "ditau_bin2_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    weight *= fMT2tree->pileUp.Weight;
	  }
	
	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;
	if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->doubleTau[0].GetMT2() >= 50 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	  continue;
	if((fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT) + (fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT) <= 250)     
	  continue;

	myQuantity = fMT2tree->doubleTau[0].GetMT2();

      }     



  
  else if (status == "etau_ees_up" || status == "etau_ees_down") {

    //-----------------ees-------------------
    
    double px_all_up_e = 0;
    double py_all_up_e = 0;
    double pz_all_up_e = 0;
    double E_all_up_e  = 0;

    for(int l=0; l <fMT2tree->NEles; l++)
      {
	if (fabs(fMT2tree->ele[l].lv.Eta()) < 1.479)
	  {     
	    px_all_up_e += 0.01*fMT2tree->ele[l].lv.Px();
	    py_all_up_e += 0.01*fMT2tree->ele[l].lv.Py();
	    pz_all_up_e += 0.01*fMT2tree->ele[l].lv.Pz();
	    E_all_up_e  += 0.01*fMT2tree->ele[l].lv.E();
	  }
	else
	  {
	    px_all_up_e += 0.025*fMT2tree->ele[l].lv.Px();
	    py_all_up_e += 0.025*fMT2tree->ele[l].lv.Py();
	    pz_all_up_e += 0.025*fMT2tree->ele[l].lv.Pz();
	    E_all_up_e  += 0.025*fMT2tree->ele[l].lv.E();
	  }
      }
    
    TLorentzVector ele_delta_up    (px_all_up_e ,  py_all_up_e,  pz_all_up_e,  E_all_up_e);
    TLorentzVector ele_delta_down  (-px_all_up_e, -py_all_up_e, -pz_all_up_e, -E_all_up_e);
    
    TLorentzVector MET_up_e   = fMT2tree->pfmet[0] - ele_delta_up;
    TLorentzVector MET_down_e = fMT2tree->pfmet[0] - ele_delta_down;
    
    TLorentzVector ele_etau_up  (0,0,0,0) ;
    TLorentzVector ele_etau_down(0,0,0,0) ;

    if (fabs(fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Eta()) < 1.479)
      {
        ele_etau_up.SetPxPyPzE(1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(), 
			       1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(), 
			       1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(), 
			       1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E()) ;
	
	ele_etau_down.SetPxPyPzE(0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(),
				 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(),
				 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(),
				 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());
      }

    else
      {
	ele_etau_up.SetPxPyPzE(1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(), 
			       1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(), 
			       1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(), 
			       1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());
	
	ele_etau_down.SetPxPyPzE(0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());
      }

    double mt2_etau_up_e = fMT2tree->CalcMT2(0, false, fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv, ele_etau_up, MET_up_e);
    double mt2_etau_down_e = fMT2tree->CalcMT2(0, false, fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv, ele_etau_down, MET_down_e);
    double invMass_etau_up_e = (ele_etau_up   + fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv).M(); 
    double invMass_etau_down_e = (ele_etau_down + fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv).M(); 
    double MinMetJetDPhi_up_e = fMT2tree->MinMetJetDPhi(0,40,5.0,200); 
    double MinMetJetDPhi_down_e = fMT2tree->MinMetJetDPhi(0,40,5.0,-200);  

    //-----------------etau_ees_up--------------------------
    
    if(status == "etau_ees_up")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    
	    if(fPUReweight)
	      weight *= fMT2tree->pileUp.Weight;
	    weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	    if(Sample.sname == "Wtolnu")
	      weight *= fMT2tree->eleTau[0].tauWjetsSF;
       }

     
	if(MET_up_e.Pt() <= 30)
	  continue;
	if(MinMetJetDPhi_up_e <= 1)
	  continue;
	if(invMass_etau_up_e <= 15 || (invMass_etau_up_e >= 45 && invMass_etau_up_e <= 75))
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20)     
	  continue;
	if(ele_etau_up.Pt() <= 25)
	  continue;
	if(mt2_etau_up_e <= 90)
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	  continue;
	myQuantity =mt2_etau_up_e;
      }

    //-----------------etau_ees_down--------------------------
    
    else if(status == "etau_ees_down")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    
	    if(fPUReweight)
	      weight *= fMT2tree->pileUp.Weight;
	    weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	    if(Sample.sname == "Wtolnu")
	      weight *= fMT2tree->eleTau[0].tauWjetsSF;
	  }
	

	if(MET_down_e.Pt() <= 30)
	  continue;
	if(MinMetJetDPhi_down_e <= 1)
	  continue;
	if(invMass_etau_down_e <= 15 || (invMass_etau_down_e >= 45 && invMass_etau_down_e <= 75))
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20)     
	  continue;
	if(ele_etau_down.Pt() <= 25)
	  continue;
	if(mt2_etau_down_e <= 90)
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	  continue;

	myQuantity =mt2_etau_down_e;
      }
  }



    //----------------------mutau_nominal----------------------------

    else if(status == "mutau_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    
	    if(fPUReweight)
	      weight *= fMT2tree->pileUp.Weight;
	    weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();
	    
	    if(Sample.sname == "Wtolnu")
	      weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	  }

	if (fMT2tree->misc.MET <= 30)
	  continue;
	if (fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->muTau[0].GetMT2() <= 90)
	  continue;
	if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	  continue;

	  myQuantity = fMT2tree->muTau[0].GetMT2();


      }     

    //----------------------etau_nominal----------------------------
    
    else if(status == "etau_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    
	    if(fPUReweight)
	      weight *= fMT2tree->pileUp.Weight;
	    weight *=  fMT2tree->eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	    if(Sample.sname == "Wtolnu")
	      weight *= fMT2tree->eleTau[0].tauWjetsSF;
	  }
	
	if (fMT2tree->misc.MET <= 30)
	  continue;
	if (fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pt() <= 25 )
	  continue;
	if (fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->eleTau[0].GetMT2() <= 90)
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	  continue;

	myQuantity = fMT2tree->eleTau[0].GetMT2();

      }     

   //----------------------ditau_bin1_pu_up----------------------------

else if(status == "ditau_bin1_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;
	}

	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;
	if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->doubleTau[0].GetMT2() <= 40)
	  continue;
	//	if(fMT2tree->doubleTau[0].GetMT2() <= 90)
	//	  continue;
	
	//	     myQuantity = fMT2tree->doubleTau[0].GetMT2();
       myQuantity = fMT2tree->pileUp.NVertices; 
     }


   //----------------------------ditau_bin1_pu_down---------------------

else if(status == "ditau_bin1_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{

		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;

	}


	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	if(fMT2tree->doubleTau[0].GetMT2() <= 40)
	  continue;


	//	   if(fMT2tree->doubleTau[0].GetMT2() <= 90)
	//   continue;
       myQuantity = fMT2tree->pileUp.NVertices; 
       //	   myQuantity = fMT2tree->doubleTau[0].GetMT2();
     }


   //----------------------ditau_bin2_pu_up----------------------------

else if(status == "ditau_bin2_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;

	}

	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if(fMT2tree->doubleTau[0].GetMT2() >= 90 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	     continue;
	   if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT <= 250)     
	     continue;

	     myQuantity = fMT2tree->doubleTau[0].GetMT2();

     }


   //----------------------------ditau_bin2_pu_down---------------------

else if(status == "ditau_bin2_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;


	}

	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if(fMT2tree->doubleTau[0].GetMT2() >= 90 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	     continue;
	   if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT <= 250)     
	     continue;

	   myQuantity = fMT2tree->doubleTau[0].GetMT2();
     }


   //----------------------mutau_pu_up----------------------------

else if(status == "mutau_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;
	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}

	if (fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Pt() <= 20 )
	  continue;

	   if ( fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree-> misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->muTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	     continue;


     myQuantity = fMT2tree->muTau[0].GetMT2();

     }


   //----------------------------mutau_pu_down---------------------

else if(status == "mutau_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;
	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}

	if (fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Pt() <= 20 )
	  continue;

	   if (fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->muTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	     continue;

	   myQuantity = fMT2tree->muTau[0].GetMT2();
    
     }


   //----------------------etau_pu_up----------------------------

else if(status == "etau_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;

	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;


	}

	if (fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pt() <= 25 )
	  continue;

	   if ( fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree-> misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->eleTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	     continue;


     myQuantity = fMT2tree->eleTau[0].GetMT2();

     }


   //----------------------------etau_pu_down---------------------

 else //(status == "etau_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;

	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;

	}


	if (fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pt() <= 25 )
	  continue;

	   if (fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->eleTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	     continue;

	   myQuantity = fMT2tree->eleTau[0].GetMT2();
     }
      
      

      if(data == 1){

  	MT2[8]->Fill(myQuantity, weight);//data
      
       }else{
	if(Sample.sname == "SUSY")
	  MT2[7]->Fill(myQuantity, weight);
	else

	  MT2[6]->Fill(myQuantity, weight);
      
	if(Sample.sname == "Top")
	  MT2[3]->Fill(myQuantity, weight);
	else
	  if(Sample.sname == "DY")	
	    MT2[2]->Fill(myQuantity, weight);
	  else
	    if(Sample.sname == "Wtolnu"){
	      MT2[1]->Fill(myQuantity, weight);}
	    else
	      if(Sample.sname == "QCD")
		MT2[0]->Fill(myQuantity, weight);
	      else
		if(Sample.sname == "VV")
		  MT2[4]->Fill(myQuantity, weight);
		else
		  if(Sample.sname == "Higgs")
		    MT2[5]->Fill(myQuantity, weight);
      }

	  AddOverAndUnderFlow(MT2[0] , true , true);
	  AddOverAndUnderFlow(MT2[1] , true , true);
	  AddOverAndUnderFlow(MT2[2] , true , true);
	  AddOverAndUnderFlow(MT2[3] , true , true);
	  AddOverAndUnderFlow(MT2[4] , true , true);
	  AddOverAndUnderFlow(MT2[5] , true , true);
	  AddOverAndUnderFlow(MT2[6] , true , true);
	  AddOverAndUnderFlow(MT2[7] , true , true);
	  AddOverAndUnderFlow(MT2[8] , true , true);

  



//  		  if(Sample.type=="susy"){

//  		    TString SUSYselection = TString::Format("(%s)", weight, myCuts.Data());
//  		    TString SUSYvariable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_PN_MLSP_MChi->GetName());
//  		    Sample.tree->Draw(SUSYvariable, SUSYselection, "goff");

//  		    SUSYselection = TString::Format(" (%s)", myCuts.Data()); 
//  		    SUSYvariable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_N_MLSP_MChi->GetName());
//  		    Sample.tree->Draw(SUSYvariable, SUSYselection, "goff");
		  
//  		  }

//               if(Sample.name == "DYToLL_M10To50")
// 		  MT2[0]->Fill(myQuantity, weight);
//               else 
// 	        if(Sample.name == "DYToLL_M50")
// 		  MT2[1]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY1ToLL_M50")	
// 		  MT2[2]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY2ToLL_M50")	
// 		  MT2[3]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY3ToLL_M50")	
// 		  MT2[4]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY4ToLL_M50")	
// 		  MT2[5]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo2L2Nu")	
// 		  MT2[6]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo2L2Q")	
// 		  MT2[7]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo4L")	
// 		  MT2[8]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "WZJetsTo3LNu")	
// 		  MT2[9]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "WZJetsTo2L2Q")	
// 		  MT2[10]->Fill(myQuantity, weight);

	
	

    
      }
//For(events)

    double e1 =0;
   for(int m = 0; m <= NumberOfSamples; m++){

   e1 = MT2[m]->Integral();
   std::cout << setfill('#') << std::setw(50) << "" << std::endl;
   cout << "sample:" << cnames[m] << endl;
   cout << "all_content:..."<< e1 << endl;   

   } 
   

				 }
    }
  //for(samples)



//  for(int j = 0; j <= NumberOfSamples; j++){
//    AddOverAndUnderFlow(MT2[j], true, true);
    //     AddOverAndUnderFlow(MT2[j], true, true);
     //          AddOverAndUnderFlow(MT2[j]);
//        MT2[j]->Multiply(pileup_data_up_histo);
//   }


    THStack* h_stack = new THStack(varname, "");
    for(int j = 0; j <= NumberOfSamples; j++){
      // MT2[j]->Rebin(3);
       if(j <= (NumberOfSamples - 3))
        h_stack -> Add(MT2[j]);
    }

   //      Double a ;
   //      double b ;


//     for(int jj = 0; jj < NumberOfSamples; jj++){
//     std::cout << setfill('#') << std::setw(70) << "" << std::endl;
//      cout << "sample:" << cnames[jj] << endl;

//       a = 0;
//       b = 0;


//    for(int k = 0; k < 35; k++){
//     std::cout << setfill('=') << std::setw(70) << "" << std::endl;
//     //    cout << "bin" << k <<":" << bins[k]<<"-"<<bins[k+1] <<"..."<<"its content:" <<MT2[jj]->GetBinContent(k)<< "its error: " << MT2[jj]->GetBinError(k)<<endl;


//      a += MT2[jj]->GetBinContent(k);

//      b += MT2[jj]->GetBinError(k)*MT2[jj]->GetBinError(k);

//  } 
//    double c =0;
//    c = sqrt(b);
//    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
//    cout << "sample:" << cnames[jj] << endl;
//    cout << "all_bin_content:..."<< a << endl;   
//    cout << "all_bin_error:..."<< c << endl;   
//    }

      double e ;
      double f ;


      for(int kk = 0; kk <= NumberOfSamples; kk++){

      e = 0;
      f = 0;

   for(int mm = 0; mm < nbins+1 ; mm++){


     e += MT2[kk]->GetBinContent(mm);

     f += MT2[kk]->GetBinError(mm)*MT2[kk]->GetBinError(mm);

 } 
   double g = sqrt(f);

   std::cout << setfill('-') << std::setw(50) << "" << std::endl;
   cout << "sample:" << cnames[kk] << endl;
   cout << "all_bin_content:..."<< e << endl;   
   cout << "all_bin_error:..."<< g << endl;   
   }


  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W", "f");
  Legend1->AddEntry(MT2[2], "ZX", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW", "f");
  Legend1->AddEntry(MT2[5], "Higgs", "f"); 
  Legend1->AddEntry(MT2[7], "SMS", "l");
  Legend1->AddEntry(MT2[8], "data", "l");


//   TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
//   Legend1->AddEntry(MT2[0], "QCD", "f");
//   Legend1->AddEntry(MT2[1], "W+jets", "f");
//   Legend1->AddEntry(MT2[2], "Z+jets", "f");
//   Legend1->AddEntry(MT2[3], "Top", "f");
//   Legend1->AddEntry(MT2[4], "WWto2L2Nu", "f");
//   Legend1->AddEntry(MT2[6], "Susy", "l");
//   Legend1->AddEntry(MT2[7], "data", "l");


  TString fileName2 = fOutputDir;
  if(!fileName2.EndsWith("/")) fileName2 += "/";
  Util::MakeOutputDir(fileName2);
  fileName2 = fileName2 + myfileName +"_ee"+".root";
  TFile *savefile2 = new TFile(fileName2.Data(), "RECREATE");
  savefile2 ->cd();
  h_stack->Write();
  MT2[0]->Write();
  MT2[1]->Write();
  MT2[2]->Write();
  MT2[3]->Write();
  MT2[4]->Write();
  MT2[5]->Write();
  MT2[6]->Write();
  MT2[7]->Write();
  MT2[8]->Write();
  Legend1->Write();
  savefile2->Close();
  std::cout << "Saved histograms in " << savefile2->GetName() << std::endl;

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;



//   TH1F *h_PN_Bkg  = (TH1F*)MT2[NumberOfSamples-2]->Clone();
//   h_PN_Bkg->SetName("h_PN_Bkg");
//   h_PN_Bkg->Rebin(h_PN_Bkg->GetNbinsX());
  
//   TH1F *h_PN_Data  = (TH1F*)MT2[NumberOfSamples]->Clone();
//   h_PN_Data->SetName("h_PN_Data");
//   h_PN_Data->Rebin(h_PN_Data->GetNbinsX());

//   h_SMSEvents->Rebin2D(4, 4);
//   h_PN_MLSP_MChi->Divide(h_SMSEvents);
  
//   TH2* hXsec = (TH2*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");
//   h_PN_MLSP_MChi->Multiply(hXsec);

//   TString fileName1 = fOutputDir;
//   if(!fileName1.EndsWith("/")) fileName1 += "/";
//   Util::MakeOutputDir(fileName1);
//   fileName1 = fileName1 + "countingForExclusion_" + myfileName +"_Histos.root";
//   TFile *savefile1 = new TFile(fileName1.Data(), "RECREATE");
//   savefile1 ->cd();
  
//   h_stack->SetName("h_stack");
  
//   h_stack->Write();
//   h_PN_Bkg->Write();
//   h_PN_Data->Write();
//   h_PN_MLSP_MChi->Write();
//   h_N_MLSP_MChi->Write();
//   savefile1->Close();
	
//   std::cout << "Saved histograms in " << savefile1->GetName() << std::endl;


  printHisto(h_stack, MT2[8], MT2[6], MT2[7], Legend1 , myfileName, "hist", true, myfileName, "Events", -10, 0, -10, true);

  plotRatioStack(h_stack, MT2[6], MT2[8], MT2[7], true, false, myfileName, Legend1, myfileName, "Events", -10, 0, -10, true);


}

void MassPlotter::eeAnalysisTESpUsys1(TString cuts, TString trigger, unsigned int nevents, TString myfileName, TString giveStatus, TString sampleName, TString variable, double *xbin, int nbins){


  TH1::SetDefaultSumw2();


  TH2D *h_PN_MLSP_MChi = new TH2D("h_PN_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
  h_PN_MLSP_MChi->Sumw2();

  TH2D *h_N_MLSP_MChi = new TH2D("h_N_MLSP_MChi", "", 125, 0, 2500, 125, 0, 2500);
  h_N_MLSP_MChi->Sumw2();


  TFile* pileup_data = new TFile("pileupData.root","READ");
  
  TH1F* pileup_data_histo = (TH1F*) pileup_data->Get("pileup");

  TFile* pileup_data_up = new TFile("pileupDataUp.root","READ");
  
  TH1F* pileup_data_up_histo = (TH1F*) pileup_data_up->Get("pileup");
 
  TFile* pileup_data_down = new TFile("pileupDataDown.root","READ");
  
  TH1F* pileup_data_down_histo = (TH1F*) pileup_data_down->Get("pileup");

  
  pileup_data_histo->Scale(1.0/pileup_data_histo->Integral());
  pileup_data_up_histo->Scale(1.0/pileup_data_up_histo->Integral());
  pileup_data_down_histo->Scale(1.0/pileup_data_down_histo->Integral());

  pileup_data_up_histo->Divide(pileup_data_histo);
  pileup_data_down_histo->Divide(pileup_data_histo);

//   for(int i = 0; i < pileup_data_up_histo->GetNbinsX();i++){
//     pileup_data_up_histo->SetBinContent(i+1,  pileup_data_up_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
//   }

//   for(int i = 0; i < pileup_data_down_histo->GetNbinsX();i++){
//     pileup_data_down_histo->SetBinContent(i+1,  pileup_data_down_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
//   }


  //-------------------added-------------------
  TString  cnames[NumberOfSamples+1] = {"QCD",     "W",    "ZX",  "Top",       "WW",  "Higgs",   "MC", "susy","data"};
  int      ccolor[NumberOfSamples+1] = { 401,      417,     419,   855,         603,     kRed,    603,      1, 632};
  //-------------------added-------------------


  //  TString  cnames[NumberOfSamples+1] = { "QCD", "Wjets", "Zjets", "Top", "WWjets",   "MC", "susy","data" };
  //  int      ccolor[NumberOfSamples+1] = { 401  ,     417,     419,   855,      603,    603,      1,  632  };

//  TString  cnames[11] = {"DYToLL_M10To50", "DYToLL_M50", "DY1ToLL_M50", "DY2ToLL_M50", "DY3ToLL_M50", "DY4ToLL_M50", "ZZJetsTo2L2Nu","ZZJetsTo2L2Q", "ZZJetsTo4L", "WZJetsTo3LNu", "WZJetsTo2L2Q"};
//   int      ccolor[11] = { 401, 417, 419, 855, 603, 650, 670, 840, 230, 584, 632};


//  static const int nbins = NumberOfBins;
  //  double bins[nbins+1] = {0.0};
  //  bins[nbins+1] = xbin;
  //{-2000.0, 0.0, 20.0, 40.0, 60.0, 90.0, 120.0, 150.0, 200.0, 500.0 , 2000.0};

  TString varname = "MT2";
  for (int i=0; i <= NumberOfSamples ; i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i],"",nbins, xbin);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[8] -> SetMarkerStyle(20);
  MT2[8] -> SetMarkerColor(kBlack);
  MT2[8] -> SetLineColor(kBlack);
  
  MT2[6] -> SetFillStyle(3004);
  MT2[6] -> SetFillColor(kBlack);

  MT2[7] -> SetMarkerStyle(20);
  MT2[7] -> SetMarkerColor(kRed-4);
  MT2[7] -> SetLineColor(kRed-4);
  MT2[7] -> SetLineWidth(3);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++)
    {


    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];

      if (sampleName != "")
	
	{
       if(Sample.sname != sampleName)
          continue;
	}

      else


    /////////////////////////
    //    if(Sample.sname != "DY")
    //      continue;
    //////////////////////////

    if(Sample.type == "data")
      {
      data = 1;
      myCuts += " && " + trigger;
      }

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();


    if (Sample.name == "QCD-Pt-15-20-MuEnriched")
      fPUReweight = false;
    else
      fPUReweight = true;

    float Weight=0.0;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    for ( unsigned int jentry=0 ; jentry<min(nentries , nevents);jentry++ )
      {

      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 )
	{
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	}

      if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")))
	{
	Weight = Sample.lumi;
	if(fPUReweight) 
        Weight /= Sample.PU_avg_weight;
	}


   TString status = giveStatus ;
   //   TString mutau_pu_up, mutau_pu_down, mutau_tes_up, mutau_tes_down, ditau_bin1_pu_up, ditau_bin1_pu_down, ditau_bin1_tes_up, ditau_bin1_tes_down, ditau_bin2_pu_up, ditau_bin2_pu_down ,ditau_bin2_tes_up, ditau_bin2_tes_down, mutau_nominal, ditau_bin1_nominal, ditau_bin2_nominal; 

 
	  double myQuantity=0; 


// 	  if(variable == "MT2")
// 	    myQuantity = mt2_ditau_up;       
// 	  if(variable == "MT2")
// 	    myQuantity = mt2_ditau_down;
//      	  if(variable == "MT2")
// 	    myQuantity = mt2_ditau_up;
//        	  if(variable == "MT2")
// 	    myQuantity = mt2_ditau_down;
//      	  if(variable == "MT2")
// 	    myQuantity = mt2_mutau_up;
//        	  if(variable == "MT2")
// 	    myQuantity = mt2_mutau_down;
//        	  if(variable == "MT2")
// 	    myQuantity = mt2_etau_up;
//        	  if(variable == "MT2")
// 	    myQuantity = mt2_etau_down;
// 	  if(variable == "MT2")
// 	    myQuantity = mt2_etau_up_e;
// 	  if(variable == "MT2")
// 	    myQuantity = mt2_etau_down_e;
// 	  if(variable == "MT2")
// 	    myQuantity = mt2_ee_up;
// 	  if(variable == "MT2")
// 	    myQuantity = mt2_ee_down;
// 	  if(variable == "MT2")
// 	    myQuantity = mt2_ee_up;
// 	  if(variable == "MT2")
// 	    myQuantity = mt2_ee_down;
// 	  if(variable == "MT2")
// 	    myQuantity = fMT2tree->doubleEle[0].MT2 ;
// 	  if(variable == "MT2")
// 	    myQuantity = fMT2tree->doubleTau[0].GetMT2();//
// 	  if(variable == "MT2")
// 	    myQuantity = fMT2tree->doubleTau[0].GetMT2();
// 	  if(variable == "MT2")
// 	    myQuantity = fMT2tree->muTau[0].GetMT2();
// 	  if(variable == "MT2")
// 	    myQuantity = fMT2tree->eleTau[0].GetMT2();

		    
		    

 
      float weight = Weight;


  if(status == "ditau_bin1_tes_up" || status == "ditau_bin1_tes_down" || status == "ditau_bin2_tes_up" || status == "ditau_bin2_tes_down" || status == "mutau_tes_up" || status == "mutau_tes_down"  || status == "etau_tes_up" || status == "etau_tes_down")
    {
      double px_all_up = 0;
      double py_all_up = 0;
      double pz_all_up = 0;
      double E_all_up  = 0;

   for(int l=0; l <fMT2tree->NTaus; l++)
     {
       px_all_up += 0.03*fMT2tree->tau[l].lv.Px();
       py_all_up += 0.03*fMT2tree->tau[l].lv.Py();
       pz_all_up += 0.03*fMT2tree->tau[l].lv.Pz();
       E_all_up  += 0.03*fMT2tree->tau[l].lv.E();
     }

   TLorentzVector tau_delta_up    (px_all_up ,  py_all_up,  pz_all_up,  E_all_up);
   TLorentzVector tau_delta_down  (-px_all_up, -py_all_up, -pz_all_up, -E_all_up);


   TLorentzVector MET_up   = fMT2tree->pfmet[0] - tau_delta_up;
   TLorentzVector MET_down = fMT2tree->pfmet[0] - tau_delta_down;

   TLorentzVector tau_mutau_up  (1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Px(), 
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pz(),  
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_mutau_down (0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Px(),
                                  0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Py(),
	                          0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pz(),
                                  0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_etau_up  (1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Px(), 
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pz(),  
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_etau_down (0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Px(),
                                  0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Py(),
	                          0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pz(),
                                  0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.E());


   TLorentzVector tau0_ditau_up (1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Px(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pz(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau0_ditau_down (0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Px(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Py(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pz(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau1_ditau_up (1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Px(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pz(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.E());

   TLorentzVector tau1_ditau_down (0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Px(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Py(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pz(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.E());

   double tau_MT_mutau_up  = fMT2tree->GetMT( tau_mutau_up , 0., MET_up , 0.);
   double tau_MT_etau_up  = fMT2tree->GetMT( tau_etau_up , 0., MET_up , 0.);
   double tau0_MT_ditau_up = fMT2tree->GetMT( tau0_ditau_up , 0., MET_up , 0.); 
   double tau1_MT_ditau_up = fMT2tree->GetMT( tau1_ditau_up , 0., MET_up , 0.);                  

   double tau_MT_mutau_down  = fMT2tree->GetMT( tau_mutau_down , 0., MET_down , 0.);                            
   double tau_MT_etau_down  = fMT2tree->GetMT( tau_etau_down , 0., MET_down , 0.);                            
   double tau0_MT_ditau_down = fMT2tree->GetMT( tau0_ditau_down , 0., MET_down , 0.);
   double tau1_MT_ditau_down = fMT2tree->GetMT( tau1_ditau_down , 0., MET_down , 0.);                      
            
   double maxMT_up   = max(tau0_MT_ditau_up   , tau1_MT_ditau_up);
   double maxMT_down = max(tau0_MT_ditau_down , tau1_MT_ditau_down);
   double sumMT_up   = tau0_MT_ditau_up   + tau1_MT_ditau_up;
   double sumMT_down = tau0_MT_ditau_down + tau1_MT_ditau_down;

   double mt2_mutau_up = fMT2tree->CalcMT2(0, false, fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv, tau_mutau_up, MET_up);
   double mt2_mutau_down = fMT2tree->CalcMT2(0, false, fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv, tau_mutau_down, MET_down);
   double invMass_mutau_up  = (tau_mutau_up   + fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv).M(); 
   double invMass_mutau_down= (tau_mutau_down + fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv).M(); 

   double mt2_etau_up = fMT2tree->CalcMT2(0, false, fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv, tau_etau_up, MET_up);
   double mt2_etau_down = fMT2tree->CalcMT2(0, false, fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv, tau_etau_down, MET_down);
   double invMass_etau_up  = (tau_etau_up   + fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv).M(); 
   double invMass_etau_down= (tau_etau_down + fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv).M(); 

   double mt2_ditau_up   = fMT2tree->CalcMT2(0, false, tau0_ditau_up, tau1_ditau_up, MET_up);
   double mt2_ditau_down = fMT2tree->CalcMT2(0, false, tau0_ditau_down, tau1_ditau_down, MET_down);
   double invMass_ditau_up    = (tau0_ditau_up + tau1_ditau_up).M();
   double invMass_ditau_down  = (tau0_ditau_down + tau1_ditau_down).M();

   double  MinMetJetDPhi_up   = fMT2tree->MinMetJetDPhi(0,40,5.0,100); 
   double  MinMetJetDPhi_down = fMT2tree->MinMetJetDPhi(0,40,5.0,-100);  

   //-----------------ditau_bin1_tes_up--------------------------

   if(status == "ditau_bin1_tes_up")
     {      
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->pileUp.Weight;
	 }
       
       if (tau0_ditau_up.Pt() <= 45 || tau1_ditau_up.Pt() <= 45)
	 continue;
       if (invMass_ditau_up <= 15 || (invMass_ditau_up >= 55  && invMass_ditau_up <= 85  ))
       	 continue;
       if (MinMetJetDPhi_up <= 1)
	 continue;
       if (mt2_ditau_up <= 90)
       	 continue;


     }

 //-----------------ditau_bin1_tes_down--------------------------

   else if(status == "ditau_bin1_tes_down")
     {      
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->pileUp.Weight;
	 }
       
       if (tau0_ditau_down.Pt() <= 45 || tau1_ditau_down.Pt() <= 45)
	 continue;
       if (invMass_ditau_down <= 15 || (invMass_ditau_down >= 55  && invMass_ditau_down <= 85  ))
	 continue;
       if (MinMetJetDPhi_down <= 1)
	 continue;
       if (mt2_ditau_down <= 90)
	 continue;
     }

 //-----------------ditau_bin2_tes_up--------------------------

   else if(status == "ditau_bin2_tes_up")
     {      
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   weight *= fMT2tree->pileUp.Weight;
	 }
       
       if (tau0_ditau_up.Pt() <= 45 || tau1_ditau_up.Pt() <= 45)
	 continue;
       if (invMass_ditau_up <= 15 || (invMass_ditau_up >= 55  && invMass_ditau_up <= 85  ))
	 continue;
       if (MinMetJetDPhi_up <= 1)
	 continue;
       if(mt2_ditau_up >= 90 || mt2_ditau_up <= 40)
	 continue;
       if (sumMT_up <= 250)
	 continue;
     }

   //-----------------ditau_bin2_tes_down--------------------------

   else if(status == "ditau_bin2_tes_down")
     {      
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   weight *= fMT2tree->pileUp.Weight;
	 }

       if (tau0_ditau_down.Pt() <= 45 || tau1_ditau_down.Pt() <= 45)
	 continue;
       if (invMass_ditau_down <= 15 || (invMass_ditau_down >= 55  && invMass_ditau_down <= 85  ))
	 continue;
       if (MinMetJetDPhi_down <= 1)
	 continue;
       if(mt2_ditau_down >= 90 || mt2_ditau_down <= 40)
	 continue;
       if (sumMT_down <= 250)
	 continue;
     }

   //-----------------mutau_tes_up--------------------------
   
   else if(status == "mutau_tes_up")
     {
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   if(fPUReweight)
	     weight *= fMT2tree->pileUp.Weight;
	   weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

	   if(Sample.sname == "Wtolnu")
	     weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	 }
       
       if (MET_up.Pt() <= 30)
	 continue;
       if (tau_mutau_up.Pt() <= 20)
         continue;
       if (MinMetJetDPhi_up <= 1)
	 continue;
       if ( invMass_mutau_up <= 15 || (invMass_mutau_up >= 45 && invMass_mutau_up <= 75))
	 continue;
       if (mt2_mutau_up <= 90)
	 continue;
       if (tau_MT_mutau_up <= 200)
	 continue;

     }

   //-----------------mutau_tes_down--------------------------
   
   else if(status == "mutau_tes_down")
     {
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   
	   if(fPUReweight)
	     weight *= fMT2tree->pileUp.Weight;
	   weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

	   if(Sample.sname == "Wtolnu")
	     weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	}

       if (MET_down.Pt() <= 30)
	 continue;
       if (tau_mutau_down.Pt() <= 20)
	 continue;
       if (MinMetJetDPhi_down <= 1)
	 continue;
       if ( invMass_mutau_down <= 15 || (invMass_mutau_down >= 45 & invMass_mutau_down <= 75))
	 continue;
       if (mt2_mutau_down <= 90)
	 continue;
       if (tau_MT_mutau_down <= 200)
	 continue;
     }

   //-----------------etau_tes_up--------------------------

   else if(status == "etau_tes_up")
     {
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   
	   if(fPUReweight)
	     weight *= fMT2tree->pileUp.Weight;
	   weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	   
	   if(Sample.sname == "Wtolnu")
	     weight *= fMT2tree->eleTau[0].tauWjetsSF;
	 }
       
       if (MET_up.Pt() <= 30)
	 continue;
       if (MinMetJetDPhi_up <= 1)
	 continue;
       if ( invMass_etau_up <= 15 || (invMass_etau_up >= 45 && invMass_etau_up <= 75))
	 continue;
       if (tau_etau_up.Pt() <= 20)
         continue;
       if (mt2_etau_up <= 90)
         continue;
       if (tau_MT_etau_up <= 200)
	 continue;
     }

   //-----------------etau_tes_down--------------------------

   else 
     //if(status == "etau_tes_down")
     {
       if(data == 1)
	 weight = 1.0;
       else
	 {
	   weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	   
	   if(fPUReweight)
	     weight *= fMT2tree->pileUp.Weight;
	   weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	   
	   if(Sample.sname == "Wtolnu")
	     weight *= fMT2tree->eleTau[0].tauWjetsSF;
	 }
       
       if (MET_down.Pt() <= 30)
	 continue;
       if (MinMetJetDPhi_down <= 1)
	 continue;
       if ( invMass_etau_down <= 15 || (invMass_etau_down >= 45 & invMass_etau_down <= 75))
	 continue;
       if (tau_etau_down.Pt() <= 20)
	 continue;
       if (mt2_etau_down <= 90)
	 continue;
       if (tau_MT_etau_down <= 200)
	 continue;


     }
    }
  
  else if (status == "etau_ees_up" || status == "etau_ees_down") {

    //-----------------ees-------------------
    
    double px_all_up_e = 0;
    double py_all_up_e = 0;
    double pz_all_up_e = 0;
    double E_all_up_e  = 0;

    for(int l=0; l <fMT2tree->NEles; l++)
      {
	if (fabs(fMT2tree->ele[l].lv.Eta()) < 1.479)
	  {     
	    px_all_up_e += 0.01*fMT2tree->ele[l].lv.Px();
	    py_all_up_e += 0.01*fMT2tree->ele[l].lv.Py();
	    pz_all_up_e += 0.01*fMT2tree->ele[l].lv.Pz();
	    E_all_up_e  += 0.01*fMT2tree->ele[l].lv.E();
	  }
	else
	  {
	    px_all_up_e += 0.025*fMT2tree->ele[l].lv.Px();
	    py_all_up_e += 0.025*fMT2tree->ele[l].lv.Py();
	    pz_all_up_e += 0.025*fMT2tree->ele[l].lv.Pz();
	    E_all_up_e  += 0.025*fMT2tree->ele[l].lv.E();
	  }
      }
    
    TLorentzVector ele_delta_up    (px_all_up_e ,  py_all_up_e,  pz_all_up_e,  E_all_up_e);
    TLorentzVector ele_delta_down  (-px_all_up_e, -py_all_up_e, -pz_all_up_e, -E_all_up_e);
    
    TLorentzVector MET_up_e   = fMT2tree->pfmet[0] - ele_delta_up;
    TLorentzVector MET_down_e = fMT2tree->pfmet[0] - ele_delta_down;
    
    TLorentzVector ele_etau_up  (0,0,0,0) ;
    TLorentzVector ele_etau_down(0,0,0,0) ;

    if (fabs(fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Eta()) < 1.479)
      {
        ele_etau_up.SetPxPyPzE(1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(), 
			       1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(), 
			       1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(), 
			       1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E()) ;
	
	ele_etau_down.SetPxPyPzE(0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(),
				 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(),
				 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(),
				 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());
      }

    else
      {
	ele_etau_up.SetPxPyPzE(1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(), 
			       1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(), 
			       1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(), 
			       1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());
	
	ele_etau_down.SetPxPyPzE(0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());
      }

    double mt2_etau_up_e = fMT2tree->CalcMT2(0, false, fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv, ele_etau_up, MET_up_e);
    double mt2_etau_down_e = fMT2tree->CalcMT2(0, false, fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv, ele_etau_down, MET_down_e);
    double invMass_etau_up_e = (ele_etau_up   + fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv).M(); 
    double invMass_etau_down_e = (ele_etau_down + fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv).M(); 
    double MinMetJetDPhi_up_e = fMT2tree->MinMetJetDPhi(0,40,5.0,200); 
    double MinMetJetDPhi_down_e = fMT2tree->MinMetJetDPhi(0,40,5.0,-200);  

    //-----------------etau_ees_up--------------------------
    
    if(status == "etau_ees_up")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    
	    if(fPUReweight)
	      weight *= fMT2tree->pileUp.Weight;
	    weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	    if(Sample.sname == "Wtolnu")
	      weight *= fMT2tree->eleTau[0].tauWjetsSF;
       }

     
	if(MET_up_e.Pt() <= 30)
	  continue;
	if(MinMetJetDPhi_up_e <= 1)
	  continue;
	if(invMass_etau_up_e <= 15 || (invMass_etau_up_e >= 45 && invMass_etau_up_e <= 75))
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20)     
	  continue;
	if(ele_etau_up.Pt() <= 25)
	  continue;
	if(mt2_etau_up_e <= 90)
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	  continue;
   }

    //-----------------etau_ees_down--------------------------
    
    else if(status == "etau_ees_down")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    
	    if(fPUReweight)
	      weight *= fMT2tree->pileUp.Weight;
	    weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	    if(Sample.sname == "Wtolnu")
	      weight *= fMT2tree->eleTau[0].tauWjetsSF;
	  }
	

	if(MET_down_e.Pt() <= 30)
	  continue;
	if(MinMetJetDPhi_down_e <= 1)
	  continue;
	if(invMass_etau_down_e <= 15 || (invMass_etau_down_e >= 45 && invMass_etau_down_e <= 75))
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20)     
	  continue;
	if(ele_etau_down.Pt() <= 25)
	  continue;
	if(mt2_etau_down_e <= 90)
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	  continue;


      }
  }

  else if (status == "ee_binI_ees_up" || status == "ee_binI_ees_down" || status == "ee_binII_ees_up" || status == "ee_binII_ees_down") {

    //-----------------ees-------------------
    
    double px_all_up_e = 0;
    double py_all_up_e = 0;
    double pz_all_up_e = 0;
    double E_all_up_e  = 0;

    for(int l=0; l <fMT2tree->NEles; l++)
      {
	if (fabs(fMT2tree->ele[l].lv.Eta()) < 1.479)
	  {     
	    px_all_up_e += 0.01*fMT2tree->ele[l].lv.Px();
	    py_all_up_e += 0.01*fMT2tree->ele[l].lv.Py();
	    pz_all_up_e += 0.01*fMT2tree->ele[l].lv.Pz();
	    E_all_up_e  += 0.01*fMT2tree->ele[l].lv.E();
	  }
	else
	  {
	    px_all_up_e += 0.025*fMT2tree->ele[l].lv.Px();
	    py_all_up_e += 0.025*fMT2tree->ele[l].lv.Py();
	    pz_all_up_e += 0.025*fMT2tree->ele[l].lv.Pz();
	    E_all_up_e  += 0.025*fMT2tree->ele[l].lv.E();
	  }
      }
    
    TLorentzVector ele_delta_up    (px_all_up_e ,  py_all_up_e,  pz_all_up_e,  E_all_up_e);
    TLorentzVector ele_delta_down  (-px_all_up_e, -py_all_up_e, -pz_all_up_e, -E_all_up_e);
    
    TLorentzVector MET_up   = fMT2tree->pfmet[0] - ele_delta_up;
    TLorentzVector MET_down = fMT2tree->pfmet[0] - ele_delta_down;
    
    TLorentzVector ele0_ee_up  (0,0,0,0) ;
    TLorentzVector ele0_ee_down(0,0,0,0) ;

    TLorentzVector ele1_ee_up  (0,0,0,0) ;
    TLorentzVector ele1_ee_down(0,0,0,0) ;

    if (fabs(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Eta()) < 1.479)
      {
        ele0_ee_up.SetPxPyPzE(1.01 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Px(), 
			       1.01 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Py(), 
			       1.01 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Pz(), 
			       1.01 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.E()) ;
	
	ele0_ee_down.SetPxPyPzE(0.99 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Px(),
				 0.99 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Py(),
				 0.99 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Pz(),
				 0.99 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.E());
      }

    else
      {
	ele0_ee_up.SetPxPyPzE(1.025 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Px(), 
			       1.025 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Py(), 
			       1.025 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Pz(), 
			       1.025 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.E());
	
	ele0_ee_down.SetPxPyPzE(0.975 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Px(),
				 0.975 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Py(),
				 0.975 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Pz(),
				 0.975 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.E());
      }


    if (fabs(fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Eta()) < 1.479)
      {
        ele1_ee_up.SetPxPyPzE(1.01 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Px(), 
			       1.01 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Py(), 
			       1.01 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Pz(), 
			       1.01 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.E()) ;
	
	ele1_ee_down.SetPxPyPzE(0.99 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Px(),
				 0.99 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Py(),
				 0.99 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Pz(),
				 0.99 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.E());
      }

    else
      {
	ele1_ee_up.SetPxPyPzE(1.025 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Px(), 
			       1.025 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Py(), 
			       1.025 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Pz(), 
			       1.025 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.E());
	
	ele1_ee_down.SetPxPyPzE(0.975 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Px(),
				 0.975 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Py(),
				 0.975 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Pz(),
				 0.975 * fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.E());
      }


    double mt2_ee_up   = fMT2tree->CalcMT2(0, false, ele0_ee_up  , ele1_ee_up   , MET_up  );
    double mt2_ee_down = fMT2tree->CalcMT2(0, false, ele0_ee_down, ele1_ee_down , MET_down);

    double invMass_ee_up   = (ele0_ee_up   + ele1_ee_up  ).M(); 
    double invMass_ee_down = (ele0_ee_down + ele1_ee_down).M(); 

    double JZB_ee_up   = (-MET_up   - ele0_ee_up   - ele1_ee_up  ).Pt() - (ele0_ee_up   + ele1_ee_up  ).Pt() ;
    double JZB_ee_down = (-MET_down - ele0_ee_down - ele1_ee_down).Pt() - (ele0_ee_down + ele1_ee_down).Pt() ;

    double sumMT_ee_up = fMT2tree->GetMT(ele0_ee_up , ele0_ee_up.M(), MET_up , 0.0) + fMT2tree->GetMT(ele1_ee_up , ele1_ee_up.M(), MET_up , 0.0);
    double sumMT_ee_down = fMT2tree->GetMT(ele0_ee_down , ele0_ee_down.M(), MET_down , 0.0) + fMT2tree->GetMT(ele1_ee_down , ele1_ee_down.M(), MET_down , 0.0);

    double sumMT_ee_up_fixedmass = fMT2tree->GetMT(ele0_ee_up , 0.5 , MET_up , 0.0) + fMT2tree->GetMT(ele1_ee_up , ele1_ee_up.M(), MET_up , 0.0);
    double sumMT_ee_down_fixedmass = fMT2tree->GetMT(ele0_ee_down , 0.5, MET_down , 0.0) + fMT2tree->GetMT(ele1_ee_down , ele1_ee_down.M(), MET_down , 0.0);

	  if(variable == "MET")
	    myQuantity = fMT2tree->misc.MET;
	  if(variable == "MET_up")
	    myQuantity = MET_up.Pt();
	  if(variable == "MET_down")
	    myQuantity = MET_down.Pt();

    	  if(variable == "INVMASS")
	    myQuantity = fMT2tree->doubleEle[0].lv.M();
	  if(variable == "INVMASS_up")
	    myQuantity = invMass_ee_up;
	  if(variable == "INVMASS_down")
	    myQuantity = invMass_ee_down;

	  if(variable == "JZB")
	    myQuantity = fMT2tree->eeJZBInDirect();
	  if(variable == "JZB_up")
	    myQuantity = JZB_ee_up;
	  if(variable == "JZB_down")
	    myQuantity = JZB_ee_down;

	  if(variable == "MT2")
	   myQuantity = fMT2tree->doubleEle[0].MT2;
	  if(variable == "MT2_up")
	   myQuantity = mt2_ee_up;
	  if(variable == "MT2_down")
	   myQuantity = mt2_ee_down;

	  if(variable == "SUMMT-massive")
	    myQuantity = fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv , fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.M(), fMT2tree->pfmet[0] , 0.0) + fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv , fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.M(), fMT2tree->pfmet[0] , 0.0);
	  if(variable == "SUMMT-fixedmass")
	    myQuantity = fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv , 0.511, fMT2tree->pfmet[0] , 0.0) + fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv , 0.511, fMT2tree->pfmet[0] , 0.0);
          if(variable == "SUMMT-massless")
	    myQuantity = fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].MT + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].MT;

	  if(variable == "SUMMT_up")
	    myQuantity = sumMT_ee_up;
	  if(variable == "SUMMT_down")
	    myQuantity = sumMT_ee_down;
	  if(variable == "SUMMT_up_fixedmass")
	    myQuantity = sumMT_ee_up_fixedmass;
	  if(variable == "SUMMT_down_fixedmass")
	    myQuantity = sumMT_ee_down_fixedmass;



    if(status == "ee_binI_ees_up")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    weight *= fMT2tree->pileUp.Weight;
	    weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      
	  }

	if(ele0_ee_up.Pt() <= 20)
	  continue;
	if(ele1_ee_up.Pt() <= 10)
	  continue;

	if(MET_up.Pt() <= 30)
	  continue;
	if(invMass_ee_up <= 15 || (invMass_ee_up >= 71 && invMass_ee_up <= 111))
	  continue;
	if(JZB_ee_up >= -50)
	  continue;
	if(mt2_ee_up <= 90)
	  continue;
	if(sumMT_ee_up <= 250 || sumMT_ee_up >= 400 )     
 	  continue;
	
   }

    //-----------------ee_binI_ees_down--------------------------

   else if(status == "ee_binI_ees_down")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    weight *= fMT2tree->pileUp.Weight;
	    weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      
	  }

	if(ele0_ee_down.Pt() <= 20)
 	  continue;
 	if(ele1_ee_down.Pt() <= 10)
 	  continue;

 	if(MET_down.Pt() <= 30)
 	  continue;
 	if(invMass_ee_down <= 15 || (invMass_ee_down >= 71 && invMass_ee_down <= 111))
 	  continue;
	if(JZB_ee_down >= -50)
	  continue;
	if(mt2_ee_down <= 90)
	  continue;
 	if(sumMT_ee_down <= 250 || sumMT_ee_down >= 400 )     
	  continue;

      }

    //-----------------ee_binII_ees_up--------------------------
    
   else if(status == "ee_binII_ees_up")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    weight *= fMT2tree->pileUp.Weight;
	    weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      
	  }

 	if(ele0_ee_up.Pt() <= 20)
 	  continue;
 	if(ele1_ee_up.Pt() <= 10)
 	  continue;

 	if(MET_up.Pt() <= 30)
 	  continue;
 	if(invMass_ee_up <= 15 || (invMass_ee_up >= 71 && invMass_ee_up <= 111))
 	  continue;
	if(JZB_ee_up >= -50)
	  continue;
	if(mt2_ee_up <= 90)
	  continue;
 	if(sumMT_ee_up <= 400 )     
 	  continue;

   }

    //-----------------ee_binII_ees_down--------------------------

   else if(status == "ee_binII_ees_down")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    weight *= fMT2tree->pileUp.Weight;
	    weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      
	  }


 	if(ele0_ee_down.Pt() <= 20)
 	  continue;
 	if(ele1_ee_down.Pt() <= 10)
 	  continue;

 	if(MET_down.Pt() <= 30)
 	  continue;
 	if(invMass_ee_down <= 15 || (invMass_ee_down >= 71 && invMass_ee_down <= 111))
 	  continue;
	if(JZB_ee_down >= -50)
	  continue;
	if(mt2_ee_down <= 90)
	  continue;
 	if(sumMT_ee_down <= 400 )     
 	  continue;
      }


  }
  
  else
    {
    //----------------------ee_binI_nominal----------------------------
    if(status == "ee_binI_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    weight *= fMT2tree->pileUp.Weight;
	    weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      
	  }
	
 	if(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Pt() <=20)
 	  continue; 
 	if(fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Pt() <=10)
 	  continue; 

	if (fMT2tree->misc.MET <= 30)
 	  continue;
 	if(fMT2tree->doubleEle[0].lv.M() <= 15 || (fMT2tree->doubleEle[0].lv.M() >= 71 &&  fMT2tree->doubleEle[0].lv.M() <= 111))
 	  continue;
	if(fMT2tree->eeJZBInDirect() >= -50)
	  continue;
	if(fMT2tree->doubleEle[0].MT2 <= 90)                                                                     
	  continue;
	if((fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].MT + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].MT<=250 ) || (fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].MT + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].MT >=400))
	//	if((fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv , fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.M(), fMT2tree->pfmet[0] , 0.0) + fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv , fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.M(), fMT2tree->pfmet[0] , 0.0) <=250) ||(fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv , fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.M(), fMT2tree->pfmet[0] , 0.0) + fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv , fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.M(), fMT2tree->pfmet[0] , 0.0) >= 400))         
	  continue;
	
	
      }
    

    //----------------------ee_binII_nominal----------------------------
    else if(status == "ee_binII_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    weight *= fMT2tree->pileUp.Weight;
	    weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      
	  }
	
 	if(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.Pt() <=20)
 	  continue; 
 	if(fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.Pt() <=10)
 	  continue; 

 	if (fMT2tree->misc.MET <= 30)
 	  continue;
	if(fMT2tree->doubleEle[0].lv.M() <= 15 || (fMT2tree->doubleEle[0].lv.M() >= 71 &&  fMT2tree->doubleEle[0].lv.M() <= 111))
	  continue;
	if(fMT2tree->eeJZBInDirect()  >= -50)
	  continue;
	if(fMT2tree->doubleEle[0].MT2 <=  90)                                                                     
	  continue;
	//	if(fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv , fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].lv.M(), fMT2tree->pfmet[0] , 0.0) + fMT2tree->GetMT(fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv , fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].lv.M(), fMT2tree->pfmet[0] , 0.0) <= 400) 	
	if(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].MT + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].MT <= 400) 
          continue;

      }
    

    if(variable == "")
      myQuantity = 0;



    //----------------------ditau_bin1_nominal----------------------------

    else if(status == "ditau_bin1_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    
	    weight *= fMT2tree->pileUp.Weight;
	  }
	

	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;
	if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->doubleTau[0].GetMT2() <= 90)
	  continue;

	myQuantity = fMT2tree->doubleTau[0].GetMT2();

     }
    
    //----------------------ditau_bin2_nominal----------------------------
    
    else if(status == "ditau_bin2_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    weight *= fMT2tree->pileUp.Weight;
	  }
	
	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;
	if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->doubleTau[0].GetMT2() >= 90 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	  continue;
	if((fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT) + (fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT) <= 250)     
	  continue;

	myQuantity = fMT2tree->doubleTau[0].GetMT2();

      }     

    //----------------------mutau_nominal----------------------------

    else if(status == "mutau_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    
	    if(fPUReweight)
	      weight *= fMT2tree->pileUp.Weight;
	    weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();
	    
	    if(Sample.sname == "Wtolnu")
	      weight *= fMT2tree->muTau[0].GetTauWjetsSF();
	  }

	if (fMT2tree->misc.MET <= 30)
	  continue;
	if (fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->muTau[0].GetMT2() <= 90)
	  continue;
	if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	  continue;

	  myQuantity = fMT2tree->muTau[0].GetMT2();


      }     

    //----------------------etau_nominal----------------------------
    
    else if(status == "etau_nominal")
      {
	if(data == 1)
	  weight = 1.0;
	else
	  {
	    weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	    
	    if(fPUReweight)
	      weight *= fMT2tree->pileUp.Weight;
	    weight *=  fMT2tree->eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;
	    if(Sample.sname == "Wtolnu")
	      weight *= fMT2tree->eleTau[0].tauWjetsSF;
	  }
	
	if (fMT2tree->misc.MET <= 30)
	  continue;
	if (fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pt() <= 25 )
	  continue;
	if (fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->eleTau[0].GetMT2() <= 90)
	  continue;
	if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	  continue;

	myQuantity = fMT2tree->eleTau[0].GetMT2();

      }     

   //----------------------ditau_bin1_pu_up----------------------------

else if(status == "ditau_bin1_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;
	}

	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;
	if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	  continue;
	if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	  continue;
	if(fMT2tree->doubleTau[0].GetMT2() <= 90)
	  continue;
	
	     myQuantity = fMT2tree->doubleTau[0].GetMT2();

     }


   //----------------------------ditau_bin1_pu_down---------------------

else if(status == "ditau_bin1_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{

		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;

	}


	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;

	   if(fMT2tree->doubleTau[0].GetMT2() <= 90)
	     continue;

	   myQuantity = fMT2tree->doubleTau[0].GetMT2();
     }


   //----------------------ditau_bin2_pu_up----------------------------

else if(status == "ditau_bin2_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;

	}

	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if(fMT2tree->doubleTau[0].GetMT2() >= 90 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	     continue;
	   if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT <= 250)     
	     continue;

	     myQuantity = fMT2tree->doubleTau[0].GetMT2();

     }


   //----------------------------ditau_bin2_pu_down---------------------

else if(status == "ditau_bin2_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;


	}

	if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	  continue;

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if(fMT2tree->doubleTau[0].GetMT2() >= 90 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	     continue;
	   if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT <= 250)     
	     continue;

	   myQuantity = fMT2tree->doubleTau[0].GetMT2();
     }


   //----------------------mutau_pu_up----------------------------

else if(status == "mutau_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;
	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}

	if (fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Pt() <= 20 )
	  continue;

	   if ( fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree-> misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->muTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	     continue;


     myQuantity = fMT2tree->muTau[0].GetMT2();

     }


   //----------------------------mutau_pu_down---------------------

else if(status == "mutau_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;
	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}

	if (fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv.Pt() <= 20 )
	  continue;

	   if (fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->muTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	     continue;

	   myQuantity = fMT2tree->muTau[0].GetMT2();
    
     }


   //----------------------etau_pu_up----------------------------

else if(status == "etau_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;

	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;


	}

	if (fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pt() <= 25 )
	  continue;

	   if ( fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree-> misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->eleTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	     continue;


     myQuantity = fMT2tree->eleTau[0].GetMT2();

     }


   //----------------------------etau_pu_down---------------------

else if(status == "etau_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;

	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;

	}


	if (fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20 )
	  continue;
	if (fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pt() <= 25 )
	  continue;

	   if (fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->eleTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	     continue;

	   myQuantity = fMT2tree->eleTau[0].GetMT2();
     }
    
    }
//---------------------------------------------------------------------


//-------------------added----------------------------


      if(data == 1){
      
	MT2[8]->Fill(myQuantity, weight);//data
      
      }else{
	if(Sample.sname == "SUSY")
	  MT2[7]->Fill(myQuantity, weight);
	else
	  MT2[6]->Fill(myQuantity, weight);
      
	if(Sample.sname == "Top")
	  MT2[3]->Fill(myQuantity, weight);
	else
	  if(Sample.sname == "DY")	
	    MT2[2]->Fill(myQuantity, weight);
	  else
	    if(Sample.sname == "Wtolnu"){
	      MT2[1]->Fill(myQuantity, weight);}
	    else
	      if(Sample.sname == "QCD")
		MT2[0]->Fill(myQuantity, weight);
	      else
		if(Sample.sname == "VV")
		  MT2[4]->Fill(myQuantity, weight);
		else
		  if(Sample.sname == "Higgs")
		    MT2[5]->Fill(myQuantity, weight);
      }

	  AddOverAndUnderFlow(MT2[0] , true , true);
	  AddOverAndUnderFlow(MT2[1] , true , true);
	  AddOverAndUnderFlow(MT2[2] , true , true);
	  AddOverAndUnderFlow(MT2[3] , true , true);
	  AddOverAndUnderFlow(MT2[4] , true , true);
	  AddOverAndUnderFlow(MT2[5] , true , true);
	  AddOverAndUnderFlow(MT2[6] , true , true);
	  AddOverAndUnderFlow(MT2[7] , true , true);
	  AddOverAndUnderFlow(MT2[8] , true , true);

  

//  		  if(Sample.type=="susy"){

//  		    TString SUSYselection = TString::Format("(%s)", weight, myCuts.Data());
//  		    TString SUSYvariable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_PN_MLSP_MChi->GetName());
//  		    Sample.tree->Draw(SUSYvariable, SUSYselection, "goff");

//  		    SUSYselection = TString::Format(" (%s)", myCuts.Data()); 
//  		    SUSYvariable = TString::Format("Susy.MassLSP:Susy.MassGlu>>+%s", h_N_MLSP_MChi->GetName());
//  		    Sample.tree->Draw(SUSYvariable, SUSYselection, "goff");
		  
//  		  }


//-------------------added----------------------------^


//      if(data == 1)
// 	{
// 	  MT2[7]->Fill(myQuantity, weight);//data
// 	}
//      else
// 	{
//            if(Sample.sname == "SUSY")
//  	    MT2[6]->Fill(myQuantity, weight);
//  	  else
//  	    MT2[5]->Fill(myQuantity, weight);//MC

//               if(Sample.sname == "Top")
//  	       MT2[3]->Fill(myQuantity, weight);
// 	     else 
// 	       if(Sample.sname == "DY")	
// 		 MT2[2]->Fill(myQuantity, weight);
// 	       else
// 		 if(Sample.sname == "Wtolnu")
// 		   MT2[1]->Fill(myQuantity, weight);
// 		 else
// 		   if(Sample.sname == "QCD")
// 		     MT2[0]->Fill(myQuantity, weight);
// 		   else
// 		     if(Sample.sname == "VV")
// 		       MT2[4]->Fill(myQuantity, weight);




//               if(Sample.name == "DYToLL_M10To50")
// 		  MT2[0]->Fill(myQuantity, weight);
//               else 
// 	        if(Sample.name == "DYToLL_M50")
// 		  MT2[1]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY1ToLL_M50")	
// 		  MT2[2]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY2ToLL_M50")	
// 		  MT2[3]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY3ToLL_M50")	
// 		  MT2[4]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY4ToLL_M50")	
// 		  MT2[5]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo2L2Nu")	
// 		  MT2[6]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo2L2Q")	
// 		  MT2[7]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo4L")	
// 		  MT2[8]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "WZJetsTo3LNu")	
// 		  MT2[9]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "WZJetsTo2L2Q")	
// 		  MT2[10]->Fill(myQuantity, weight);

	
	

    
      
      }//For(events)

    double e1 =0;
   for(int m = 0; m <= NumberOfSamples; m++){

   e1 = MT2[m]->Integral();
   std::cout << setfill('#') << std::setw(50) << "" << std::endl;
   cout << "sample:" << cnames[m] << endl;
   cout << "all_content:..."<< e1 << endl;   

   } 
   
    }
//for(samples)



  for(int j = 0; j <= NumberOfSamples; j++){
    AddOverAndUnderFlow(MT2[j], true, true);
    //     AddOverAndUnderFlow(MT2[j], true, true);
     //          AddOverAndUnderFlow(MT2[j]);
//        MT2[j]->Multiply(pileup_data_up_histo);
   }


    THStack* h_stack = new THStack(varname, "");
    for(int j = 0; j <= NumberOfSamples; j++){
      // MT2[j]->Rebin(3);
       if(j <= (NumberOfSamples - 3))
        h_stack -> Add(MT2[j]);
    }

   //      Double a ;
   //      double b ;


//     for(int jj = 0; jj < NumberOfSamples; jj++){
//     std::cout << setfill('#') << std::setw(70) << "" << std::endl;
//      cout << "sample:" << cnames[jj] << endl;

//       a = 0;
//       b = 0;


//    for(int k = 0; k < 35; k++){
//     std::cout << setfill('=') << std::setw(70) << "" << std::endl;
//     //    cout << "bin" << k <<":" << bins[k]<<"-"<<bins[k+1] <<"..."<<"its content:" <<MT2[jj]->GetBinContent(k)<< "its error: " << MT2[jj]->GetBinError(k)<<endl;


//      a += MT2[jj]->GetBinContent(k);

//      b += MT2[jj]->GetBinError(k)*MT2[jj]->GetBinError(k);

//  } 
//    double c =0;
//    c = sqrt(b);
//    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
//    cout << "sample:" << cnames[jj] << endl;
//    cout << "all_bin_content:..."<< a << endl;   
//    cout << "all_bin_error:..."<< c << endl;   
//    }

      double e ;
      double f ;


      for(int kk = 0; kk <= NumberOfSamples; kk++){

      e = 0;
      f = 0;

   for(int mm = 0; mm < nbins+1 ; mm++){


     e += MT2[kk]->GetBinContent(mm);

     f += MT2[kk]->GetBinError(mm)*MT2[kk]->GetBinError(mm);

 } 
   double g = sqrt(f);

   std::cout << setfill('-') << std::setw(50) << "" << std::endl;
   cout << "sample:" << cnames[kk] << endl;
   cout << "all_bin_content:..."<< e << endl;   
   cout << "all_bin_error:..."<< g << endl;   
   }


  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W", "f");
  Legend1->AddEntry(MT2[2], "ZX", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW", "f");
  Legend1->AddEntry(MT2[5], "Higgs", "f"); 
  Legend1->AddEntry(MT2[7], "SMS", "l");
  Legend1->AddEntry(MT2[8], "data", "l");


//   TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
//   Legend1->AddEntry(MT2[0], "QCD", "f");
//   Legend1->AddEntry(MT2[1], "W+jets", "f");
//   Legend1->AddEntry(MT2[2], "Z+jets", "f");
//   Legend1->AddEntry(MT2[3], "Top", "f");
//   Legend1->AddEntry(MT2[4], "WWto2L2Nu", "f");
//   Legend1->AddEntry(MT2[6], "Susy", "l");
//   Legend1->AddEntry(MT2[7], "data", "l");


  TString fileName2 = fOutputDir;
  if(!fileName2.EndsWith("/")) fileName2 += "/";
  Util::MakeOutputDir(fileName2);
  fileName2 = fileName2 + myfileName +"_ee"+".root";
  TFile *savefile2 = new TFile(fileName2.Data(), "RECREATE");
  savefile2 ->cd();
  h_stack->Write();
  MT2[0]->Write();
  MT2[1]->Write();
  MT2[2]->Write();
  MT2[3]->Write();
  MT2[4]->Write();
  MT2[5]->Write();
  MT2[6]->Write();
  MT2[7]->Write();
  MT2[8]->Write();
  Legend1->Write();
  savefile2->Close();
  std::cout << "Saved histograms in " << savefile2->GetName() << std::endl;

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;



//   TH1F *h_PN_Bkg  = (TH1F*)MT2[NumberOfSamples-2]->Clone();
//   h_PN_Bkg->SetName("h_PN_Bkg");
//   h_PN_Bkg->Rebin(h_PN_Bkg->GetNbinsX());
  
//   TH1F *h_PN_Data  = (TH1F*)MT2[NumberOfSamples]->Clone();
//   h_PN_Data->SetName("h_PN_Data");
//   h_PN_Data->Rebin(h_PN_Data->GetNbinsX());

//   h_SMSEvents->Rebin2D(4, 4);
//   h_PN_MLSP_MChi->Divide(h_SMSEvents);
  
//   TH2* hXsec = (TH2*) TFile::Open("referenceXSecs.root")->Get("C1C1_8TeV_NLONLL_LSP");
//   h_PN_MLSP_MChi->Multiply(hXsec);

//   TString fileName1 = fOutputDir;
//   if(!fileName1.EndsWith("/")) fileName1 += "/";
//   Util::MakeOutputDir(fileName1);
//   fileName1 = fileName1 + "countingForExclusion_" + myfileName +"_Histos.root";
//   TFile *savefile1 = new TFile(fileName1.Data(), "RECREATE");
//   savefile1 ->cd();
  
//   h_stack->SetName("h_stack");
  
//   h_stack->Write();
//   h_PN_Bkg->Write();
//   h_PN_Data->Write();
//   h_PN_MLSP_MChi->Write();
//   h_N_MLSP_MChi->Write();
//   savefile1->Close();
	
//   std::cout << "Saved histograms in " << savefile1->GetName() << std::endl;


  printHisto(h_stack, MT2[8], MT2[6], MT2[7], Legend1 , myfileName, "hist", true, myfileName, "Events", -10, 0, -10, true);

  plotRatioStack(h_stack, MT2[6], MT2[8], MT2[7], true, false, myfileName, Legend1, myfileName, "Events", -10, 0, -10, true);


}



void MassPlotter::eeAnalysisTESpUsys2(TString cuts, TString trigger, unsigned int nevents, TString myfileName, TString giveStatus ){


  TH1::SetDefaultSumw2();

  TFile* pileup_data = new TFile("pileupData.root","READ");
  
  TH1F* pileup_data_histo = (TH1F*) pileup_data->Get("pileup");

  TFile* pileup_data_up = new TFile("pileupDataUp.root","READ");
  
  TH1F* pileup_data_up_histo = (TH1F*) pileup_data_up->Get("pileup");
 
  TFile* pileup_data_down = new TFile("pileupDataDown.root","READ");
  
  TH1F* pileup_data_down_histo = (TH1F*) pileup_data_down->Get("pileup");

  
  pileup_data_histo->Scale(1.0/pileup_data_histo->Integral());
  pileup_data_up_histo->Scale(1.0/pileup_data_up_histo->Integral());
  pileup_data_down_histo->Scale(1.0/pileup_data_down_histo->Integral());


  pileup_data_up_histo->Divide(pileup_data_histo);
  pileup_data_down_histo->Divide(pileup_data_histo);

  //   for(int i = 0; i < pileup_data_up_histo->GetNbinsX();i++){
  //     pileup_data_up_histo->SetBinContent(i+1,  pileup_data_up_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
//    }

//    for(int i = 0; i < pileup_data_down_histo->GetNbinsX();i++){
//      pileup_data_down_histo->SetBinContent(i+1,  pileup_data_down_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
//    }
  //-------------------added-------------------
  TString  cnames[NumberOfSamples+1] = {"QCD",     "W",    "ZX",  "Top",       "WW",  "Higgs",   "MC", "susy","data"};
  int      ccolor[NumberOfSamples+1] = { 401,      417,     419,   855,         603,     kRed,    603,      1, 632};
  //-------------------added-------------------

  //  TString  cnames[NumberOfSamples+1] = { "QCD", "Wjets", "Zjets", "Top", "WWjets",   "MC", "susy","data" };
  //  int      ccolor[NumberOfSamples+1] = { 401  ,     417,     419,   855,      603,    603,      1,  632  };

//  TString  cnames[11] = {"DYToLL_M10To50", "DYToLL_M50", "DY1ToLL_M50", "DY2ToLL_M50", "DY3ToLL_M50", "DY4ToLL_M50", "ZZJetsTo2L2Nu","ZZJetsTo2L2Q", "ZZJetsTo4L", "WZJetsTo3LNu", "WZJetsTo2L2Q"};
//   int      ccolor[11] = { 401, 417, 419, 855, 603, 650, 670, 840, 230, 584, 632};


  static const int nbins = 10;
  double bins[nbins+1] = {-2000.0, 0.0, 20.0, 40.0, 60.0, 90.0, 120.0, 150.0, 200.0, 500.0 , 2000.0};

  TString varname = "MT2";
  for (int i=0; i <= (NumberOfSamples) ; i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i],"",nbins, bins);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }
  //-------------------added------------------------
  MT2[NumberOfSamples] -> SetMarkerStyle(20);
  MT2[NumberOfSamples] -> SetMarkerColor(kBlack);
  MT2[NumberOfSamples] -> SetLineColor(kBlack);
  
  MT2[NumberOfSamples-2] -> SetFillStyle(3004);
  MT2[NumberOfSamples-2] -> SetFillColor(kBlack);
  //-------------------added-----------------------

//   MT2[7] -> SetMarkerStyle(20);
//   MT2[7] -> SetMarkerColor(kBlack);
//   MT2[7] -> SetLineColor(kBlack);
  
//   MT2[5] -> SetFillStyle(3004);
//   MT2[5] -> SetFillColor(kBlack);

   MT2[NumberOfSamples-1] -> SetMarkerStyle(20);
   MT2[NumberOfSamples-1] -> SetMarkerColor(kRed-4);
   MT2[NumberOfSamples-1] -> SetLineColor(kRed-4);
   MT2[NumberOfSamples-1] -> SetLineWidth(3);


  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++)
    {
 

    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];

    /////////////////////////
    //    if(Sample.sname != "DY")
    //      continue;
    //////////////////////////

    if(Sample.type == "data")
      {
      data = 1;
      myCuts += " && " + trigger;
      }

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();


    if (Sample.name == "QCD-Pt-15-20-MuEnriched")
      fPUReweight = false;
    else
      fPUReweight = true;

    float Weight=0.0;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    for ( unsigned int jentry=0 ; jentry<min(nentries , nevents);jentry++ )
      {

      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 )
	{
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	}

      if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")))
	{
	Weight = Sample.lumi;
	if(fPUReweight) 
        Weight /= Sample.PU_avg_weight;
	}

 

 
      float weight = Weight;

      double px_all_up = 0;
      double py_all_up = 0;
      double pz_all_up = 0;
      double E_all_up  = 0;

      double px_all_up_e = 0;
      double py_all_up_e = 0;
      double pz_all_up_e = 0;
      double E_all_up_e  = 0;


   for(int l=0; l <fMT2tree->NTaus; l++)
     {
       px_all_up += 0.03*fMT2tree->tau[l].lv.Px();
       py_all_up += 0.03*fMT2tree->tau[l].lv.Py();
       pz_all_up += 0.03*fMT2tree->tau[l].lv.Pz();
       E_all_up  += 0.03*fMT2tree->tau[l].lv.E();
     }

   for(int l=0; l <fMT2tree->NEles; l++)
     {
       if ( fabs(fMT2tree->ele[l].lv.Eta()) < 1.479)

	 {     
	   px_all_up_e += 0.01*fMT2tree->ele[l].lv.Px();
	   py_all_up_e += 0.01*fMT2tree->ele[l].lv.Py();
	   pz_all_up_e += 0.01*fMT2tree->ele[l].lv.Pz();
	   E_all_up_e  += 0.01*fMT2tree->ele[l].lv.E();
	 }
       else
	 {
	   px_all_up_e += 0.025*fMT2tree->ele[l].lv.Px();
	   py_all_up_e += 0.025*fMT2tree->ele[l].lv.Py();
	   pz_all_up_e += 0.025*fMT2tree->ele[l].lv.Pz();
	   E_all_up_e  += 0.025*fMT2tree->ele[l].lv.E();
	 }

     }



   TLorentzVector tau_delta_up    (px_all_up ,  py_all_up,  pz_all_up,  E_all_up);
   TLorentzVector tau_delta_down  (-px_all_up, -py_all_up, -pz_all_up, -E_all_up);

   TLorentzVector ele_delta_up    (px_all_up_e ,  py_all_up_e,  pz_all_up_e,  E_all_up_e);
   TLorentzVector ele_delta_down  (-px_all_up_e, -py_all_up_e, -pz_all_up_e, -E_all_up_e);
   

   TLorentzVector MET_up   = fMT2tree->pfmet[0] - tau_delta_up;
   TLorentzVector MET_down = fMT2tree->pfmet[0] - tau_delta_down;

   TLorentzVector MET_up_e   = fMT2tree->pfmet[0] - ele_delta_up;
   TLorentzVector MET_down_e = fMT2tree->pfmet[0] - ele_delta_down;


   TLorentzVector tau_mutau_up  (1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Px(), 
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pz(),  
                                 1.03 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_mutau_down (0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Px(),
                                  0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Py(),
	                          0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pz(),
                                  0.97 * fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_etau_up  (1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Px(), 
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pz(),  
                                 1.03 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau_etau_down (0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Px(),
                                  0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Py(),
	                          0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pz(),
                                  0.97 * fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector ele_etau_up ;
   TLorentzVector ele_etau_down ;
   if ( fabs(fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Eta()) < 1.479)
     {
       ele_etau_up.SetPxPyPzE   (1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(), 
                                 1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(), 
                                 1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(), 
			         1.01 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E()) ;

       ele_etau_down.SetPxPyPzE (0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(),
                                 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(),
                                 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(),
				 0.99 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());
     }

   else
     {
       ele_etau_up.SetPxPyPzE   (1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(), 
                                 1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(), 
                                 1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(), 
				 1.025 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());

       ele_etau_down.SetPxPyPzE (0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Px(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Py(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.Pz(),
				 0.975 * fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv.E());
     }


   TLorentzVector tau0_ditau_up (1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Px(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pz(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau0_ditau_down (0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Px(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Py(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pz(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.E());

   TLorentzVector tau1_ditau_up (1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Px(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Py(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pz(),
                                 1.03 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.E());

   TLorentzVector tau1_ditau_down (0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Px(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Py(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pz(),
                                   0.97 * fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.E());

   double tau_MT_mutau_up  = fMT2tree->GetMT( tau_mutau_up , 0., MET_up , 0.);
   double tau_MT_etau_up  = fMT2tree->GetMT( tau_etau_up , 0., MET_up , 0.);
   double tau0_MT_ditau_up = fMT2tree->GetMT( tau0_ditau_up , 0., MET_up , 0.); 
   double tau1_MT_ditau_up = fMT2tree->GetMT( tau1_ditau_up , 0., MET_up , 0.);                  

   double tau_MT_mutau_down  = fMT2tree->GetMT( tau_mutau_down , 0., MET_down , 0.);                            
   double tau_MT_etau_down  = fMT2tree->GetMT( tau_etau_down , 0., MET_down , 0.);                            
   double tau0_MT_ditau_down = fMT2tree->GetMT( tau0_ditau_down , 0., MET_down , 0.);
   double tau1_MT_ditau_down = fMT2tree->GetMT( tau1_ditau_down , 0., MET_down , 0.);                      
            
   double maxMT_up   = max(tau0_MT_ditau_up   , tau1_MT_ditau_up);
   double maxMT_down = max(tau0_MT_ditau_down , tau1_MT_ditau_down);
   double sumMT_up   = tau0_MT_ditau_up   + tau1_MT_ditau_up;
   double sumMT_down = tau0_MT_ditau_down + tau1_MT_ditau_down;

   double mt2_mutau_up = fMT2tree->CalcMT2(0, false, fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv, tau_mutau_up, MET_up);
   double mt2_mutau_down = fMT2tree->CalcMT2(0, false, fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv, tau_mutau_down, MET_down);
   double invMass_mutau_up  = (tau_mutau_up   + fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv).M(); 
   double invMass_mutau_down= (tau_mutau_down + fMT2tree->muo[fMT2tree->muTau[0].GetMuIndex0()].lv).M(); 

   double mt2_etau_up = fMT2tree->CalcMT2(0, false, fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv, tau_etau_up, MET_up);
   double mt2_etau_down = fMT2tree->CalcMT2(0, false, fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv, tau_etau_down, MET_down);
   double invMass_etau_up  = (tau_etau_up   + fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv).M(); 
   double invMass_etau_down= (tau_etau_down + fMT2tree->ele[fMT2tree->eleTau[0].GetEleIndex0()].lv).M(); 

   double mt2_etau_up_e = fMT2tree->CalcMT2(0, false, fMT2tree->tau[fMT2tree->eleTau[0].GetEleIndex0()].lv, ele_etau_up, MET_up_e);
   double mt2_etau_down_e = fMT2tree->CalcMT2(0, false, fMT2tree->tau[fMT2tree->eleTau[0].GetEleIndex0()].lv, ele_etau_down, MET_down_e);
   double invMass_etau_up_e  = (ele_etau_up   + fMT2tree->tau[fMT2tree->eleTau[0].GetEleIndex0()].lv).M(); 
   double invMass_etau_down_e= (ele_etau_down + fMT2tree->tau[fMT2tree->eleTau[0].GetEleIndex0()].lv).M(); 
   double MinMetJetDPhi_up_e   = fMT2tree->MinMetJetDPhi(0,40,5.0,200); 
   double MinMetJetDPhi_down_e = fMT2tree->MinMetJetDPhi(0,40,5.0,-200);  

   double mt2_ditau_up   = fMT2tree->CalcMT2(0, false, tau0_ditau_up, tau1_ditau_up, MET_up);
   double mt2_ditau_down = fMT2tree->CalcMT2(0, false, tau0_ditau_down, tau1_ditau_down, MET_down);
   double invMass_ditau_up    = (tau0_ditau_up + tau1_ditau_up).M();
   double invMass_ditau_down  = (tau0_ditau_down + tau1_ditau_down).M();

   double MinMetJetDPhi_up   = fMT2tree->MinMetJetDPhi(0,40,5.0,100); 
   double MinMetJetDPhi_down = fMT2tree->MinMetJetDPhi(0,40,5.0,-100);  


   TString status = giveStatus ;
   //   TString mutau_pu_up, mutau_pu_down, mutau_tes_up, mutau_tes_down, ditau_bin1_pu_up, ditau_bin1_pu_down, ditau_bin1_tes_up, ditau_bin1_tes_down, ditau_bin2_pu_up, ditau_bin2_pu_down ,ditau_bin2_tes_up, ditau_bin2_tes_down, mutau_nominal, ditau_bin1_nominal, ditau_bin2_nominal; 

   double myQuantity = 0;

   //----------------------mutau_nominal----------------------------

   if(status == "mutau_nominal")
     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
	    weight *= fMT2tree->pileUp.Weight;
	  weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}

	   if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].lv.Pt() <= 20)
             continue;
	   if (fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->muTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	     continue;


     myQuantity = fMT2tree->muTau[0].GetMT2();
     }     


   //-----------------mutau_tes_up--------------------------

else if(status == "mutau_tes_up")

     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight;
	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}

	   if (MET_up.Pt() <= 30)
	     continue;
	   if (MinMetJetDPhi_up <= 1)
	     continue;
	   if (invMass_mutau_up <= 15 || (invMass_mutau_up >= 45 && invMass_mutau_up <= 75))
	     continue;
	   if (tau_mutau_up.Pt() <= 20)
	     continue;
	   if (mt2_mutau_up <= 90)
	     continue;
	   if (tau_MT_mutau_up <= 200)
	     continue;
	   myQuantity = mt2_mutau_up;

     }

   //-----------------mutau_tes_down--------------------------

else if(status == "mutau_tes_down")
     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight;
	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}


	   if (MET_down.Pt() <= 30)
	     continue;
	   if (MinMetJetDPhi_down <= 1)
	     continue;
	   if ( invMass_mutau_down <= 15 || (invMass_mutau_down >= 45 & invMass_mutau_down <= 75))
	     continue;
	   if (tau_mutau_down.Pt() <= 20)
	     continue;
	   if (mt2_mutau_down <= 90)
	     continue;
	   if (tau_MT_mutau_down <= 200)
	     continue;
	   
	   myQuantity = mt2_mutau_down;

     }

   //----------------------ditau_bin1_nominal----------------------------

else if(status == "ditau_bin1_nominal")
     {

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	    weight *= fMT2tree->pileUp.Weight;
	}

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   
	   //   if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	   //	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
 	     continue;
 	   if(fMT2tree->doubleTau[0].GetMT2() <= 90)
 	     continue;

	   myQuantity = fMT2tree->doubleTau[0].GetMT2();//

	   //	   if(max(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT , fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT) <= 200)     	   //	     continue;
     
     }


   //-----------------ditau_bin1_tes_up--------------------------

else if(status == "ditau_bin1_tes_up")
     {      

       if(data == 1)
	 weight = 1.0;
       else
	 {
	
		    weight *= fMT2tree->pileUp.Weight;
	 }


       if (MinMetJetDPhi_up <= 1)
	 continue;
       if (tau0_ditau_up.Pt() <= 45 || tau1_ditau_up.Pt() <= 45)
	 continue;
       if (invMass_ditau_up <= 15 || (invMass_ditau_up >= 55  && invMass_ditau_up <= 85  ))
	 continue;
       if(mt2_ditau_up <= 90)
	 continue;

      myQuantity = mt2_ditau_up;

     }

   //-----------------ditau_bin1_tes_down--------------------------

else if(status == "ditau_bin1_tes_down")
     {      


            if(data == 1)
      	weight = 1.0;
            else
	      {
		weight *= fMT2tree->pileUp.Weight;
	      }


	    //     if (MinMetJetDPhi_down <= 1)
	    //       continue;
     
     if (invMass_ditau_down <= 15 || (invMass_ditau_down >= 55  && invMass_ditau_down <= 85  ))
       continue;
	    
     if (tau0_ditau_down.Pt() <= 45 || tau1_ditau_down.Pt() <= 45)
       continue;
	    
     if(mt2_ditau_down <= 90)
	      continue;

     myQuantity = mt2_ditau_down;

     }

   //----------------------ditau_bin2_nominal----------------------------

else if(status == "ditau_bin2_nominal")
     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  //	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	    weight *= fMT2tree->pileUp.Weight;
	}

    if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].lv.Pt() <= 45 || fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].lv.Pt() <= 45 )
	     continue;

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   //	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	   //	     continue;
	   if(fMT2tree->doubleTau[0].GetMT2() >= 90 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	     continue;
	   if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT <= 250)     
	     continue;


     myQuantity = fMT2tree->doubleTau[0].GetMT2();
     }     




   //-----------------ditau_bin2_tes_up--------------------------

else if(status == "ditau_bin2_tes_up")
     {      

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  //	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

		    weight *= fMT2tree->pileUp.Weight;
     }

	   //     if (MinMetJetDPhi_up <= 1)
	   //	continue;

     if (invMass_ditau_up <= 15 || (invMass_ditau_up >= 55  && invMass_ditau_up <= 85  ))
       continue;
     
     if (tau0_ditau_up.Pt() <= 45 || tau1_ditau_up.Pt() <= 45)
       continue;

     if(mt2_ditau_up >= 90 || mt2_ditau_up <= 40)
       continue;

      if (sumMT_up <= 250)
	continue;

     myQuantity = mt2_ditau_up;


     }

   //-----------------ditau_bin2_tes_down--------------------------

else if(status == "ditau_bin2_tes_down")
     {      

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  //		  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

		    weight *= fMT2tree->pileUp.Weight;
     }

	   //     if (MinMetJetDPhi_down <= 1)
	   //        continue;

     if (invMass_ditau_down <= 15 || (invMass_ditau_down >= 55  && invMass_ditau_down <= 85  ))
       continue;
     
     if (tau0_ditau_down.Pt() <= 45 || tau1_ditau_down.Pt() <= 45)
       continue;
	   
     if(mt2_ditau_down >= 90 || mt2_ditau_down <= 40)
       continue;
	   
     if (sumMT_down <= 250)
       continue;
	   
     myQuantity = mt2_ditau_down;

     }

   //----------------------mutau_pu_up----------------------------

else if(status == "mutau_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;
	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}

	   if ( fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree-> misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->muTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	     continue;


     myQuantity = fMT2tree->muTau[0].GetMT2();

     }


   //----------------------------mutau_pu_down---------------------

else if(status == "mutau_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;
	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

        if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();


	}

	   if (fMT2tree->muTau[0].GetLV().M() <= 15 || (fMT2tree->muTau[0].GetLV().M() >= 45 && fMT2tree->muTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->muTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->muTau[0].GetTauIndex0()].MT <= 200)     
	     continue;

	   myQuantity = fMT2tree->muTau[0].GetMT2();
     }


   //----------------------ditau_bin1_pu_up----------------------------

else if(status == "ditau_bin1_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;
	}


	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if(fMT2tree->doubleTau[0].GetMT2() <= 90)
	     continue;

	     myQuantity = fMT2tree->doubleTau[0].GetMT2();

     }


   //----------------------------ditau_bin1_pu_down---------------------

else if(status == "ditau_bin1_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{

		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;

	}


	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;

	   if(fMT2tree->doubleTau[0].GetMT2() <= 90)
	     continue;

	   myQuantity = fMT2tree->doubleTau[0].GetMT2();
     }


   //----------------------ditau_bin2_pu_up----------------------------

else if(status == "ditau_bin2_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;

	}

	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if(fMT2tree->doubleTau[0].GetMT2() >= 90 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	     continue;
	   if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT <= 250)     
	     continue;

	     myQuantity = fMT2tree->doubleTau[0].GetMT2();

     }


   //----------------------------ditau_bin2_pu_down---------------------

else if(status == "ditau_bin2_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;


	}


	   if (fMT2tree->doubleTau[0].GetLV().M() <= 15 || (fMT2tree->doubleTau[0].GetLV().M() >= 55 && fMT2tree->doubleTau[0].GetLV().M() <= 85))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if(fMT2tree->doubleTau[0].GetMT2() >= 90 || fMT2tree->doubleTau[0].GetMT2() <= 40)
	     continue;
	   if(fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT <= 250)     
	     continue;

	   myQuantity = fMT2tree->doubleTau[0].GetMT2();
     }


   //---etauj----
   //----------------------etau_nominal----------------------------

   if(status == "etau_nominal")
     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
	    weight *= fMT2tree->pileUp.Weight;

 	weight *=  fMT2tree->eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;

	}

	   if (fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;

	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->eleTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	     continue;
	   if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20)     
	     continue;
	   if(fMT2tree->ele[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 25)     
	     continue;



     myQuantity = fMT2tree->eleTau[0].GetMT2();
     }     

   //-----------------etau_tes_up--------------------------

else if(status == "etau_tes_up")

     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  //	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight;
	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;
	}

     if (MET_up.Pt() <= 30)
       continue;

     //     if (MinMetJetDPhi_up <= 1)
     //       continue;

     //     if ( invMass_etau_up <= 15 || (invMass_etau_up >= 45 && invMass_etau_up <= 75))
     //       continue;

     if (tau_etau_up.Pt() <= 20)
         continue;
     if (mt2_etau_up <= 90)
         continue;

     if (tau_MT_etau_up <= 200)
      continue;

     myQuantity = mt2_etau_up;

     }

   //-----------------etau_tes_down--------------------------

else if(status == "etau_tes_down")
     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  //	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight;

	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;

	}


     if (MET_down.Pt() <= 30)
       continue;

     //     if (MinMetJetDPhi_down <= 1)
     //       continue;
	   
     //     if ( invMass_etau_down <= 15 || (invMass_etau_down >= 45 & invMass_etau_down <= 75))
     //       continue;
	   
     if (tau_etau_down.Pt() <= 20)
       continue;
	   
     if (mt2_etau_down <= 90)
       continue;
	   
     if (tau_MT_etau_down <= 200)
       continue;
	   
     myQuantity = mt2_etau_down;

     }


   //-----------------etau_ees_up--------------------------

else if(status == "etau_ees_up")

     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight;
	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;
	}

     if (MET_up_e.Pt() <= 30)
       continue;

     if (MinMetJetDPhi_up_e <= 1)
       continue;

     if (invMass_etau_up_e <= 15 || (invMass_etau_up_e >= 45 && invMass_etau_up_e <= 75))
        continue;
     if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20)     
       continue;
     if (ele_etau_up.Pt() <= 25)
       continue;
     if (mt2_etau_up_e <= 90)
         continue;
     if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
       continue;

     myQuantity = mt2_etau_up_e;

     }

   //-----------------etau_ees_down--------------------------

else if(status == "etau_ees_down")
     {
           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight;

	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;

	}


     if (MET_down_e.Pt() <= 30)
       continue;

     if (MinMetJetDPhi_down_e <= 1)
       continue;

     if (invMass_etau_down_e <= 15 || (invMass_etau_down_e >= 45 && invMass_etau_down_e <= 75))
       continue;
     if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].lv.Pt() <= 20)     
       continue;
     if (ele_etau_down.Pt() <= 25)
       continue;
     if (mt2_etau_down_e <= 90)
         continue;
     if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
       continue;

     myQuantity = mt2_etau_down_e;


     }


   //----------------------etau_pu_up----------------------------

else if(status == "etau_pu_up")
     {
      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_up_Weight;

	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;


	}

	   if ( fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree-> misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->eleTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	     continue;


     myQuantity = fMT2tree->eleTau[0].GetMT2();

     }


   //----------------------------etau_pu_down---------------------

else if(status == "etau_pu_down")
     {
      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight * pu_down_Weight;

	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	if(Sample.sname == "Wtolnu")
	  weight *= fMT2tree->eleTau[0].tauWjetsSF;

	}

	   if (fMT2tree->eleTau[0].GetLV().M() <= 15 || (fMT2tree->eleTau[0].GetLV().M() >= 45 && fMT2tree->eleTau[0].GetLV().M() <= 75))
	     continue;
	   if (fMT2tree->misc.MinMetJetDPhiPt40 <= 1.0)
	     continue;
	   if (fMT2tree->misc.MET <= 30)
	     continue;
	   if(fMT2tree->eleTau[0].GetMT2() <= 90)
	     continue;
	   if(fMT2tree->tau[fMT2tree->eleTau[0].GetTauIndex0()].MT <= 200)     
	     continue;

	   myQuantity = fMT2tree->eleTau[0].GetMT2();
     }



//-------------------added----------------------------


      if(data == 1){
      
	MT2[NumberOfSamples]->Fill(myQuantity, weight);//data
      
      }else{
	if(Sample.sname == "SUSY")
	  MT2[NumberOfSamples-1]->Fill(myQuantity, weight);
	else
	  MT2[NumberOfSamples-2]->Fill(myQuantity, weight);
      
	if(Sample.sname == "Top")
	  MT2[3]->Fill(myQuantity, weight);
	else
	  if(Sample.sname == "DY")	
	    MT2[2]->Fill(myQuantity, weight);
	  else
	    if(Sample.sname == "Wtolnu"){
	      MT2[1]->Fill(myQuantity, weight);}
	    else
	      if(Sample.sname == "QCD")
		MT2[0]->Fill(myQuantity, weight);
	      else
		if(Sample.sname == "VV")
		  MT2[4]->Fill(myQuantity, weight);
		else
		  if(Sample.sname == "Higgs")
		    MT2[5]->Fill(myQuantity, weight);
      }

//-------------------added----------------------------^

//      if(data == 1)
// 	{
// 	  MT2[7]->Fill(myQuantity, weight);//data
// 	}
//      else
// 	{
//            if(Sample.sname == "SUSY")
//  	    MT2[6]->Fill(myQuantity, weight);
//  	  else
//  	    MT2[5]->Fill(myQuantity, weight);//MC

//               if(Sample.sname == "Top")
//  	       MT2[3]->Fill(myQuantity, weight);
// 	     else 
// 	       if(Sample.sname == "DY")	
// 		 MT2[2]->Fill(myQuantity, weight);
// 	       else
// 		 if(Sample.sname == "Wtolnu")
// 		   MT2[1]->Fill(myQuantity, weight);
// 		 else
// 		   if(Sample.sname == "QCD")
// 		     MT2[0]->Fill(myQuantity, weight);
// 		   else
// 		     if(Sample.sname == "VV")
// 		       MT2[4]->Fill(myQuantity, weight);


//               if(Sample.name == "DYToLL_M10To50")
// 		  MT2[0]->Fill(myQuantity, weight);
//               else 
// 	        if(Sample.name == "DYToLL_M50")
// 		  MT2[1]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY1ToLL_M50")	
// 		  MT2[2]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY2ToLL_M50")	
// 		  MT2[3]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY3ToLL_M50")	
// 		  MT2[4]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "DY4ToLL_M50")	
// 		  MT2[5]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo2L2Nu")	
// 		  MT2[6]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo2L2Q")	
// 		  MT2[7]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "ZZJetsTo4L")	
// 		  MT2[8]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "WZJetsTo3LNu")	
// 		  MT2[9]->Fill(myQuantity, weight);
// 	      else 
// 		if(Sample.name == "WZJetsTo2L2Q")	
// 		  MT2[10]->Fill(myQuantity, weight);

	
	

      }//For(events)

    double e1 =0;
    for(int m = 0; m <= (NumberOfSamples); m++){

   e1 = MT2[m]->Integral();
   std::cout << setfill('#') << std::setw(50) << "" << std::endl;
   cout << "sample:" << cnames[m] << endl;
   cout << "all_content:..."<< e1 << endl;   

   } 
   
}
//for(samples)


//-----------------added---------------------

  for(int j = 0; j <= (NumberOfSamples); j++){
    AddOverAndUnderFlow(MT2[j], true, true);
  }
//   printYield();

//----------------added---------------------^
      double e ;
      double f ;


      for(int kk = 0; kk <= (NumberOfSamples); kk++){

      e = 0;
      f = 0;

   for(int mm = 0; mm < nbins; mm++){


     e += MT2[kk]->GetBinContent(mm);

     f += MT2[kk]->GetBinError(mm)*MT2[kk]->GetBinError(mm);

 } 
   double g = sqrt(f);

   std::cout << setfill('-') << std::setw(50) << "" << std::endl;
   cout << "sample:" << cnames[kk] << endl;
   cout << "all_bin_content:..."<< e << endl;   
   cout << "all_bin_error:..."<< g << endl;   
   }

//----------------------added-------------------------------
  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j <= (NumberOfSamples); j++){
    // MT2[j]->Rebin(3);
    if(j <= (NumberOfSamples - 3))
      h_stack -> Add(MT2[j]);
  }





  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W", "f");
  Legend1->AddEntry(MT2[2], "ZX", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WW", "f");
  Legend1->AddEntry(MT2[5], "Higgs", "f"); 
  Legend1->AddEntry(MT2[NumberOfSamples-1], "SMS", "l");
  Legend1->AddEntry(MT2[NumberOfSamples], "data", "l");


  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + myfileName +"_Histos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  h_stack->Write();

  for(int j = 0; j <= (NumberOfSamples); j++)
    MT2[j]->Write();
  Legend1->Write();

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  printHisto(h_stack, MT2[NumberOfSamples], MT2[NumberOfSamples-2], MT2[NumberOfSamples-1], Legend1 , "MTC", "hist", true, "MT2", "Events", 0, -10, 2, true);

  plotRatioStack(h_stack, MT2[NumberOfSamples-2], MT2[NumberOfSamples], MT2[NumberOfSamples-1], true, false, "MT2_ratio", Legend1, "MT2", "Events", 0, -10, 2, true, "C");


//   TCanvas *MyC = new TCanvas("MyC", "MyC");
//   MyC->Divide(2,1);
//   MyC->cd(1);
//   TGraphAsymmErrors* sig1 = plotSig(MT2[NumberOfSamples-1], MT2[NumberOfSamples-2], "MT2", "Lower Cut", /*type = */ 2, /*sys = */ 0.1);
//   sig1->Draw("ACP");

//   MyC->cd(2);
//   TGraphAsymmErrors* sig2 = plotSig(MT2[NumberOfSamples-1], MT2[NumberOfSamples-2], "MT2", "Upper Cut", /*type = */ 2, /*sys = */ 0.1);
//   sig2->Draw("ACP");

//-----------------added---------------------^




//   for(int j = 0; j < NumberOfSamples; j++){
//     //     AddOverAndUnderFlow(MT2[j], true, true);
//           AddOverAndUnderFlow(MT2[j]);
// //        MT2[j]->Multiply(pileup_data_up_histo);
//    }


//     THStack* h_stack = new THStack(varname, "");
//     for(int j = 0; j < NumberOfSamples; j++){
//       // MT2[j]->Rebin(3);
//        if(j < (NumberOfSamples - 3))
//         h_stack -> Add(MT2[j]);
//     }

//    //      Double a ;
//    //      double b ;


// //     for(int jj = 0; jj < NumberOfSamples; jj++){
// //     std::cout << setfill('#') << std::setw(70) << "" << std::endl;
// //      cout << "sample:" << cnames[jj] << endl;

// //       a = 0;
// //       b = 0;


// //    for(int k = 0; k < 35; k++){
// //     std::cout << setfill('=') << std::setw(70) << "" << std::endl;
// //     //    cout << "bin" << k <<":" << bins[k]<<"-"<<bins[k+1] <<"..."<<"its content:" <<MT2[jj]->GetBinContent(k)<< "its error: " << MT2[jj]->GetBinError(k)<<endl;


// //      a += MT2[jj]->GetBinContent(k);

// //      b += MT2[jj]->GetBinError(k)*MT2[jj]->GetBinError(k);

// //  } 
// //    double c =0;
// //    c = sqrt(b);
// //    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
// //    cout << "sample:" << cnames[jj] << endl;
// //    cout << "all_bin_content:..."<< a << endl;   
// //    cout << "all_bin_error:..."<< c << endl;   
// //    }

//       double e ;
//       double f ;


//       for(int kk = 0; kk < NumberOfSamples; kk++){

//       e = 0;
//       f = 0;

//    for(int mm = 0; mm < nbins; mm++){


//      e += MT2[kk]->GetBinContent(mm);

//      f += MT2[kk]->GetBinError(mm)*MT2[kk]->GetBinError(mm);

//  } 
//    double g = sqrt(f);

//    std::cout << setfill('-') << std::setw(50) << "" << std::endl;
//    cout << "sample:" << cnames[kk] << endl;
//    cout << "all_bin_content:..."<< e << endl;   
//    cout << "all_bin_error:..."<< g << endl;   
//    }


//   TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
//   Legend1->AddEntry(MT2[0], "QCD", "f");
//   Legend1->AddEntry(MT2[1], "W+jets", "f");
//   Legend1->AddEntry(MT2[2], "Z+jets", "f");
//   Legend1->AddEntry(MT2[3], "Top", "f");
//   Legend1->AddEntry(MT2[4], "WWto2L2Nu", "f");
//   Legend1->AddEntry(MT2[6], "Susy", "l");
//   Legend1->AddEntry(MT2[7], "data", "l");


//   TString fileName2 = fOutputDir;
//   if(!fileName2.EndsWith("/")) fileName2 += "/";
//   Util::MakeOutputDir(fileName2);
//   fileName2 = fileName2 + myfileName +"_eeA"+".root";
//   TFile *savefile2 = new TFile(fileName2.Data(), "RECREATE");
//   savefile2 ->cd();
//   h_stack->Write();
//   MT2[0]->Write();
//   MT2[1]->Write();
//   MT2[2]->Write();
//   MT2[3]->Write();
//   MT2[4]->Write();
//   MT2[5]->Write();
//   MT2[6]->Write();
//   MT2[7]->Write();
//   Legend1->Write();
//   savefile2->Close();
//   std::cout << "Saved histograms in " << savefile2->GetName() << std::endl;

//   cout<<" trigger "<<trigger<<endl;
//   cout<<" cuts "<<cuts<<endl;

//   printHisto(h_stack, MT2[7], MT2[5], MT2[6], Legend1 , myfileName, "hist", true, myfileName, "Events", -10, 0, -10, true);

//   plotRatioStack(h_stack, MT2[5], MT2[7], MT2[6], true, false, myfileName, Legend1, myfileName, "Events", -10, 0, -10, true);


}


void MassPlotter::eeAnalysisPUsys(TString cuts, TString trigger, unsigned int nevents, TString myfileName){

  TH1::SetDefaultSumw2();

  TFile* pileup_data = new TFile("pileupData.root","READ");
  
  TH1F* pileup_data_histo = (TH1F*) pileup_data->Get("pileup");

  TFile* pileup_data_up = new TFile("pileupDataUp.root","READ");
  
  TH1F* pileup_data_up_histo = (TH1F*) pileup_data_up->Get("pileup");
 
  TFile* pileup_data_down = new TFile("pileupDataDown.root","READ");
  
  TH1F* pileup_data_down_histo = (TH1F*) pileup_data_down->Get("pileup");

  
  pileup_data_histo->Scale(1.0/pileup_data_histo->Integral());
  pileup_data_up_histo->Scale(1.0/pileup_data_up_histo->Integral());
  pileup_data_down_histo->Scale(1.0/pileup_data_down_histo->Integral());


  for(int i = 0; i < pileup_data_up_histo->GetNbinsX();i++){
    pileup_data_up_histo->SetBinContent(i+1,  pileup_data_up_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
  }

  for(int i = 0; i < pileup_data_down_histo->GetNbinsX();i++){
    pileup_data_down_histo->SetBinContent(i+1,  pileup_data_down_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
  }


  TString  cnames[NumberOfSamples+1] = {"QCD", "Wjets", "Zjets", "Top", "WWjets", "MC", "susy","data"};
  int      ccolor[NumberOfSamples+1] = { 401,       417,    419,   855,       603,  603,      1, 632};

  static const int nbins = 8;
  double bins[nbins+1] = {0.0,20.0,40.0,60.0,90.0,120.0,150.0,200.0,500.0};      //MT2

  TString varname = "MT2";
  for (int i=0; i <= (NumberOfSamples); i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i], "",nbins, bins);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[7] -> SetMarkerStyle(20);
  MT2[7] -> SetMarkerColor(kBlack);
  MT2[7] -> SetLineColor(kBlack);
  
  MT2[5] -> SetFillStyle(3004);
  MT2[5] -> SetFillColor(kBlack);

  MT2[6] -> SetMarkerStyle(20);
  MT2[6] -> SetMarkerColor(kRed-4);
  MT2[6] -> SetLineColor(kRed-4);
  MT2[6] -> SetLineWidth(3);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++)
    {
    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];

    if(Sample.type == "data")
      {
      data = 1;
      myCuts += " && " + trigger;
      }

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();


    if (Sample.name == "QCD-Pt-15-20-MuEnriched")
      fPUReweight = false;
    else
      fPUReweight = true;

    float Weight=0.0;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    for ( unsigned int jentry=0 ; jentry<min(nentries , nevents);jentry++ )
      {
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 )
	{
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	}

      if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")))
	{
	Weight = Sample.lumi;
	if(fPUReweight) 
        Weight /= Sample.PU_avg_weight;
	}

      int nbin_up = pileup_data_up_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_up_Weight = pileup_data_up_histo->GetBinContent(nbin_up);
	 //         cout << "pu_Weight" << fMT2tree->pileUp.Weight << endl;
	 //         cout << "pu_up_Weight" << pu_up_Weight << endl;

      int nbin_down = pileup_data_down_histo->FindBin(fMT2tree->pileUp.PUtrueNumInt);
      double pu_down_Weight = pileup_data_down_histo->GetBinContent(nbin_down);
	 //         cout << "pu_up_Weight" << pu_up_Weight << endl;
 
      float weight = Weight;

           if(data == 1)
     	weight = 1.0;
           else
     	{
	
	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  if(fPUReweight)
		    weight *= fMT2tree->pileUp.Weight;
	  //		weight *= fMT2tree->pileUp.Weight * pu_up_Weight;
	//			weight *= fMT2tree->pileUp.Weight * pu_down_Weight;

	//--------------------mutau weight-----------------------------

	//	weight *=  fMT2tree->muTau[0].GetTauEnergySF() * fMT2tree->muTau[0].GetMuIdSF() * fMT2tree->muTau[0].GetMuIsoSF() * fMT2tree->muTau[0].GetMuTrgSF() * fMT2tree->muTau[0].GetTauTrgSF();

	//        if(Sample.sname == "Wtolnu")
	  //	  weight *= fMT2tree->muTau[0].GetTauWjetsSF();

	//--------------------eletau weight--------------------
	//	weight *=  fMT2tree-> eleTau[0].tauTrgSF * fMT2tree->eleTau[0].eleTrgSF * fMT2tree->eleTau[0].eleIdIsoSF * fMT2tree->eleTau[0].tauEnergySF;

	//        if(Sample.sname == "Wtolnu")
	//	  weight *= fMT2tree->eleTau[0].tauWjetsSF;




			    //   weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;
	}




//-----------------------------------------------------------

//      double myQuantity = fMT2tree->muTau[0].GetMT2();
   //         float myQuantity = fMT2tree->doubleEle[0].MT2;

           double myQuantity =fMT2tree->doubleTau[0].GetMT2(); //fMT2tree->doubleTau[0].MT2;
  //      float myQuantity = fMT2tree->eleTau[0].GetMT2();

      if(data == 1)
	{
	  MT2[7]->Fill(myQuantity, weight);//data
	}
      else
	{
          if(Sample.sname == "SUSY")
	    MT2[6]->Fill(myQuantity, weight);
	  else
	    MT2[5]->Fill(myQuantity, weight);//MC

	     if(Sample.sname == "QCD")
	       MT2[0]->Fill(myQuantity, weight);
             else
	     if(Sample.sname == "Wtolnu")
	       MT2[1]->Fill(myQuantity, weight);
	     else
	     if(Sample.sname == "DY")	
	       MT2[2]->Fill(myQuantity, weight);
	     else
             if(Sample.sname == "Top")
	       MT2[3]->Fill(myQuantity, weight);
	     else
	     if(Sample.sname == "VV")
	       MT2[4]->Fill(myQuantity, weight);

      }
	
      }//for(events)


    }//for(samples)


   for(int j = 0; j <= NumberOfSamples; j++){
     AddOverAndUnderFlow(MT2[j]);
   }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[4], "WWto2L2Nu", "f");
  Legend1->AddEntry(MT2[6], "Susy", "l");
  Legend1->AddEntry(MT2[7], "data", "l");


   THStack* h_stack = new THStack(varname, "");
   for(int j = 0; j <= NumberOfSamples; j++){
     if (j==5)
       continue;
     //      if(j < (NumberOfSamples - 3))
//      MT2[j] -> SetFillColor (ccolor[j]);
//      MT2[j] -> SetLineColor (ccolor[j]);
//      MT2[j] -> SetLineWidth (2);
//      MT2[j] -> SetMarkerColor(ccolor[j]);
//      MT2[j] -> SetStats(false);
     h_stack -> Add(MT2[j]);
   }

     h_stack->Draw();
     Legend1->Draw("same");


      double a ;
      double b ;

      double e ;
      double f ;

//     for(int jj = 0; jj < NumberOfSamples; jj++){
//     std::cout << setfill('#') << std::setw(70) << "" << std::endl;
//      cout << "sample:" << cnames[jj] << endl;

//       a = 0;
//       b = 0;


//    for(int k = 0; k < 35; k++){
//     std::cout << setfill('=') << std::setw(70) << "" << std::endl;
//     //    cout << "bin" << k <<":" << bins[k]<<"-"<<bins[k+1] <<"..."<<"its content:" <<MT2[jj]->GetBinContent(k)<< "its error: " << MT2[jj]->GetBinError(k)<<endl;


//      a += MT2[jj]->GetBinContent(k);

//      b += MT2[jj]->GetBinError(k)*MT2[jj]->GetBinError(k);

//  } 
//    double c =0;
//    c = sqrt(b);
//    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
//    cout << "sample:" << cnames[jj] << endl;
//    cout << "all_bin_content:..."<< a << endl;   
//    cout << "all_bin_error:..."<< c << endl;   
//    }



      for(int kk = 0; kk <= (NumberOfSamples); kk++){

      e = 0;
      f = 0;

   for(int mm = 0; mm < 35; mm++){


     e += MT2[kk]->GetBinContent(mm);

     f += MT2[kk]->GetBinError(mm)*MT2[kk]->GetBinError(mm);

 } 
   double g = sqrt(f);

   std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   cout << "sample:" << cnames[kk] << endl;
   cout << "all_bin_content:..."<< e << endl;   
   cout << "all_bin_error:..."<< g << endl;   
   }



  TString fileName2 = fOutputDir;
  if(!fileName2.EndsWith("/")) fileName2 += "/";
  Util::MakeOutputDir(fileName2);
  fileName2 = fileName2 + myfileName +"_eeA"+".root";
  TFile *savefile2 = new TFile(fileName2.Data(), "RECREATE");
  savefile2 ->cd();
  h_stack->Write();
  MT2[0]->Write();
  MT2[1]->Write();
  MT2[2]->Write();
  MT2[3]->Write();
  MT2[4]->Write();
  MT2[5]->Write();
  MT2[6]->Write();
  MT2[7]->Write();
  Legend1->Write();
  savefile2->Close();
  std::cout << "Saved histograms in " << savefile2->GetName() << std::endl;

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  printHisto(h_stack, MT2[7], MT2[5], MT2[6], Legend1 , myfileName, "hist", true, myfileName, "Events", -10, 0, -10, true);

  plotRatioStack(h_stack, MT2[5], MT2[7], MT2[6], true, false, myfileName, Legend1, myfileName, "Events", -10, 0, -10, true);


}

void MassPlotter::eeAnalysisPUsys(){

  TH1::SetDefaultSumw2();

  //  TFile* susy_nominal = new TFile("results/etau-nominal_Yields.root","READ");
  TFile* susy_nominal = new TFile("results/etau-nominal-invMassAdded_Yields.root","READ");
  TH1F* susy_nominal_histo = (TH1F*) susy_nominal->Get("h_PN_MLSP_MChi");

  //  TFile* susy_tes_up = new TFile("results/etau-tes-up_Yields.root","READ");
  TFile* susy_tes_up = new TFile("results/etau-tes-up-invMassAdded_Yields.root","READ");
  TH1F* susy_tes_up_histo = (TH1F*) susy_tes_up->Get("h_PN_MLSP_MChi");
 

  //  TFile* susy_tes_down = new TFile("results/etau-tes-down_Yields.root","READ");
  TFile* susy_tes_down = new TFile("results/etau-tes-down-invMassAdded_Yields.root","READ");  
  TH1F* susy_tes_down_histo = (TH1F*) susy_tes_down->Get("h_PN_MLSP_MChi");


  susy_tes_up_histo->Add(susy_nominal_histo, -1);
  susy_tes_up_histo->Divide(susy_nominal_histo);

  susy_tes_down_histo->Add(susy_nominal_histo, -1);
  susy_tes_down_histo->Divide(susy_nominal_histo);

  TCanvas *Myc = new  TCanvas("a", "a");
  Myc->Divide(1,2);

  Myc->cd(1);
  susy_tes_up_histo->GetXaxis()->SetTitle("up/nominal");
  susy_tes_up_histo->Draw();

  Myc->cd(2);
  susy_tes_down_histo->GetXaxis()->SetTitle("down/nominal");
  susy_tes_down_histo->Draw();
  
  //  pileup_data_histo->Scale(1.0/pileup_data_histo->Integral());
  //  pileup_data_up_histo->Scale(1.0/pileup_data_up_histo->Integral());
  //  pileup_data_down_histo->Scale(1.0/pileup_data_down_histo->Integral());


  //  for(int i = 0; i < pileup_data_up_histo->GetNbinsX();i++){
  //    pileup_data_up_histo->SetBinContent(i+1,  pileup_data_up_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
  //  }

  //  for(int i = 0; i < pileup_data_down_histo->GetNbinsX();i++){
  //    pileup_data_down_histo->SetBinContent(i+1,  pileup_data_down_histo->GetBinContent(i+1)/pileup_data_histo->GetBinContent(i+1));
  //  }

}

void MassPlotter::eeAnalysisCDF(TString cuts, TString trigger, unsigned int nevents, TString myfileName){

  TH1::SetDefaultSumw2();

  TString  cnames[NumberOfSamples+1] = {"QCD", "Wjets", "Zjets", "Top", "WWjets", "MC", "susy","data"};
  int      ccolor[NumberOfSamples+1] = { 401,       417,    419,   855,       603,  603,      1, 632};

  TString varname = "MT2";
  for (int i=0; i<= (NumberOfSamples); i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i],"",35,0,3.5);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }
  
  MT2[7] -> SetMarkerStyle(20);
  MT2[7] -> SetMarkerColor(kBlack);
  MT2[7] -> SetLineColor(kBlack);
  
  MT2[5] -> SetFillStyle(3004);
  MT2[5] -> SetFillColor(kBlack);
  
  MT2[6] -> SetMarkerStyle(20);
  MT2[6] -> SetMarkerColor(kRed-4);
  MT2[6] -> SetLineColor(kRed-4);
  MT2[6] -> SetLineWidth(3);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++)
    {
    TString myCuts = cuts;
    
    int data = 0;
    sample Sample = fSamples[ii];

    if(Sample.type == "data")
      {
	data = 1;
      myCuts += " && " + trigger;
      }
    
    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
    
    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();
    

    if (Sample.name == "QCD-Pt-15-20-MuEnriched")
      fPUReweight = false;
    else
      fPUReweight = true;
    
    float Weight=0.0;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);
    
    for ( unsigned int jentry=0 ; jentry<min(nentries , nevents);jentry++ )
      {
	//Sample.tree->GetEntry(jentry);
	Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

	if ( fVerbose>2 && jentry % 100000 == 0 )
	  {
	    fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	  }
	
	if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")))
	  {
	Weight = Sample.lumi;
	if(fPUReweight) 
	  Weight /= Sample.PU_avg_weight;
	  }

	float weight = Weight;
	
	if(data == 1)
     	weight = 1.0;
	else
     	{
	  
	  //	  	  weight *= fMT2tree->SFWeight.BTagCSV40eq0; 

	  weight *= fMT2tree->pileUp.Weight;
	  
	  
	}
	
	

	double myQuantity =fMT2tree->misc.MinMetJetDPhiPt40;
	
	   if(data == 1)
	     {
	       MT2[7]->Fill(myQuantity, weight);//data
	     }
	   else
	     {
	       if(Sample.sname == "SUSY")
		 MT2[6]->Fill(myQuantity, weight);
	       else
		 MT2[5]->Fill(myQuantity, weight);//MC
	       
	       if(Sample.sname == "Top")
		 MT2[3]->Fill(myQuantity, weight);
	       else
		 if(Sample.sname == "DY")	
		   MT2[2]->Fill(myQuantity, weight);
		 else
		   if(Sample.sname == "Wtolnu")
		     MT2[1]->Fill(myQuantity, weight);
		   else
		     if(Sample.sname == "QCD")
		       MT2[0]->Fill(myQuantity, weight);
		     else
		       if(Sample.sname == "VV")
			 MT2[4]->Fill(myQuantity, weight);
	     }
	   
	   
      }//for(events)

    
    }//for(samples)

  

  for(int j = 0; j <= NumberOfSamples; j++){
    AddOverAndUnderFlow(MT2[j]);
  }
  
  MT2[6]->Scale(1.0/MT2[6]->Integral());
  MT2[4]->Scale(1.0/MT2[4]->Integral());
  

  TCanvas *c1 = new TCanvas("c1","");
  c1->Divide(1,2);
  
  TH1F *CDF_susy = new TH1F("CDF_susy","", 35, 0,3.5);
  TH1F *CDF_MC = new TH1F("CDF_MC","", 35,0,3.5);
  float susy_event = 0;
  float MC_event = 0;
  
  
  
  for (int i=0;i< 35;i++)
    {
      susy_event += MT2[6]->GetBinContent(35-i);
      CDF_susy->SetBinContent(35-i,  susy_event);
      
      MC_event += MT2[4]->GetBinContent(35-i);
      CDF_MC->SetBinContent(35-i,  MC_event);
      
    }
  
  c1->cd(1);
  MT2[6]->Draw();
  MT2[4]->Draw("same");
  
  c1->cd(2);
  CDF_susy->SetLineColor(kRed-4);
  CDF_susy->Draw();
  CDF_MC->Draw("same");
  

}


void MassPlotter::eeVS(Long64_t nevents, TString cuts, TString trigger){ 
   
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  //  TH2F *SignalDeltaM100 = new TH2F("SignalDeltaM100" , "SignalDeltaM100" , 100, -150,1000, 100,-150,1000);
  //  TH2F *SignalDeltaM250 = new TH2F("SignalDeltaM250" , "SignalDeltaM250" , 100, -150,1000, 100,-150,1000);
  TH2F *data  = new TH2F("data"     , "data"    , 30, 0, 3, 30, 0, 3);
  TH2F *mc    = new TH2F("mc"       , "mc"      , 30, 0, 3, 30, 0, 3);
  TH2F *susy  = new TH2F("susy"     , "susy"    , 30, 0, 3, 30, 0, 3);
  //  TH2F *Data    = new TH2F("d"     , "Bkg"    , 100, -150,1000, 100,-150,1000);
  //  TH2F *AvsB2Bkg    = new TH2F("AvsB2Bkg"   , "AvsB2Bkg"   , 100, -150,1000, 100,-150,1000);

  for(int i = 0; i <(int) fSamples.size(); i++){
    sample Sample = fSamples[i];
   
    //    if (Sample.type != "data")
    //      continue;
    TString myCuts = cuts;
    if( Sample.type=="data") myCuts += " && " + trigger;//just for completeness

    //    if((Sample.sname != "Wtolnu") && (Sample.type != "susy"))
    //      continue;
	   
    fMT2tree = new MT2tree();
 
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

   std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    Sample.tree->SetEventList(myEvtList);
    Long64_t nentries =  myEvtList->GetN();
    
    for (Long64_t jentry = 0; jentry < min(nentries, nevents); jentry++) {
     
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
      
      // float sumMT=fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].MT+fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].MT;

      //      if(Sample.type == "mc"){
      if(Sample.type == "data"){

	//	Bkg->Fill(fMT2tree->doubleEle[0].MT2, sumMT);
	//	Bkg->Fill(fMT2tree->PVisibleZetaMuTau(), fMT2tree->PZetaImbalancedMuTau());
	//	data->Fill(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Iso04, fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Iso04);
	data->Fill(fMT2tree->ele[0].Iso04, fMT2tree->ele[1].Iso04);
      }

      else if(Sample.type != "susy"){

	//	Bkg->Fill(fMT2tree->doubleEle[0].MT2, sumMT);
	//	Bkg->Fill(fMT2tree->PVisibleZetaMuTau(), fMT2tree->PZetaImbalancedMuTau());
	//	mc->Fill(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Iso04, fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Iso04);
	mc->Fill(fMT2tree->ele[0].Iso04, fMT2tree->ele[1].Iso04);
      }

      else {
	//	susy->Fill(fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Iso04, fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Iso04);
	susy->Fill(fMT2tree->ele[0].Iso04, fMT2tree->ele[1].Iso04);
      }


//       if(Sample.type == "susy"){
// 	 if(fMT2tree->misc.ProcessID!=10 || ((fMT2tree->Susy.MassGlu - fMT2tree->Susy.MassLSP > 75)  && (fMT2tree->Susy.MassGlu - fMT2tree->Susy.MassLSP) < 125))
// 	   SignalDeltaM100->Fill(fMT2tree->doubleEle[0].MT2, sumMT);
//          if(fMT2tree->misc.ProcessID!=10 || ((fMT2tree->Susy.MassGlu - fMT2tree->Susy.MassLSP > 225)  && (fMT2tree->Susy.MassGlu - fMT2tree->Susy.MassLSP) < 275))
// 	   SignalDeltaM250->Fill(fMT2tree->doubleEle[0].MT2, sumMT);

    }
  }
  

  TCanvas *MyC = new TCanvas("MyC", "MyC");
  MyC->Divide(2,2);
  MyC->cd(1);
  data->Draw();
  MyC->cd(2);
  mc->Draw();
  MyC->cd(3);
  susy->Draw();

  //  MyC->cd(2);
  //  SignalDeltaM100->Draw();
  //  MyC->cd(3);
  //  SignalDeltaM250->Draw();
  //  MyC->cd(4);
  //  AvsB2Signal->Draw();

}

void MassPlotter::miscEfficiency(int sample_index, unsigned int nevents){

  TH1::SetDefaultSumw2();
 	sample Sample = fSamples[sample_index];

	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	cout << "printing misc efficiency: \n"
	     << Sample.name << endl;	
	std::cout << setfill('-') << std::setw(70) << "" << std::endl;

	enum counters_t { count_begin, met0to30=count_begin, met30to60, met60to100, met100to150, met150to200, met200to1000,  count_end };
	Monitor counters[count_end];
	TString lablx[count_end] = {"met0to30", "met30to60" , "met60to100", "met100to150", "met150to200", "met200to1000"};
	//	TString lablx[count_end] = {"0<=MT2<40", "40<=MT2<80" , "80<=MT2<120", "120<=MT2"};

	fMT2tree = new MT2tree();
	Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
	unsigned int nentries =  Sample.tree->GetEntries();
	unsigned int nbytes = 0, nb = 0;
	for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
		nb = Sample.tree->GetEntry(jentry);   nbytes += nb;
		Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
		if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;

		bool acceptance(false);

		if (fMT2tree->misc.MinMetJetDPhiPt40 > 1) acceptance=true;

		Double_t weight = fMT2tree->pileUp.Weight;


		if(Sample.sname=="VV")
		  
		  {

		if(fMT2tree->misc.MET                     < 30    )
		  {counters[met0to30].fill("met0to30", weight);
		     if (acceptance) 
		       counters[met0to30].fill("acceptance", weight);}
		  
                if(fMT2tree->misc.MET                     >= 30    && fMT2tree->misc.MET                     < 60)		
		  {counters[met30to60].fill("met30to60", weight);
		     if (acceptance) 
		       counters[met30to60].fill("acceptance", weight);}
		  

                if(fMT2tree->misc.MET                     >= 60    && fMT2tree->misc.MET                     < 100)		
		  {counters[met60to100].fill("met60to100", weight);
		     if (acceptance) 
		       counters[met60to100].fill("acceptance", weight);}

                if(fMT2tree->misc.MET                     >= 100    && fMT2tree->misc.MET                     < 150)		
		  {counters[met100to150].fill("met100to150", weight);
		     if (acceptance) 
		       counters[met100to150].fill("acceptance", weight);}

                if(fMT2tree->misc.MET                     >= 150    && fMT2tree->misc.MET                     < 200)		
     		  {counters[met150to200].fill("met150to200", weight);
                     if (acceptance) 
		       counters[met150to200].fill("acceptance", weight);}

                if(fMT2tree->misc.MET                     >= 200    && fMT2tree->misc.MET                     < 1000)		
		  {counters[met200to1000].fill("met200to1000", weight);
		     if (acceptance) 
		       counters[met200to1000].fill("acceptance", weight);}

  }

      
else  {
  if(fMT2tree->misc.ProcessID!=10 || (abs(fMT2tree->Susy.MassGlu - 100) <= 5.0 && abs(fMT2tree->Susy.MassLSP - 0) <= 5.0))
       {
		  
		if(fMT2tree->misc.MET                     < 30    )
		  {counters[met0to30].fill("met0to30", weight);
		     if (acceptance) 
		       counters[met0to30].fill("acceptance", weight);}
		  
                if(fMT2tree->misc.MET                     >= 30    && fMT2tree->misc.MET                     < 60)		
		  {counters[met30to60].fill("met30to60", weight);
		     if (acceptance) 
		       counters[met30to60].fill("acceptance", weight);}
		  

                if(fMT2tree->misc.MET                     >= 60    && fMT2tree->misc.MET                     < 100)		
		  {counters[met60to100].fill("met60to100", weight);
		     if (acceptance) 
		       counters[met60to100].fill("acceptance", weight);}

                if(fMT2tree->misc.MET                     >= 100    && fMT2tree->misc.MET                     < 150)		
		  {counters[met100to150].fill("met100to150", weight);
		     if (acceptance) 
		       counters[met100to150].fill("acceptance", weight);}

                if(fMT2tree->misc.MET                     >= 150    && fMT2tree->misc.MET                     < 200)		
     		  {counters[met150to200].fill("met150to200", weight);
                     if (acceptance) 
		       counters[met150to200].fill("acceptance", weight);}

                if(fMT2tree->misc.MET                     >= 200    && fMT2tree->misc.MET                     < 1000)		
		  {counters[met200to1000].fill("met200to1000", weight);
		     if (acceptance) 
		       counters[met200to1000].fill("acceptance", weight);}
		//     }
       }
 }
	


	}

	// print stats     
	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	std::cout << "Statistics" << std::endl;
	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
	for ( counters_t iCount=count_begin; iCount<count_end; iCount = counters_t(iCount+1) ) {
		counters[iCount].print();
		std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;  
	}
	std::cout << setfill('=') << std::setw(70) << "" << std::endl;

	// fill histo
	TH1D* h_all = new TH1D("h_all", "", count_end, 0., (double) count_end );
	TH1D* h_acc = new TH1D("h_acc", "", count_end, 0., (double) count_end );

	TH1D* h_1   = new TH1D("h_1"  , "", count_end, 0., (double) count_end );

	for(int i=0; i<count_end; ++i){
		h_1    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);

		h_all->SetBinContent(i+1,counters[i].counts((string) lablx[i]));
		h_acc->SetBinContent(i+1,counters[i].counts("acceptance"));

	}

	h_1    ->Divide(h_acc,h_all);	

	TCanvas *col  = new TCanvas("miscEfficiecncy", "miscEfficiecncy");
	col -> cd();
// 	h_1    ->SetLineColor(kBlue);
// 	h_1    ->SetMarkerStyle(20);
// 	h_1    ->SetMinimum(0);
// 	h_1    ->SetMaximum(1);
// 	h_1    ->SetDrawOption("E");
	h_1    ->Draw();
	//_________________________________


	//  TCanvas *MyC = new TCanvas("MyC", "MyC");
	//  MyC->cd();  

  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";

  Util::MakeOutputDir(fileName);

  if(sample_index == 0)
  fileName = fileName + "miscEffvsMET_susy.root";

  if(sample_index == 1)
  fileName = fileName + "miscEffvsMET_Higgs.root";

  if(sample_index == 2)
  fileName = fileName + "miscEffvsMET_WW.root";

  if(sample_index == 3)
  fileName = fileName + "miscEffvsMET_ZZ2L2Nu.root";

  if(sample_index == 4)
  fileName = fileName + "miscEffvsMET_ZZTo2L2Q.root";

  if(sample_index == 5)
  fileName = fileName + "miscEffvsMET_ZZTo4L.root";

  if(sample_index == 6)
  fileName = fileName + "miscEffvsMET_WZTo3LNu.root";

  if(sample_index == 7)
  fileName = fileName + "miscEffvsMET_WZTo2L2Q.root";

  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile->cd();
  h_1->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;


// 	delete h_1;
// 	delete h_all;
// 	delete h_acc;
// 	delete col;
}

void MassPlotter::eeZInOut(TString cuts, TString trigger, unsigned int nevents, TString myfileName){

  TH1::SetDefaultSumw2();

  //  TString  cnames[NumberOfSamples+1] = {"QCD", "W", "ZX", "Top", "WW", "Higgs", "MC", "SUSY", "data" };
  //  int      ccolor[NumberOfSamples+1] = { 401 , 417,  419,   855,  603,    kRed,  603,      1,    632 };

  //    TString  cnames[3] = {"sumMT>250", "All", "Efficiency"};
  //    TString  cnames[3] = {"Zin", "Zout", "Outin_Efficiency"};
  //  int      ccolor1[NumberOfSamples+1] = { 401 , 417,  419,   855,  603,    kRed,  603,      1,    632 };
      //    double bins[nbins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 45.0, 50.0, 55.0, 60.0, 75.0, 90.0};

  static const int nbins = 6;
  double bins[nbins+1] = {0.0, 40, 50, 80, 110, 140, 200};


    TH1D* in  = new TH1D("in","in", nbins, bins);
    TH1D* out = new TH1D("out","out", nbins, bins);
    //    TH1D* out_a = new TH1D("efficiency","efficiency", nbins, bins);

  //  TString varname = "MT2";
    //  for (int i=0; i <= 12 ; i++){
    //    MT2[i] = new TH1D("","",nbins, bins);
//     MT2[i] -> SetFillColor (ccolor[i]);
//     MT2[i] -> SetLineColor (ccolor[i]);
//     MT2[i] -> SetLineWidth (2);
//     MT2[i] -> SetMarkerColor(ccolor[i]);
//     MT2[i] -> SetStats(false);
//  }

    MT2[2] = new TH1D("","",nbins, bins);
    MT2[12] = new TH1D("","",nbins, bins);
//   MT2[NumberOfSamples] -> SetMarkerStyle(20);
//   MT2[NumberOfSamples] -> SetMarkerColor(kBlack);
//   MT2[NumberOfSamples] -> SetLineColor(kBlack);
  
//   MT2[NumberOfSamples-2] -> SetFillStyle(3004);
//   MT2[NumberOfSamples-2] -> SetFillColor(kBlack);
  
//   MT2[NumberOfSamples-1] -> SetMarkerStyle(20);
//   MT2[NumberOfSamples-1] -> SetMarkerColor(kRed-4);
//   MT2[NumberOfSamples-1] -> SetLineColor(kRed-4);
//   MT2[NumberOfSamples-1] -> SetLineWidth(3);

//   TString varname1 = "MT21";
//   vector<TH1D*> MT21; 
//   for (int i=0; i <= (NumberOfSamples) ; i++){
//     MT21[i] = new TH1D(varname1+"_"+cnames[i],"",nbins, bins);
//     MT21[i] -> SetFillColor (ccolor[i]);
//     MT21[i] -> SetLineColor (ccolor[i]);
//     MT21[i] -> SetLineWidth (2);
//     MT21[i] -> SetMarkerColor(ccolor[i]);
//     MT21[i] -> SetStats(false);
//   }

//   MT21[NumberOfSamples] -> SetMarkerStyle(20);
//   MT21[NumberOfSamples] -> SetMarkerColor(kBlack);
//   MT21[NumberOfSamples] -> SetLineColor(kBlack);
  
//   MT21[NumberOfSamples-2] -> SetFillStyle(3004);
//   MT21[NumberOfSamples-2] -> SetFillColor(kBlack);
  
//   MT21[NumberOfSamples-1] -> SetMarkerStyle(20);
//   MT21[NumberOfSamples-1] -> SetMarkerColor(kRed-4);
//   MT21[NumberOfSamples-1] -> SetLineColor(kRed-4);
//   MT21[NumberOfSamples-1] -> SetLineWidth(3);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++)
    {



      TString myCuts = cuts;
      
      int data = 0;
      sample Sample = fSamples[ii];

    /////////////////////////
    if(Sample.sname != "DY")
      continue;
    //////////////////////////


      if(Sample.type == "data")
	{
	  data = 1;
	  myCuts += " && " + trigger;
	}

      fMT2tree = new MT2tree();
      Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      std::cout << setfill('=') << std::setw(70) << "" << std::endl;
      cout << "looping over :     " <<endl;	
      cout << "   Name:           " << Sample.name << endl;
      cout << "   File:           " << (Sample.file)->GetName() << endl;
      cout << "   Events:         " << Sample.nevents  << endl;
      cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
      cout << "   Xsection:       " << Sample.xsection << endl;
      cout << "   kfactor:        " << Sample.kfact << endl;
      cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
      cout << "   Weight:         " << Weight <<endl;
      std::cout << setfill('-') << std::setw(70) << "" << std::endl;
      
      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

      unsigned int nentries = myEvtList->GetN();
      float Weight=0.0;

	Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);

      for ( unsigned int jentry=0 ; jentry<min(nentries , nevents);jentry++ )
	{

	  Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

	  if (fVerbose>2 && jentry % 100000 == 0 )
	    {
	      fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	      fflush(stdout);
	    }

	  if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")))
	    {
	      Weight = Sample.lumi;
	      Weight /= Sample.PU_avg_weight;
	    }

	  float weight = Weight;
	  
	  if(data == 1)
	    weight = 1.0;
	  else
	    {
	
	      weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	      weight *= fMT2tree->pileUp.Weight;
	      weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      

	    }


	  double myQuantity = fMT2tree->doubleEle[0].MT2;


	  //	  if("((fMT2tree->doubleEle[0].Isolated == 1) && ((fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Charge + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Charge) == 0)) ")                                                           


  //	  if((fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT + fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT) > 250.0)
	    //	    if(fMT2tree->misc.MinMetJetDPhiPt40 > 1)

	  	  if(fMT2tree->doubleEle[0].lv.M() >= 86 && fMT2tree->doubleEle[0].lv.M() <= 96)     

	    {
		 
//  	      if(data == 1){
		
//  	      MT2[0]->Fill(myQuantity, weight);//data

//  	      }
//  	      else{
		
//  		if(Sample.sname != "SUSY")
//  		  MT2[1]->Fill(myQuantity, weight);
		
//		if(Sample.sname == "DY")	
		  MT2[2]->Fill(myQuantity, weight);
		
		//	      }
	    }


		  //	  if("(MT2tree->ee_e0Loose01to04() &&  ((fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Charge + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Charge) ==  2 ||  ((fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Charge + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Charge) == -2)))");

	  if((fMT2tree->doubleEle[0].lv.M() > 20 && fMT2tree->doubleEle[0].lv.M() < 86) || fMT2tree->doubleEle[0].lv.M() > 96)
	    {
	      
// 	      if(data == 1){
		
// 		MT2[10]->Fill(myQuantity, weight);//data
		
// 	      }
// 	      else{
// 		if(Sample.sname != "SUSY")
// 		  MT2[11]->Fill(myQuantity, weight);
		
//		if(Sample.sname == "DY")	
		  MT2[12]->Fill(myQuantity, weight);
		
		//	      }
	    }
		

// 	  AddOverAndUnderFlow(MT2[0] , true , true);
// 	  AddOverAndUnderFlow(MT2[1] , true , true);
 	  AddOverAndUnderFlow(MT2[2] , true , true);
// 	  AddOverAndUnderFlow(MT2[10] , true , true);
// 	  AddOverAndUnderFlow(MT2[11] , true , true);
 	  AddOverAndUnderFlow(MT2[12] , true , true);

	}
      
    }  //for(samples)
  

  //  in->Add(MT2[0], +1);
  //  in->Add(MT2[1], -1);
  in->Add(MT2[2], +1);

  //  out_b->Add(MT2[10], +1);
  //  out_b->Add(MT2[11], -1);
  out->Add(MT2[12], +1);


  out->Divide(in);



//  AddOverAndUnderFlow(out, true, true);
//  AddOverAndUnderFlow(in, true, true);

//  AddOverAndUnderFlow(MT2[0], true, true);
//  AddOverAndUnderFlow(MT2[1], true, true);



  //  TFile *f = new TFile("hsimple.root");
 
  //  TH1F *hpx = (TH1F*)f->Get("hpx");
 
  //create a function with 3 parameters in the range [-3,3]
  //  TF1 *func = new TF1("fit",fitf,-3,3,3);
  //  func->SetParameters(500,hpx->GetMean(),hpx->GetRMS());
  //  func->SetParNames("Constant","Mean_value","Sigma");


//  double out1 = out->Integral();
//  double in1 = in->Integral();
//  cout << "sumMtgt250 / " << out1 / in1 << endl ;
//  MT2[2]->Divide(MT2[0],MT2[1]);

  //  MT2[3]->Multiply(MT2[1],MT2[2]);

  //  cout << "Wjets estimation:.." << MT2[3]->Integral() << endl;

//     double e = 0;
//     double f = 0;
//     double g = 0;
//     double h = 0;
//     double i = 0;
//     double j = 0;

//   for(int kk = 0; kk <= 1 ; kk++){
    
    
//     for(int mm = 0; mm < nbins+1 ; mm++){
      
//       e = in->GetBinContent(mm);      
//       f = out_b->GetBinContent(mm);

//       g = in->GetBinError(mm);
//       h = out_b->GetBinError(mm);
      
//       cout << cnames[kk]<<": " <<bins[mm] <<"-" << bins[mm+1] <<":.. " << e <<"__error:.."<< g << endl;
//       cout << cnames[kk]<<": " <<bins[mm] <<"-" << bins[mm+1] <<":.. " << f <<"__error:.."<< h << endl;

//     } 

//   }


//   double in_yield, out_b_yield, eff_yield = 0.0;
//   double in_err, out_b_err, eff_err = 0.0 ;

//   in_yield = in -> IntegralAndError(0, nbins, in_err) ;
//   out_b_yield = out_b -> IntegralAndError(0, nbins, out_b_err) ;

//   cout << cnames[0] << "-yield:.." << in_yield  <<"__error:..." << in_err  << endl ;
//   cout << cnames[1] << "-yield:.." << out_b_yield <<"__error:.."  << out_b_err << endl ;

//  out_a->Divide(in, out_b,1,1);

//     for(int mm = 0; mm < nbins+1 ; mm++){
      
//       i = out_a->GetBinContent(mm);
//       j = out_a->GetBinError(mm)  ;

//       cout << cnames[2]<<": " <<bins[mm] <<"-" << bins[mm+1] <<":.. " << i <<"__error:.."<< j << endl;

//     } 

//     eff_yield = out_a -> IntegralAndError(0, nbins, eff_err); 
//   cout << cnames[2] << "-yield:.." << eff_yield <<"__error:.."  << eff_err << endl ;

  out->Fit("pol2");
  out->Fit("pol1");
  out->Fit("pol0");


  out->Draw();

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;
    
  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + myfileName +".root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  in->Write();
  out->Write();
  //  out_a->Write();

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;


}

void MassPlotter::eeQCDCtoBRatio(TString cuts, TString trigger, unsigned int nevents, TString myfileName){

  TH1::SetDefaultSumw2();

  //  TString  cnames[NumberOfSamples+1] = {"QCD", "W", "ZX", "Top", "WW", "Higgs", "MC", "SUSY", "data" };
  //  int      ccolor[NumberOfSamples+1] = { 401 , 417,  419,   855,  603,    kRed,  603,      1,    632 };

  //    TString  cnames[3] = {"sumMT>250", "All", "Efficiency"};
    TString  cnames[3] = {"B", "C", "CtoB_Ratio"};
  //  int      ccolor1[NumberOfSamples+1] = { 401 , 417,  419,   855,  603,    kRed,  603,      1,    632 };
      //    double bins[nbins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 45.0, 50.0, 55.0, 60.0, 75.0, 90.0};

    //  static const int nbins = 8;
    //  double bins[nbins+1] = {40.0, 50.0 , 60.0, 70, 80.0};

  static const int nbins = 5;
  //  double bins[nbins+1] = {40, 60, 90, 120, 150, 200};
  double bins[nbins+1] = {40, 50, 60, 70, 80, 90};

    TH1D* B_region  = new TH1D("B_region","B_region", nbins, bins);
    TH1D* C_region  = new TH1D("C_region","C_region", nbins, bins);
   //    TH1D* out_a = new TH1D("efficiency","efficiency", nbins, bins);

  //  TString varname = "MT2";
  for (int i=0; i <= 12 ; i++){
    MT2[i] = new TH1D("","",nbins, bins);
//     MT2[i] -> SetFillColor (ccolor[i]);
//     MT2[i] -> SetLineColor (ccolor[i]);
//     MT2[i] -> SetLineWidth (2);
//     MT2[i] -> SetMarkerColor(ccolor[i]);
//     MT2[i] -> SetStats(false);
  }

//   MT2[NumberOfSamples] -> SetMarkerStyle(20);
//   MT2[NumberOfSamples] -> SetMarkerColor(kBlack);
//   MT2[NumberOfSamples] -> SetLineColor(kBlack);
  
//   MT2[NumberOfSamples-2] -> SetFillStyle(3004);
//   MT2[NumberOfSamples-2] -> SetFillColor(kBlack);
  
//   MT2[NumberOfSamples-1] -> SetMarkerStyle(20);
//   MT2[NumberOfSamples-1] -> SetMarkerColor(kRed-4);
//   MT2[NumberOfSamples-1] -> SetLineColor(kRed-4);
//   MT2[NumberOfSamples-1] -> SetLineWidth(3);

//   TString varname1 = "MT21";
//   vector<TH1D*> MT21; 
//   for (int i=0; i <= (NumberOfSamples) ; i++){
//     MT21[i] = new TH1D(varname1+"_"+cnames[i],"",nbins, bins);
//     MT21[i] -> SetFillColor (ccolor[i]);
//     MT21[i] -> SetLineColor (ccolor[i]);
//     MT21[i] -> SetLineWidth (2);
//     MT21[i] -> SetMarkerColor(ccolor[i]);
//     MT21[i] -> SetStats(false);
//   }

//   MT21[NumberOfSamples] -> SetMarkerStyle(20);
//   MT21[NumberOfSamples] -> SetMarkerColor(kBlack);
//   MT21[NumberOfSamples] -> SetLineColor(kBlack);
  
//   MT21[NumberOfSamples-2] -> SetFillStyle(3004);
//   MT21[NumberOfSamples-2] -> SetFillColor(kBlack);
  
//   MT21[NumberOfSamples-1] -> SetMarkerStyle(20);
//   MT21[NumberOfSamples-1] -> SetMarkerColor(kRed-4);
//   MT21[NumberOfSamples-1] -> SetLineColor(kRed-4);
//   MT21[NumberOfSamples-1] -> SetLineWidth(3);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(unsigned int ii = 0; ii < fSamples.size(); ii++)
    {



      TString myCuts = cuts;
      
      int data = 0;
      sample Sample = fSamples[ii];

    /////////////////////////
      //    if(Sample.sname != "Wtolnu")
      //      continue;
    //////////////////////////


      if(Sample.type == "data")
	{
	  data = 1;
	  myCuts += " && " + trigger;
	}

      fMT2tree = new MT2tree();
      Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      std::cout << setfill('=') << std::setw(70) << "" << std::endl;
      cout << "looping over :     " <<endl;	
      cout << "   Name:           " << Sample.name << endl;
      cout << "   File:           " << (Sample.file)->GetName() << endl;
      cout << "   Events:         " << Sample.nevents  << endl;
      cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
      cout << "   Xsection:       " << Sample.xsection << endl;
      cout << "   kfactor:        " << Sample.kfact << endl;
      cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
      cout << "   Weight:         " << Weight <<endl;
      std::cout << setfill('-') << std::setw(70) << "" << std::endl;
      
      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

      unsigned int nentries = myEvtList->GetN();
      float Weight=0.0;

	Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);

      for ( unsigned int jentry=0 ; jentry<min(nentries , nevents);jentry++ )
	{

	  Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

	  if (fVerbose>2 && jentry % 100000 == 0 )
	    {
	      fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	      fflush(stdout);
	    }

	  if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50")))
	    {
	      Weight = Sample.lumi;
	      Weight /= Sample.PU_avg_weight;
	    }

	  float weight = Weight;
	  
	  if(data == 1)
	    weight = 1.0;
	  else
	    {
	
	      weight *= fMT2tree->SFWeight.BTagCSV40eq0; 
	      weight *= fMT2tree->pileUp.Weight;
	      weight *= fMT2tree->doubleEle[0].Ele0IdIsoSF * fMT2tree->doubleEle[0].Ele1IdIsoSF * fMT2tree->doubleEle[0].DiEleTrgSF;	      

	    }


	  double myQuantity = fMT2tree->doubleEle[0].MT2;


	  if(((fMT2tree->doubleEle[0].Isolated == -1) && ((fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Charge + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Charge) == 0)))
                                                           

  //	  if((fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].MT + fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].MT) > 250.0)
	    //	    if(fMT2tree->misc.MinMetJetDPhiPt40 > 1)

	  //	  if(fMT2tree->doubleEle[0].lv.M() >= 81 && fMT2tree->doubleEle[0].lv.M() <= 101)     

	    {
		 
	      if(data == 1){
		
	      MT2[0]->Fill(myQuantity, weight);//data

	      }
	      else{
		
		if(Sample.sname != "SUSY")
		  MT2[1]->Fill(myQuantity, weight);
		
		if(Sample.sname == "QCD")	
		  MT2[2]->Fill(myQuantity, weight);
		
	      }
	    }


	  if(((fMT2tree->doubleEle[0].Isolated == -1) &&  (((fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Charge + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Charge) ==  2) ||  ((fMT2tree->ele[fMT2tree->doubleEle[0].Ele0Ind].Charge + fMT2tree->ele[fMT2tree->doubleEle[0].Ele1Ind].Charge) == -2))));

	  //	  if((fMT2tree->doubleEle[0].lv.M() > 15 && fMT2tree->doubleEle[0].lv.M() < 81) || fMT2tree->doubleEle[0].lv.M() > 101)

	    {
	      
	      if(data == 1){
		
		MT2[10]->Fill(myQuantity, weight);//data
		
	      }
	      else{
		if(Sample.sname != "SUSY")
		  MT2[11]->Fill(myQuantity, weight);
		
		if(Sample.sname == "QCD")	
		  MT2[12]->Fill(myQuantity, weight);
		
	      }
	    }
		

 	  AddOverAndUnderFlow(MT2[0] , true , true);
 	  AddOverAndUnderFlow(MT2[1] , true , true);
 	  AddOverAndUnderFlow(MT2[2] , true , true);
 	  AddOverAndUnderFlow(MT2[10] , true , true);
 	  AddOverAndUnderFlow(MT2[11] , true , true);
 	  AddOverAndUnderFlow(MT2[12] , true , true);

	}
      
    }  //for(samples)
  



  C_region->Add(MT2[0], +1);
  C_region->Add(MT2[1], -1);
  C_region->Add(MT2[2], +1);

  B_region->Add(MT2[10], +1);
  B_region->Add(MT2[11], -1);
  B_region->Add(MT2[12], +1);

  //  TCanvas *con = new TCanvas ("con", "con");
//   con->Divide(2,2);
//   con->cd(1);
//   C_region->Draw();
//   con->cd(2);
//   B_region->Draw();

  //  Out_b->Divide(in);
  
//  AddOverAndUnderFlow(out, true, true);
//  AddOverAndUnderFlow(in, true, true);

//  AddOverAndUnderFlow(MT2[0], true, true);
//  AddOverAndUnderFlow(MT2[1], true, true);



  //  TFile *f = new TFile("hsimple.root");
 
  //  TH1F *hpx = (TH1F*)f->Get("hpx");
 
  //create a function with 3 parameters in the range [-3,3]
  //  TF1 *func = new TF1("fit",fitf,-3,3,3);
  //  func->SetParameters(500,hpx->GetMean(),hpx->GetRMS());
  //  func->SetParNames("Constant","Mean_value","Sigma");


//  double out1 = out->Integral();
//  double in1 = in->Integral();
//  cout << "sumMtgt250 / " << out1 / in1 << endl ;
//  MT2[2]->Divide(MT2[0],MT2[1]);

  //  MT2[3]->Multiply(MT2[1],MT2[2]);

  //  cout << "Wjets estimation:.." << MT2[3]->Integral() << endl;

//     double e = 0;
//     double f = 0;
//     double g = 0;
//     double h = 0;
//     double i = 0;
//     double j = 0;

//   for(int kk = 0; kk <= 1 ; kk++){
    
    
//     for(int mm = 0; mm < nbins+1 ; mm++){
      
//       e = in->GetBinContent(mm);      
//       f = out_b->GetBinContent(mm);

//       g = in->GetBinError(mm);
//       h = out_b->GetBinError(mm);
      
//       cout << cnames[kk]<<": " <<bins[mm] <<"-" << bins[mm+1] <<":.. " << e <<"__error:.."<< g << endl;
//       cout << cnames[kk]<<": " <<bins[mm] <<"-" << bins[mm+1] <<":.. " << f <<"__error:.."<< h << endl;

//     } 

//   }


//   double in_yield, out_b_yield, eff_yield = 0.0;
//   double in_err, out_b_err, eff_err = 0.0 ;

//   in_yield = in -> IntegralAndError(0, nbins, in_err) ;
//   out_b_yield = out_b -> IntegralAndError(0, nbins, out_b_err) ;

//   cout << cnames[0] << "-yield:.." << in_yield  <<"__error:..." << in_err  << endl ;
//   cout << cnames[1] << "-yield:.." << out_b_yield <<"__error:.."  << out_b_err << endl ;

//  out_a->Divide(in, out_b,1,1);

//     for(int mm = 0; mm < nbins+1 ; mm++){
      
//       i = out_a->GetBinContent(mm);
//       j = out_a->GetBinError(mm)  ;

//       cout << cnames[2]<<": " <<bins[mm] <<"-" << bins[mm+1] <<":.. " << i <<"__error:.."<< j << endl;

//     } 

//     eff_yield = out_a -> IntegralAndError(0, nbins, eff_err); 
//   cout << cnames[2] << "-yield:.." << eff_yield <<"__error:.."  << eff_err << endl ;

  C_region->Divide(B_region);  
  C_region->Fit("pol0");

  //  con->cd(3);
  C_region->Draw();


  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;
    
  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + myfileName +".root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();
  B_region->Write();
  C_region->Write();
  //  out_a->Write();

  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;


}

void MassPlotter::QCDEstimationinTauTau(TString cuts, TString trigger, unsigned int nevents, TString myfileName ){

TH1::SetDefaultSumw2();

int NumberOfBins =5;
//float xbin[] = {40,50,65,90,140};
float xbin[] = {100,125,150,200,250,300};



TH1D  *hdatainregionA  = new TH1D("hdatainregionA  ", "hdatainregionA ", NumberOfBins ,xbin);
TH1D  *hnonQCDinregionA  = new TH1D("hnonQCDinregionA",  "hnonQCDinregionA", NumberOfBins ,xbin);


TH1D  *hdatainregionB  = new TH1D("hdatainregionB",  "hdatainregionB", NumberOfBins ,xbin);
TH1D  *hnonQCDinregionB  = new TH1D("hnonQCDinregionB" , "hnonQCDinregionB",NumberOfBins ,xbin);


TH1D  *hdatainregionC  = new TH1D("hdatainregionC",  "hdatainregionC", NumberOfBins,xbin);
TH1D  *hnonQCDinregionC  = new TH1D("hnonQCDinregionC",  "hnonQCDinregionC", NumberOfBins,xbin);

TH1D  *hdatainregionQ1  = new TH1D("hdatainregionQ1",  "hdatainregionQ1", NumberOfBins,xbin);
TH1D  *hnonQCDinregionQ1  = new TH1D("hnonQCDinregionQ1",  "hnonQCDinregionQ1",NumberOfBins ,xbin);

TH1D  *hdatainregionQ2  = new TH1D("hdatainregionQ2",  "hdatainregionQ2",NumberOfBins ,xbin);
TH1D  *hnonQCDinregionQ2 = new TH1D("hnonQCDinregionQ2",  "hnonQCDinregionQ2", NumberOfBins,xbin);

TH1D  *hQCDinAllregion = new TH1D("hQCDinAllregion",  "hQCDinAllregion", NumberOfBins,xbin);
TH1D  *hQCDinregionC1 = new TH1D("hQCDinregionC1",  "hQCDinregionC1", NumberOfBins,xbin);
TH1D  *hQCDinregionC2 = new TH1D("hQCDinregionC2",  "hQCDinregionC2", NumberOfBins,xbin);
TH1D  *hQCDinregionC3 = new TH1D("hQCDinregionC3",  "hQCDinregionC3", NumberOfBins,xbin);


TH1D  *hQCDinregionQ3 = new TH1D("hQCDinregionQ3",  "hQCDinregionQ3", NumberOfBins,xbin);



 for(unsigned int ii = 0; ii < fSamples.size(); ii++){
    
   TString myCuts = cuts;
 
   int data = 0;
   sample Sample = fSamples[ii];
    
   if(Sample.type == "data"){
     data = 1;
     myCuts += " && " + trigger;
   }
    if(Sample.type == "susy")
    continue;


   fMT2tree = new MT2tree();
   Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

   float Weight;
    
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
        Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

    unsigned int nentries =  myEvtList->GetN();

    for (unsigned int jentry=0; jentry<min(nentries, nevents);jentry++) {
       
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

   
     if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){

	Weight = Sample.lumi;
	if(fPUReweight) Weight /= Sample.PU_avg_weight;
      
      }
 
      float weight = Weight;

      if(data == 1)
      weight = 1.0;
      else{

	
       weight *= fMT2tree->SFWeight.BTagCSV40eq0 ;

	if(fPUReweight)
	weight *= fMT2tree->pileUp.Weight;
	  }


     float myQuantity = fMT2tree->SumMTTauTauChannel();
       


   int TauIsolation= fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].Isolation3Hits+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].Isolation3Hits;

	 
             if(Sample.type == "data"  && (TauIsolation>=8) && fMT2tree->SumMTTauTauChannel()>250 && fMT2tree->doubleTau[0].GetSumCharge() !=0

		    ) {
         
          hdatainregionA  ->Fill (myQuantity,weight);
	  }

      if(Sample.type == "data"  && (TauIsolation >=8)&& fMT2tree->doubleTau[0].GetSumCharge() !=0
		   ){
          hdatainregionB  ->Fill (myQuantity,weight);
          }

       if(Sample.type == "data"  && TauIsolation<=6 &&(fMT2tree->doubleTau[0].GetSumCharge() ==0) ){
          
          hdatainregionC  ->Fill (myQuantity,weight);
	  }  


       if(Sample.type != "data"  && Sample.sname != "QCD"  &&(TauIsolation>=8)	&& fMT2tree->SumMTTauTauChannel()>250&& fMT2tree->doubleTau[0].GetSumCharge() !=0
		 ){
          
          hnonQCDinregionA->Fill (myQuantity,weight);
	   }

        if(Sample.type != "data"  && Sample.sname != "QCD" && (TauIsolation>=8)	 && fMT2tree->doubleTau[0].GetSumCharge() !=0  ){
         
          hnonQCDinregionB->Fill (myQuantity,weight);
	  }  

       if(Sample.type != "data"  && Sample.sname != "QCD" && (TauIsolation<=6 )&& (fMT2tree->doubleTau[0].GetSumCharge() ==0)){
         
          hnonQCDinregionC->Fill (myQuantity,weight);
	  }

  

        if(Sample.type == "data"  && fMT2tree->misc.MinMetJetDPhiPt40 > 1 && (TauIsolation >=8)&& fMT2tree->doubleTau[0].GetSumCharge() !=0
		   ){
          hdatainregionQ1  ->Fill (myQuantity,weight);

	  }  

        if(Sample.type != "data"  && Sample.sname != "QCD" && fMT2tree->misc.MinMetJetDPhiPt40 > 1 &&(TauIsolation >=8)&& fMT2tree->doubleTau[0].GetSumCharge() !=0)
		   {
          hnonQCDinregionQ1  ->Fill (myQuantity,weight);

	  }   
	
        if(Sample.type == "data"  && fMT2tree->misc.MinMetJetDPhiPt40 > 1 && (TauIsolation <=6)  && (fMT2tree->doubleTau[0].GetSumCharge() ==0)
		   ){
          hdatainregionQ2  ->Fill (myQuantity,weight);

	  }  

        if(Sample.type != "data"  && Sample.sname != "QCD" && fMT2tree->misc.MinMetJetDPhiPt40 > 1 &&(TauIsolation <=6) && (fMT2tree->doubleTau[0].GetSumCharge() ==0))
		   {
          hnonQCDinregionQ2  ->Fill (myQuantity,weight);

	   }
          }
        }

      
      cout<<"data"<<hdatainregionA->GetBinContent(5)<<" +- "<<hdatainregionA->GetBinError(5)<<endl;
      cout<<"hnonQCDinregionA"<<hnonQCDinregionA->GetBinContent(5)<<" +- "<<hnonQCDinregionA->GetBinError(5)<<endl;
      float errstatnonQCDA= hnonQCDinregionA->GetBinError(5);
      float nonQCDA= hnonQCDinregionA->GetBinContent(5);
      float errsysnonQCDA = nonQCDA*0.25;
      hnonQCDinregionA->SetBinError(5,sqrt(errsysnonQCDA * errsysnonQCDA + errstatnonQCDA * errstatnonQCDA));
      cout<<"hnonQCDinregionA"<<hnonQCDinregionA->GetBinContent(5)<<" +- "<<hnonQCDinregionA->GetBinError(5)<<endl;
      hdatainregionA->Add(hnonQCDinregionA, -1);
      TH1D *hQCDinregionA = (TH1D*)hdatainregionA->Clone();
      cout<<"hQCDinregionA"<<hQCDinregionA->GetBinContent(5)<<" +- "<<hQCDinregionA->GetBinError(5)<<endl;
      
      float nonQCDB=0;
    //  float err=0;
      float newerrB=0;
     float newerrB2=0;
   //  float errstatnonQCDB[NumberOfBins-1];
   //   float errsysnonQCDB[NumberOfBins-1];
      
      for(int j=0;j<NumberOfBins;j++){
      
      float errstatnonQCDB=hnonQCDinregionB->GetBinError(j);
   //  errstatnonQCDB[j-1]=hnonQCDinregionB->GetBinError(j);
      float errsysnonQCDB=0.25*(hnonQCDinregionB->GetBinContent(j));
 //    errsysnonQCDB[j-1]=0.25*(hnonQCDinregionB->GetBinContent(j));
       
      //cout<<"errsys"<<errsysnonQCDB[j-1]<<endl;
       cout<<"errsys"<<errsysnonQCDB<<endl;
  //    hnonQCDinregionB->SetBinError(j,sqrt(errsysnonQCDB[j-1] * errsysnonQCDB[j-1]+ errstatnonQCDB[j-1] * errstatnonQCDB[j-1]));
        //hnonQCDinregionB->SetBinError(j,sqrt(errsysnonQCDB * errsysnonQCDB+ errstatnonQCDB * errstatnonQCDB));
hnonQCDinregionB->SetBinError(j,sqrt(hnonQCDinregionB->GetBinError(j)*hnonQCDinregionB->GetBinError(j) + 0.25*(hnonQCDinregionB->GetBinContent(j)) * 0.25*(hnonQCDinregionB->GetBinContent(j))));


       nonQCDB += hnonQCDinregionB->GetBinContent(j);
       cout<<"nonQCDB"<<nonQCDB<<endl;
       newerrB+= hnonQCDinregionB->GetBinError(j);
       newerrB2+= (hnonQCDinregionB->GetBinError(j) *hnonQCDinregionB->GetBinError(j));
          }

       cout<<"nonQCDB"<<nonQCDB<<endl;
       cout<<"newerrB"<<newerrB<<endl;
      cout<<"newerrB2"<<sqrt(newerrB2)<<endl;
      
      hdatainregionB->Add(hnonQCDinregionB, -1);
      TH1D *hQCDinregionB = (TH1D*)hdatainregionB->Clone();
       double err;
      hQCDinregionB->IntegralAndError(1,4,err);
      cout<<"QCD in region B"<<hQCDinregionB->IntegralAndError(1,4,err)<<"+-"<<err<<endl;


     float nonQCDC=0;
     float newerrC=0;
     float newerrC2=0;
     for(int j=0;j<NumberOfBins;j++){
      
      float errstatnonQCDC=hnonQCDinregionC->GetBinError(j);
      float errsysnonQCDC=0.25*(hnonQCDinregionC->GetBinContent(j));

   hnonQCDinregionC->SetBinError(j,sqrt(hnonQCDinregionC->GetBinError(j)*hnonQCDinregionC->GetBinError(j) + 0.25*(hnonQCDinregionC->GetBinContent(j)) * 0.25*(hnonQCDinregionC->GetBinContent(j))));


       nonQCDC += hnonQCDinregionC->GetBinContent(j);
       cout<<"nonQCDC"<<nonQCDC<<endl;
       newerrC+= hnonQCDinregionC->GetBinError(j);
       newerrC+= (hnonQCDinregionC->GetBinError(j) *hnonQCDinregionC->GetBinError(j));
          }

       cout<<"nonQCDC"<<nonQCDC<<endl;
       cout<<"newerrC"<<newerrC<<endl;
      cout<<"newerrC2"<<sqrt(newerrC2)<<endl;




      hdatainregionC->Add(hnonQCDinregionC, -1);
      TH1D *hQCDinregionC = (TH1D*)hdatainregionC->Clone();

hQCDinregionC->IntegralAndError(1,4,err);
   cout<<"QCD in region C"<<hQCDinregionC->IntegralAndError(1,4,err)<<"+-"<<err<<endl;
      
       hQCDinregionC1->Divide(hQCDinregionC,hQCDinregionB);
    
      
 /*  TF1 *func = new TF1("func","pol1");
      func->SetRange(100,250);
      hQCDinregionC1->Fit("func","R");
      double b0= func->GetParameter(0);
      double Deltab0=func->GetParError(0);
      double b1= func->GetParameter(1);
      double Deltab1=func->GetParError(1);
      cout<<">>>>"<<b0<<"+-"<<Deltab0<<endl;
      cout<<">>>>"<<b1<<"+-"<<Deltab1<<endl;

      for (int j=0;j<NumberOfBins+1;j++){
      hQCDinregionC3->SetBinContent(j,func->GetParameter(0));
      hQCDinregionC3->SetBinError(j,func->GetParError(0));
      hQCDinregionC2->SetBinContent(j,func->GetParameter(1));
      hQCDinregionC2->SetBinError(j,func->GetParError(1));
        }
       hQCDinregionC2->Scale(250);
          hQCDinregionC2->Add(hQCDinregionC3);

     
      cout<<"<<<<<<<<<"<<hQCDinregionC2->GetBinContent(3)<<" +- "<<hQCDinregionC2->GetBinError(3)<<endl;*/

   TF1 *func = new TF1("func","pol0");
      func->SetRange(100,250);
      hQCDinregionC1->Fit("func","R");
  for (int j=0;j<NumberOfBins+1;j++){
     hQCDinregionC1->SetBinContent(j,func->GetParameter(0));
      hQCDinregionC1->SetBinError(j,func->GetParError(0));}
      cout<<"<<<<<<<<<"<<hQCDinregionC1->GetBinContent(3)<<" +- "<<hQCDinregionC1->GetBinError(3)<<endl;

 // for(int j=0;j<NumberOfBins+1;j++){
  //  hQCDinregionC1->SetBinContent(j,hQCDinregionC1->GetBinContent(3));
 //   hQCDinregionC1->SetBinError(j,hQCDinregionC1->GetBinError(3));
 // }
 
       cout<<"TransferFactor"<<hQCDinregionC1->GetBinContent(1)<<" +- "<<hQCDinregionC1->GetBinError(1)<<endl;
      cout<<"TransferFactor"<<hQCDinregionC1->GetBinContent(2)<<" +- "<<hQCDinregionC1->GetBinError(2)<<endl;
      cout<<"TransferFactor"<<hQCDinregionC1->GetBinContent(3)<<" +- "<<hQCDinregionC1->GetBinError(3)<<endl;
      cout<<"TransferFactor"<<hQCDinregionC1->GetBinContent(4)<<" +- "<<hQCDinregionC1->GetBinError(4)<<endl;
  //      cout<<"TransferFactor"<<hQCDinregionC1->GetBinContent(5)<<" +- "<<hQCDinregionC1->GetBinError(5)<<endl;

  
             hdatainregionQ1->Add(hnonQCDinregionQ1, -1);
             hdatainregionQ2->Add(hnonQCDinregionQ2, -1);

             TH1D *hQCDinregionQ1=(TH1D*) hdatainregionQ1->Clone();
             TH1D *hQCDinregionQ2=(TH1D*) hdatainregionQ2->Clone();
    
              hQCDinregionQ1->Add(hQCDinregionQ2);
              TH1D *hQCDinregionB1 = (TH1D*)hQCDinregionB->Clone();
       

   


              hQCDinregionB->Add(hQCDinregionC);

               
              hQCDinregionQ1->Divide(hQCDinregionB);
        
 /*         
      TF1 *func2 = new TF1("func2","pol1");
      func2->SetRange(40,90);
      hQCDinregionQ1->Fit("func2","R");
      double a0= func2->GetParameter(0);
      double Deltaa0=func2->GetParError(0);
      double a1= func2->GetParameter(1);
      double Deltaa1=func2->GetParError(1);
      cout<<">>>>"<<a0<<"+-"<<Deltaa0<<endl;
      cout<<">>>>"<<a1<<"+-"<<Deltaa1<<endl;
      
      for (int j=0;j<NumberOfBins+1;j++){
      hQCDinregionQ3->SetBinContent(j,func2->GetParameter(0));
      hQCDinregionQ3->SetBinError(j,func2->GetParError(0));
      hQCDinregionQ2->SetBinContent(j,func2->GetParameter(1));
      hQCDinregionQ2->SetBinError(j,func2->GetParError(1));
}
    
     hQCDinregionQ2->Scale(90);
     hQCDinregionQ2->Add(hQCDinregionQ3);

    */  
     // double err;
      //cout<<">>>>>>"<<hQCDinregionQ2->IntegralAndError(8,9,err);
      
  //    cout<<"<<<<<<<<<"<<hQCDinregionQ2->GetBinContent(6)<<" +- "<<hQCDinregionQ2->GetBinError(6)<<endl;

  /*  TF1 *func2 = new TF1("func2","pol0");
      func2->SetRange(40,90);
     hQCDinregionQ1->Fit("func2","R");
     for (int j=0;j<NumberOfBins+1;j++){
    hQCDinregionQ1->SetBinContent(j,func2->GetParameter(0));
      hQCDinregionQ1->SetBinError(j,func2->GetParError(0));}
      cout<<"<<<<<<<<<"<<hQCDinregionQ1->GetBinContent(3)<<" +- "<<hQCDinregionQ1->GetBinError(3)<<endl;*/





    
       
for(int j=0;j<NumberOfBins+1;j++){
  hQCDinregionQ1->SetBinContent(j,hQCDinregionQ1->GetBinContent(j));
  hQCDinregionQ1->SetBinError(j,hQCDinregionQ1->GetBinError(j));
}
         
   
cout<<"EfficiencyQ1 "<<hQCDinregionQ1->GetBinContent(1)<<" +- "<<hQCDinregionQ1->GetBinError(1)<<endl; 
cout<<"Efficiency"<<hQCDinregionQ1->GetBinContent(2)<<" +- "<<hQCDinregionQ1->GetBinError(2)<<endl; 
cout<<"Efficiency"<<hQCDinregionQ1->GetBinContent(3)<<" +- "<<hQCDinregionQ1->GetBinError(3)<<endl;       
cout<<"Efficiency"<<hQCDinregionQ1->GetBinContent(4)<<" +- "<<hQCDinregionQ1->GetBinError(4)<<endl;


for (int j=0;j<NumberOfBins+1;j++){
hQCDinregionA->SetBinContent(j,hQCDinregionA->GetBinContent(5));
hQCDinregionA->SetBinError(j,hQCDinregionA->GetBinError(5));
}

   
cout<<"QCDinregionA"<<hQCDinregionA->GetBinContent(5)<<" +- "<<hQCDinregionA->GetBinError(5)<<endl; 
 
         
      hQCDinAllregion-> Multiply(hQCDinregionB1,hQCDinregionC1);
      hQCDinAllregion-> Multiply(hQCDinregionQ1);
     
     cout<<"QCDinAllregion"<<hQCDinAllregion->GetBinContent(1)<<" +- "<<hQCDinAllregion->GetBinError(1)<<endl; 
 cout<<"QCDinAllregion"<<hQCDinAllregion->GetBinContent(2)<<" +- "<<hQCDinAllregion->GetBinError(2)<<endl; 
 cout<<"QCDinAllregion"<<hQCDinAllregion->GetBinContent(3)<<" +- "<<hQCDinAllregion->GetBinError(3)<<endl; 
 cout<<"QCDinAllregion"<<hQCDinAllregion->GetBinContent(4)<<" +- "<<hQCDinAllregion->GetBinError(4)<<endl;

  
  

     hQCDinregionA->Multiply(hQCDinregionC1);
     hQCDinregionA->Multiply(hQCDinregionQ1);



         
       cout<<"QCDinregionD"<<hQCDinregionA->GetBinContent(4)<<" +- "<<hQCDinregionA->GetBinError(4)<<endl; 
 
 
       TCanvas *MyC = new TCanvas("QCDEstimation","QCDEstimation");
        MyC->Divide(1,1);
        MyC->cd(1);
        hQCDinregionQ1->GetXaxis()->SetTitle("M_{T2}");
        hQCDinregionQ1->GetYaxis()->SetTitle("#Delta #Phi_{4}^{min}");
        hQCDinregionQ1->Draw();
  

     
  TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  Util::MakeOutputDir(fileName);
  fileName = fileName + myfileName +"_Histos.root";
  TFile *savefile = new TFile(fileName.Data(), "RECREATE");
  savefile ->cd();  
  hQCDinregionB1->Write();	  
  hQCDinregionC1->Write();	
  hQCDinregionQ1->Write();
  hQCDinAllregion->Write();
  savefile->Close();
  std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
  
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

	  }  

void MassPlotter::TauTauAnalysisbin2(TString cuts,TString trigger,unsigned int nevents, TString myfileName){
  
  TH1::SetDefaultSumw2();
  
  setFlags(10);

  /*TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  fileName = fileName + myfileName;
  TFile *myfile1 = new TFile(fileName.Data(), "READ");

  //TH1D* hdatainregionB  = (TH1D*) myfile1->Get("hdatainregionB");
//  TH1D* hdatainregionC  = (TH1D*) myfile1->Get("hdatainregionC");
  TH1D* hQCDinAllregion = (TH1D*) myfile1->Get("hQCDinAllregion");*/
  //  TFile *myfile1 = new TFile("../MassPlots/MT2_TauTau_Histos.root", "READ");
TFile *myfile1 = new TFile("../MassPlots/MT2_TauTau_QCDestimation_Bin2_Histos.root", "READ");
TH1D* hQCDinAllregion = (TH1D*) myfile1->Get("hQCDinAllregion");

   static const int nbins =5;
   //double bins[nbins+1] ={40,50,65,90,140};
double bins[nbins+1] ={100,125,150,200,250,1000};
//{"QCD", "W", "ZX", "Top", "WW", "Higgs", "MC", "susy","data"};
//int ccolor[NumberOfSamples+1] = { 401, 417, 419, 855, 603, kRed, 603, 1, 632};


  TString  cnames[NumberOfSamples+1] = {"Higgs",  "WW", "Top", "ZX","W", "QCD","MC", "susy", "data"};
  int      ccolor[NumberOfSamples+1] = { kRed,    603,    855, 419 , 417, 401, 500,     1, 632};
   
   TString varname = "MT2";
  for (int i=0; i<=(NumberOfSamples+1); i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i], "", nbins, bins);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

   MT2[NumberOfSamples] -> SetMarkerStyle(20);
   MT2[NumberOfSamples] -> SetMarkerColor(kBlack);
   MT2[NumberOfSamples] -> SetLineColor(kBlack);
   MT2[NumberOfSamples-2] -> SetFillStyle(3004);
   MT2[NumberOfSamples-2] -> SetFillColor(kBlack);

cout<<" trigger "<<trigger<<endl;
cout<<" cuts "<<cuts<<endl;
 

 
  
  vector<TH1D*> h_samples;

   

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){

    TString myCuts = cuts;

    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data"){
      data = 1;
      myCuts += " && " + trigger;
    }

   
 
    h_samples.push_back(new TH1D(varname+"_"+Sample.name, "", nbins, bins));
    h_samples[ii] -> Sumw2();
    h_samples[ii] -> SetLineColor(Sample.color);
   
    h_samples[ii] -> SetMarkerColor(Sample.color);
    h_samples[ii] -> SetStats(false);
    if(Sample.type == "susy" ){
    h_samples[ii] -> SetLineColor(kBlack);
    h_samples[ii] -> SetLineStyle(kDotted);
    }

  
    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);


    float Weight;
   
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "   looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

int TauIsolation= fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].Isolation3Hits+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].Isolation3Hits;    

 
      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
   //   Sample.tree->SetEventList(myEvtList);

      unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();

      for (unsigned int jentry=0; jentry<nentries;jentry++) {
    
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

	if ( fVerbose>2 && jentry % 100000 == 0 ){
	  fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	  fflush(stdout);
	}

	 if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){
	  
	  Weight = Sample.lumi;
	  if(fPUReweight) Weight /= Sample.PU_avg_weight;
	  
	}

        float weight = Weight;
	if(data == 1)
	weight = 1.0;
        else{
	
	
	weight *= fMT2tree->SFWeight.BTagCSV40eq0;  

	if(fPUReweight)
	  weight *= fMT2tree->pileUp.Weight;  }   
        
        float myQuantity=fMT2tree->SumMTTauTauChannel();
 
          h_samples[ii]->Fill(myQuantity, weight);
             
      } 
        AddOverAndUnderFlow(h_samples[ii]);



        
	if(data == 1){
        MT2[NumberOfSamples]->Add(h_samples[ii]);//data
}
else if(Sample.sname == "SUSY")
MT2[NumberOfSamples-1]->Add(h_samples[ii]);

else if(Sample.sname == "Top")
MT2[2]->Add(h_samples[ii]);
else
if(Sample.sname == "DY")
MT2[3]->Add(h_samples[ii]);
else
if(Sample.sname == "Wtolnu"){
MT2[4]->Add(h_samples[ii]);}

else
if(Sample.sname == "VV")
MT2[1]->Add(h_samples[ii]);
else
if(Sample.sname == "Higgs")
MT2[0]->Add(h_samples[ii]);

else if (Sample.sname  == "QCD")    {    
 for (int j=0;j<nbins+1;j++){ 
MT2[5]->SetBinContent(1,hQCDinAllregion->GetBinContent(1));
MT2[5]->SetBinError(1,hQCDinAllregion->GetBinError(1)); 
MT2[5]->SetBinContent(2,hQCDinAllregion->GetBinContent(2));
MT2[5]->SetBinError(2,hQCDinAllregion->GetBinError(2));
MT2[5]->SetBinContent(3,hQCDinAllregion->GetBinContent(3));
MT2[5]->SetBinError(3,hQCDinAllregion->GetBinError(3));
MT2[5]->SetBinContent(4,hQCDinAllregion->GetBinContent(4) );
MT2[5]->SetBinError(4,hQCDinAllregion->GetBinError(4));

MT2[5]->SetBinContent(5,0.82);
MT2[5]->SetBinError(5,0.72);
                               
                  }          
                 
             }

  }
  Cout(0,MT2[NumberOfSamples-1]);
	Cout(0,MT2[0]);    


MT2[NumberOfSamples-2]->Add(MT2[0]);
MT2[NumberOfSamples-2]->Add(MT2[1]);
MT2[NumberOfSamples-2]->Add(MT2[2]);
MT2[NumberOfSamples-2]->Add(MT2[3]);
MT2[NumberOfSamples-2]->Add(MT2[4]);
MT2[NumberOfSamples-2]->Add(MT2[5]);

        cout<<">>>>>>>>>>>"<<hQCDinAllregion->GetBinContent(1)<<" +- "<<hQCDinAllregion->GetBinError(1)<<endl;
        cout<<">>>>>>>>>>>"<<hQCDinAllregion->GetBinContent(2)<<" +- "<<hQCDinAllregion->GetBinError(2)<<endl;
        cout<<">>>>>>>>>>>"<<hQCDinAllregion->GetBinContent(3)<<" +- "<<hQCDinAllregion->GetBinError(3)<<endl;
       cout<<">>>>>>>>>>>"<<hQCDinAllregion->GetBinContent(4)<<" +- "<<hQCDinAllregion->GetBinError(4)<<endl;
     

  
        cout<<">>>>>>>>>>>"<<MT2[5]->GetBinContent(1)<<" +- "<<MT2[5]->GetBinError(1)<<endl;
        cout<<">>>>>>>>>>>"<<MT2[5]->GetBinContent(2)<<" +- "<<MT2[5]->GetBinError(2)<<endl;
        cout<<">>>>>>>>>>>"<<MT2[5]->GetBinContent(3)<<" +- "<<MT2[5]->GetBinError(3)<<endl;
        cout<<">>>>>>>>>>>"<<MT2[5]->GetBinContent(4)<<" +- "<<MT2[5]->GetBinError(4)<<endl;
   
  Cout(1,MT2[NumberOfSamples-1]);
	Cout(1,MT2[0]);    


AddOverAndUnderFlow(MT2[0], true, true);
AddOverAndUnderFlow(MT2[1], true, true);
AddOverAndUnderFlow(MT2[2], true, true);
AddOverAndUnderFlow(MT2[3], true, true);
AddOverAndUnderFlow(MT2[4], true, true);
AddOverAndUnderFlow(MT2[5], true, true);
AddOverAndUnderFlow(MT2[NumberOfSamples-1], true, true);
AddOverAndUnderFlow(MT2[NumberOfSamples-2], true, true);
AddOverAndUnderFlow(MT2[NumberOfSamples], true, true);
  
 Cout(2,MT2[NumberOfSamples-1]);
	Cout(2,MT2[0]);    


  printYield();

  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j <= (NumberOfSamples+1); j++){
  if(j <= (NumberOfSamples - 3))
    h_stack -> Add(MT2[j]);
      
  }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
Legend1->AddEntry(MT2[0], "Higgs", "f");
Legend1->AddEntry(MT2[1], "WW", "f");
Legend1->AddEntry(MT2[2], "Top", "f");
Legend1->AddEntry(MT2[3], "ZX", "f");
Legend1->AddEntry(MT2[4], "W", "f");
Legend1->AddEntry(MT2[5], "QCD", "f");
Legend1->AddEntry(MT2[NumberOfSamples-1], "SMS", "l");
Legend1->AddEntry(MT2[NumberOfSamples], "data", "l");
 
  
 
  
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

 
  printHisto(h_stack, MT2[NumberOfSamples], MT2[NumberOfSamples-2], MT2[NumberOfSamples-1], Legend1 , "MTC", "hist", true, "#Sigma m_{T}^{#tau_{i}}", "Events", 0, -10, 2, true);
plotRatioStack(h_stack, MT2[NumberOfSamples-2], MT2[NumberOfSamples], MT2[NumberOfSamples-1], true, false, "MT2_ratio", Legend1, "#Sigma m_{T}^{#tau_{i}}", "Events", 0, -10, 2, true, "C");

}

void MassPlotter::TauTauAnalysisbin1(TString cuts,TString trigger,unsigned int nevents, TString myfileName){
  
  TH1::SetDefaultSumw2();
  
  setFlags(10);

  /*TString fileName = fOutputDir;
  if(!fileName.EndsWith("/")) fileName += "/";
  fileName = fileName + myfileName;
  TFile *myfile1 = new TFile(fileName.Data(), "READ");

  //TH1D* hdatainregionB  = (TH1D*) myfile1->Get("hdatainregionB");
//  TH1D* hdatainregionC  = (TH1D*) myfile1->Get("hdatainregionC");
  TH1D* hQCDinAllregion = (TH1D*) myfile1->Get("hQCDinAllregion");*/
  //  TFile *myfile1 = new TFile("../MassPlots/MT2_TauTau_Histos.root", "READ");
TFile *myfile1 = new TFile("../MassPlots/MT2_TauTau_QCDestimation_Bin1_Histos.root", "READ");
TH1D* hQCDinAllregion = (TH1D*) myfile1->Get("hQCDinAllregion");

   static const int nbins =4;
   double bins[nbins+1] ={40,50,65,90,140};
//double bins[nbins+1] ={100,125,150,200,250,300};
//{"QCD", "W", "ZX", "Top", "WW", "Higgs", "MC", "susy","data"};
//int ccolor[NumberOfSamples+1] = { 401, 417, 419, 855, 603, kRed, 603, 1, 632};


  TString  cnames[NumberOfSamples+1] = {"Higgs",  "WW", "Top", "ZX","W", "QCD","MC", "susy", "data"};
  int      ccolor[NumberOfSamples+1] = { kRed,    603,    855, 419 , 417, 401, 500,     1, 632};
   
   TString varname = "MT2";
  for (int i=0; i<=(NumberOfSamples); i++){
    MT2[i] = new TH1D(varname+"_"+cnames[i], "", nbins, bins);
    MT2[i] -> SetFillColor (ccolor[i]);
    MT2[i] -> SetLineColor (ccolor[i]);
    MT2[i] -> SetLineWidth (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

   MT2[NumberOfSamples] -> SetMarkerStyle(20);
   MT2[NumberOfSamples] -> SetMarkerColor(kBlack);
   MT2[NumberOfSamples] -> SetLineColor(kBlack);
   MT2[NumberOfSamples-2] -> SetFillStyle(3004);
   MT2[NumberOfSamples-2] -> SetFillColor(kBlack);

cout<<" trigger "<<trigger<<endl;
cout<<" cuts "<<cuts<<endl;
 

 
  
  vector<TH1D*> h_samples;

   

  for(unsigned int ii = 0; ii < fSamples.size(); ii++){

    TString myCuts = cuts;

    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data"){
      data = 1;
      myCuts += " && " + trigger;
    }

   
 
    h_samples.push_back(new TH1D(varname+"_"+Sample.name, "", nbins, bins));
    h_samples[ii] -> Sumw2();
    h_samples[ii] -> SetLineColor(Sample.color);
   
    h_samples[ii] -> SetMarkerColor(Sample.color);
    h_samples[ii] -> SetStats(false);
    if(Sample.type == "susy" ){
    h_samples[ii] -> SetLineColor(kBlack);
    h_samples[ii] -> SetLineStyle(kDotted);
    }

  
    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);


    float Weight;
   
    if(fPUReweight) Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents*Sample.PU_avg_weight);
    else            Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "   looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

int TauIsolation= fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex0()].Isolation3Hits+fMT2tree->tau[fMT2tree->doubleTau[0].GetTauIndex1()].Isolation3Hits;    

 
      Sample.tree->Draw(">>selList", myCuts);
      TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
   //   Sample.tree->SetEventList(myEvtList);

      unsigned int nentries = myEvtList->GetN();//Sample.tree->GetEntries();

      for (unsigned int jentry=0; jentry<nentries;jentry++) {
    
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

	if ( fVerbose>2 && jentry % 100000 == 0 ){
	  fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	  fflush(stdout);
	}

	 if(fStitching && (Sample.sname == "Wtolnu" || (Sample.shapename == "ZJetsToLL" && Sample.name != "DYToLL_M10To50"))){
	  
	  Weight = Sample.lumi;
	  if(fPUReweight) Weight /= Sample.PU_avg_weight;
	  
	}

        float weight = Weight;
	if(data == 1)
	weight = 1.0;
        else{
	
	
	//weight *= fMT2tree->SFWeight.BTagCSV40eq0;  

	if(fPUReweight)
	  weight *= fMT2tree->pileUp.Weight;  }   
        
        float myQuantity=fMT2tree->doubleTau[0].GetMT2();
 
          h_samples[ii]->Fill(myQuantity, weight);
             
         } 
        AddOverAndUnderFlow(h_samples[ii]);
        
	if(data == 1){
        MT2[NumberOfSamples]->Add(h_samples[ii]);//data
}
else if(Sample.sname == "SUSY")
MT2[NumberOfSamples-1]->Add(h_samples[ii]);

else if(Sample.sname == "Top")
MT2[2]->Add(h_samples[ii]);
else
if(Sample.sname == "DY")
MT2[3]->Add(h_samples[ii]);
else
if(Sample.sname == "Wtolnu"){
MT2[4]->Add(h_samples[ii]);}

else
if(Sample.sname == "VV")
MT2[1]->Add(h_samples[ii]);
else
if(Sample.sname == "Higgs")
MT2[0]->Add(h_samples[ii]);

else if (Sample.sname  == "QCD")    {    
 for (int j=0;j<nbins+1;j++){ 
MT2[5]->SetBinContent(1,hQCDinAllregion->GetBinContent(1));
MT2[5]->SetBinError(1,hQCDinAllregion->GetBinError(1)); 
MT2[5]->SetBinContent(2,hQCDinAllregion->GetBinContent(2));
MT2[5]->SetBinError(2,hQCDinAllregion->GetBinError(2));
MT2[5]->SetBinContent(3,hQCDinAllregion->GetBinContent(3));
MT2[5]->SetBinError(3,hQCDinAllregion->GetBinError(3));
MT2[5]->SetBinContent(4,0.15 );
MT2[5]->SetBinError(4,0.35);

//MT2[0]->SetBinContent(5,0.82);
//MT2[0]->SetBinError(5,0.72);
                               
                  }          
                 
             }

  }
MT2[NumberOfSamples-2]->Add(MT2[0]);
MT2[NumberOfSamples-2]->Add(MT2[1]);
MT2[NumberOfSamples-2]->Add(MT2[2]);
MT2[NumberOfSamples-2]->Add(MT2[3]);
MT2[NumberOfSamples-2]->Add(MT2[4]);
MT2[NumberOfSamples-2]->Add(MT2[5]);

        cout<<">>>>>>>>>>>"<<hQCDinAllregion->GetBinContent(1)<<" +- "<<hQCDinAllregion->GetBinError(1)<<endl;
        cout<<">>>>>>>>>>>"<<hQCDinAllregion->GetBinContent(2)<<" +- "<<hQCDinAllregion->GetBinError(2)<<endl;
        cout<<">>>>>>>>>>>"<<hQCDinAllregion->GetBinContent(3)<<" +- "<<hQCDinAllregion->GetBinError(3)<<endl;
       cout<<">>>>>>>>>>>"<<hQCDinAllregion->GetBinContent(4)<<" +- "<<hQCDinAllregion->GetBinError(4)<<endl;
     

  
        cout<<">>>>>>>>>>>"<<MT2[5]->GetBinContent(1)<<" +- "<<MT2[5]->GetBinError(1)<<endl;
        cout<<">>>>>>>>>>>"<<MT2[5]->GetBinContent(2)<<" +- "<<MT2[5]->GetBinError(2)<<endl;
        cout<<">>>>>>>>>>>"<<MT2[5]->GetBinContent(3)<<" +- "<<MT2[5]->GetBinError(3)<<endl;
        cout<<">>>>>>>>>>>"<<MT2[5]->GetBinContent(4)<<" +- "<<MT2[5]->GetBinError(4)<<endl;
   
       

for(int j = 0; j <= (NumberOfSamples); j++){
    AddOverAndUnderFlow(MT2[j], true, true);
  }


  printYield();

  THStack* h_stack = new THStack(varname, "");
  for(int j = 0; j <= (NumberOfSamples); j++){
  if(j <= (NumberOfSamples - 3))
    h_stack -> Add(MT2[j]);
      
  }

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
Legend1->AddEntry(MT2[0], "Higgs", "f");
Legend1->AddEntry(MT2[1], "WW", "f");
Legend1->AddEntry(MT2[2], "Top", "f");
Legend1->AddEntry(MT2[3], "ZX", "f");
Legend1->AddEntry(MT2[4], "W", "f");
Legend1->AddEntry(MT2[5], "QCD", "f");
Legend1->AddEntry(MT2[NumberOfSamples-1], "SMS", "l");
Legend1->AddEntry(MT2[NumberOfSamples], "data", "l");
 
  
 
  
  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

 
  printHisto(h_stack, MT2[NumberOfSamples], MT2[NumberOfSamples-2], MT2[NumberOfSamples-1], Legend1 , "MTC", "hist", true, "M_{T2}", "Events", 0, -10, 2, true);
plotRatioStack(h_stack, MT2[NumberOfSamples-2], MT2[NumberOfSamples], MT2[NumberOfSamples-1], true, false, "MT2_ratio", Legend1, "M_{T2}", "Events", 0, -10, 2, true, "C");

}
