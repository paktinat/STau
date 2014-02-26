#include "../include/LHEFReader.h"

#include "TFile.h"
#include "TTree.h"


void LHEFReader::SetReader(string fileName){
  if(inputFileStream.is_open()){
    delete reader;
    inputFileStream.close();
  }
  inputFileStream.open( fileName.c_str() );
  reader = new LHEF::Reader( inputFileStream );
}

LHEFReader::LHEFReader()
{

}

bool LHEFReader::LoadNextEvent(){
  return reader->readEvent() ;
}

STauSTauEvent* LHEFReader::GetSTauSTauEvent(){

  const LHEF::HEPEUP &hepeup = reader->hepeup;

  int particle;

  int SLP = -1;
  int TauP = -1;
  int LSPP = -1;

  int SLN = -1;
  int TauN = -1;
  int LSPN = -1;

  for(particle = 0; particle < hepeup.NUP; ++particle)
  {
    if( hepeup.IDUP[particle] == 1000015 )
      SLP = particle;
    else if(hepeup.IDUP[particle] == -1000015 )
      SLN = particle;
    else if( abs(hepeup.IDUP[particle]) == 15 ){
      if(hepeup.MOTHUP[particle].first == SLP+1)
	TauP = particle;
      else if(hepeup.MOTHUP[particle].first == SLN+1)
	TauN = particle;
      else
	cout << "Orphan SM Particle" << endl;
    }else if( abs(hepeup.IDUP[particle]) == 1000022){
      if(hepeup.MOTHUP[particle].first == SLP+1)
	LSPP = particle;
      else if(hepeup.MOTHUP[particle].first == SLN+1)
	LSPN = particle;
      else
	cout << "Orphan LSP" << endl;
    }
  }

  particle = SLP ;
  STauSTau_Event.STauP.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
			    hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = SLN ;
  STauSTau_Event.STauM.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
			    hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = TauP ;
  STauSTau_Event.STauP.SMChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
				    hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = TauN ;
  STauSTau_Event.STauM.SMChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
				    hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = LSPN ;
  STauSTau_Event.STauM.SusyChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
				      hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = LSPP ;
  STauSTau_Event.STauP.SusyChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
				      hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );


  STauSTau_Event.CalcMT2();
  
  STauSTau_Event.STauMass = STauSTau_Event.STauP.mass;
  STauSTau_Event.LSPMass = STauSTau_Event.STauP.SusyChild.mass;

  return &STauSTau_Event;
}

CharginoChargino* LHEFReader::GetCharginoCharginoEvent(){

  const LHEF::HEPEUP &hepeup = reader->hepeup;

  int particle;

  int charginoP = -1;
  int SLP = -1;
  int SMP1 = -1;
  int SMP2 = -1;
  int LSPP = -1;

  int charginoN = -1;
  int SLN = -1;
  int SMN1 = -1;
  int SMN2 = -1;
  int LSPN = -1;

  for(particle = 0; particle < hepeup.NUP; ++particle)
  {
    //cout << particle << "," << charginoP << "," << SLP << "," << SMP1 << "," << SMP2 << "," << LSPP << " : " << charginoN << "," << SLN << "," << SMN1 << "," << SMN2 << "," << "," << LSPN << ";"  << endl;
    if( hepeup.IDUP[particle] == 1000024 )
      charginoP = particle;
    else if( hepeup.IDUP[particle] == -1000024 )
      charginoN = particle;
    else if( abs(hepeup.IDUP[particle]) == 1000016 || abs(hepeup.IDUP[particle]) == 1000015 ){
      if(hepeup.MOTHUP[particle].first == charginoN+1)
	SLN = particle;
      else if(hepeup.MOTHUP[particle].first == charginoP+1)
	SLP = particle;
      else
	cout << "Orphan SLeptong" << hepeup.IDUP[particle] << " " << hepeup.MOTHUP[particle].first << endl;
    }else if( abs(hepeup.IDUP[particle]) == 16 || abs(hepeup.IDUP[particle]) == 15 ){
      if(hepeup.MOTHUP[particle].first == charginoN+1)
	SMN1 = particle;
      else if(hepeup.MOTHUP[particle].first == charginoP+1)
	SMP1 = particle;
      else if(hepeup.MOTHUP[particle].first == SLP+1)
	SMP2 = particle;
      else if(hepeup.MOTHUP[particle].first == SLN+1)
	SMN2 = particle;
      else
	cout << "Orphan SM Particle" << endl;
    }else if( abs(hepeup.IDUP[particle]) == 1000022){
      if(hepeup.MOTHUP[particle].first == SLP+1)
	LSPP = particle;
      else if(hepeup.MOTHUP[particle].first == SLN+1)
	LSPN = particle;
      else
	cout << "Orphan LSP" << endl;
    }
  }


  particle = charginoP;
  CharginoCharginoEvent.CharginoP.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
				       hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = charginoN;
  CharginoCharginoEvent.CharginoN.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
				       hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = SLN ;
  CharginoCharginoEvent.CharginoN.SusyChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
						 hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = SLP ;
  CharginoCharginoEvent.CharginoP.SusyChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
						 hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = SMN1 ;
  CharginoCharginoEvent.CharginoN.SMChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
					       hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = SMP1 ;
  CharginoCharginoEvent.CharginoP.SMChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
					       hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = LSPN ;
  CharginoCharginoEvent.CharginoN.SusyChild.SusyChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
							   hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = LSPP ;
  CharginoCharginoEvent.CharginoP.SusyChild.SusyChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
							   hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = SMN2 ;
  CharginoCharginoEvent.CharginoN.SusyChild.SMChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
							   hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );

  particle = SMP2 ;
  CharginoCharginoEvent.CharginoP.SusyChild.SMChild.Set( hepeup.PUP[particle][0] , hepeup.PUP[particle][1] , hepeup.PUP[particle][2] , hepeup.PUP[particle][3] , hepeup.PUP[particle][4] ,
							   hepeup.ISTUP[particle]  , hepeup.MOTHUP[particle].first , hepeup.MOTHUP[particle].second , hepeup.IDUP[particle] );


  //CharginoCharginoEvent.CalcMET();
  CharginoCharginoEvent.CalcMT2();
  CharginoCharginoEvent.CalcDecayMode();
  
  CharginoCharginoEvent.CharginoMass = CharginoCharginoEvent.CharginoP.mass;
  CharginoCharginoEvent.LSPMass = CharginoCharginoEvent.CharginoP.SusyChild.SusyChild.mass;
  CharginoCharginoEvent.STauMass = CharginoCharginoEvent.CharginoP.SusyChild.mass;
  
  return &CharginoCharginoEvent;
}


int main(int argc, char *argv[])
{
  if(argc < 4)
    {
      cout << " Usage: " << argv[0] << "EventType output_file" << " input_files" << endl;
      cout << "EventType : 1 for CharginoChargino, 2 for STauSTau" << endl;
      cout << " output_file - output file in ROOT format." << endl;
      cout << " input_files - input files in LHEF format," << endl;
      return 1;
    }

  int event_type = atoi( argv[1] );
  string event_type_name = (event_type==1?"CharginoChargino":"STauSTau");
  cout << "Event Type : " << event_type_name << endl;

  LHEFReader theReader;

  TFile* rootFile = TFile::Open( TString(argv[2]) , "recreate");
  TTree tree("LHEFiles" , "");
  if(event_type == 1)
    tree.Branch("CharginoCharginoEvent" , &(theReader.CharginoCharginoEvent ) );
  else if(event_type == 2 )
    tree.Branch("STauSTauEvent" , &(theReader.STauSTau_Event) );


  for(int input_index = 3 ; input_index < argc ; input_index ++){

    theReader.SetReader( argv[input_index] );

    cout << "** Calculating number of events to process. Please wait..." << endl;
    Long64_t allEntries = theReader.reader->getNumberOfEvents();
    cout << "** Input file contains " << allEntries << " events" << endl;

    int nevent = 0;

    while( theReader.LoadNextEvent() ){
      //cout << nevent << endl;
      if(event_type == 1)
	theReader.GetCharginoCharginoEvent();
      else if(event_type == 2)
	theReader.GetSTauSTauEvent();
      tree.Fill();
      nevent++;
    }
  }
  tree.Write();
  rootFile->Close();
}
