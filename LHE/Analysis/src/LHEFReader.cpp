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
  if(argc < 3)
    {
      cout << " Usage: " << argv[0] << " output_file" << " input_files" << endl;
      cout << " output_file - output file in ROOT format." << endl;
      cout << " input_files - input files in LHEF format," << endl;
      return 1;
    }

  LHEFReader theReader;

  TFile* rootFile = TFile::Open( TString(argv[1]) , "recreate");
  TTree tree("LHEFiles" , "");
  tree.Branch("CharginoCharginoEvent" , &(theReader.CharginoCharginoEvent ) );


  for(int input_index = 2 ; input_index < argc ; input_index ++){

    theReader.SetReader( argv[input_index] );

    cout << "** Calculating number of events to process. Please wait..." << endl;
    Long64_t allEntries = theReader.reader->getNumberOfEvents();
    cout << "** Input file contains " << allEntries << " events" << endl;

    int nevent = 0;

    while( theReader.LoadNextEvent() ){
      //cout << nevent << endl;
      theReader.GetCharginoCharginoEvent();
      tree.Fill();
      nevent++;
    }
  }
  tree.Write();
  rootFile->Close();
}
