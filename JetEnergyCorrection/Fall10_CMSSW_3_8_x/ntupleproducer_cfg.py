import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("NTupleProducer")

### Message Logger #############################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append('NTP')
process.MessageLogger.cerr.NTP = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.categories.append('EcalSeverityLevelError')
process.MessageLogger.cerr.EcalSeverityLevelError = cms.untracked.PSet(
    limit = cms.untracked.int32(1),
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
                                      fileMode = cms.untracked.string("NOMERGE")
                                    )

### Parsing of command line parameters #############################################
### (type of run: data, MC; reconstruction: RECO, PAT, PF) #####################
options = VarParsing.VarParsing ('standard') # set 'standard'  options
options.register ('runon', # register 'runon' option
                  'data',  # the default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "Type of sample to run on: data (default), MC")
# get and parse the command line arguments
# set NTupleProducer defaults (override the output, files and maxEvents parameter)

#options.files= 'file:/scratch/pnef/mc/Fall11/AOD/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_007E0957-964E-E011-BA72-485B39800BBB.root'
#options.files= '/store/data/Run2011A/SingleElectron/AOD/PromptReco-v1/000/161/311/12418487-D557-E011-BCFA-001D09F24E39.root'
#options.files= '/store/data/Run2011A/MuOnia/AOD/PromptReco-v1/000/161/311/E6B512D1-CD57-E011-BC63-001D09F2423B.root'
#options.files= '/store/data/Run2010A/EG/AOD/Apr21ReReco-v1/0003/D6055361-7770-E011-BE49-00266CF32F90.root'
#options.files= 'file:/scratch/pnef/mc/Summer11-QCD_Pt-300to470_TuneZ2_7TeV_pythia6-AODSIM-PU_S3_START42_V11-v2-0000-4E79FBE4-F37C-E011-B178-003048D4607A.root'
#options.files= 'file:/shome/pnef/SUSY/reco-data/data//Run2011A/HT/AOD/PromptReco-v4/000/165/102/C49C75EC-CF80-E011-9BA4-003048F110BE.root'
options.files= 'file:/shome/pnef/SUSY/reco-data/mc/Summer11-QCD_Pt-300to470_TuneZ2_7TeV_pythia6-AODSIM-PU_S3_START42_V11-v2-0000-4E79FBE4-F37C-E011-B178-003048D4607A.root'
options.maxEvents = -1 # If it is different from -1, string "_numEventXX" will be added to the output file name

# Now parse arguments from command line (might overwrite defaults)
options.parseArguments()
options.output='NTupleProducer_42X_'+options.runon+'.root'

### Running conditions #########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
# to check what cff to use
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if options.runon=='data':
    # CMSSW_4_2
    process.GlobalTag.globaltag = "GR_R_42_V12::All"
    # CMSSW_3_8_X:
#    process.GlobalTag.globaltag = "GR_R_311_V2::All"
else:
    # CMSSW_3_8_X:
    # process.GlobalTag.globaltag = "START311_V2::All"
    process.GlobalTag.globaltag = "START41_V0::All"
    # CMSSW_4_2_X:
    # process.GlobalTag.globaltag = "START42_V12::All"


### Input/Output ###############################################################
# Input
process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(options.files)
      # Enable if you see duplicate error:
      # duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
# Output
process.TFileService = cms.Service("TFileService",
# Keep track of the type of data source and reco type in the ntuple file name
	fileName = cms.string(options.output),
	closeFileFast = cms.untracked.bool(True)
)

### Electron ID ##############################################################
process.load("DiLeptonAnalysis.NTupleProducer.simpleEleIdSequence_cff")

### Jet Corrections ##########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetEnergyCorrections
# note: this runs the L1Fast-Jet corrections for PF jets. not applied on Calo
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
# Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True
# process.kt6PFJets.Rho_EtaMax   = cms.double(4.4) # this is the default value in 4_2
# process.kt6PFJets.Ghost_EtaMax = cms.double(5.0) # this is the default value in 4_2
# Turn-on the FastJet jet area calculation 
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = process.kt6PFJets.Rho_EtaMax

### JES MET Corrections ########################################################
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet

process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
process.metMuonJESCorAK5.inputUncorJetsLabel = "ak5CaloJets"
process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"
process.metCorSequence = cms.Sequence(process.metMuonJESCorAK5)

### Cleaning ###################################################################
# flag HB/HE noise
# info: https://twiki.cern.ch/twiki/bin/view/CMS/HcalNoiseInfoLibrary
# uncomment following two lines to filer HCAL noise rather than to flag it!
# process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
# process.HBHENoiseFilter.maxRBXEMF = cms.double(0.01)
# the next two lines produce the HCAL noise summary flag with an additional cut on RBX
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.maxRBXEMF = cms.double(0.01)

# ECAL dead cells: this is not a filter. Only a flag is stored.
# Ecal gap boundary energy: specify minimal gap BE here.
### NOTE: THIS PART WAS REDUCED TO ENSURE COMPATIBILITY WITH CMSSW 3_9_7. THIS IS TEMPORARY!
#process.load("PhysicsTools/EcalAnomalousEventFilter/ecalanomalouseventfilter_cfi")
#process.EcalAnomalousEventFilter.FilterAlgo = cms.untracked.string("TuningMode")
#process.EcalAnomalousEventFilter.cutBoundEnergyGapEE = cms.untracked.double(5)
#process.EcalAnomalousEventFilter.cutBoundEnergyGapEB = cms.untracked.double(5)
#process.EcalAnomalousEventFilter.enableGap           = cms.untracked.bool(True)

# See for example DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py
# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
     applyfilter = cms.untracked.bool(True),
     debugOn = cms.untracked.bool(False),
     numtrack = cms.untracked.uint32(10),
     thresh = cms.untracked.double(0.25)
     )

process.goodVertices = cms.EDFilter("VertexSelector",
	src = cms.InputTag("offlinePrimaryVertices"),
	cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
	filter = cms.bool(False),
	)

############ PF2PAT ##########################################
# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

# Add a pro forma output module because PF2PAT complains otherwise...
process.out = cms.OutputModule("PoolOutputModule",
      # save only events passing the full path
      SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
      # save PAT Layer 1 output; you need a '*' to
      outputCommands = cms.untracked.vstring('drop *', *patEventContent )
      )

# Configure PAT to use PF2PAT instead of AOD sources
from PhysicsTools.PatAlgos.tools.pfTools import *

### Configuration in common to all collections
pfPostfixes = [ 'PF','PF2','PF3' ]
for pf in pfPostfixes:
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=(options.runon != 'data'), postfix=pf) 
    # turn to false when running on data and MC (for the moment)
    getattr(process, 'patElectrons'+pf).embedGenMatch = False
    getattr(process, 'patMuons'+pf).embedGenMatch = False

    # PatElectronID
    getattr(process, 'patElectrons'+pf).addElectronID = cms.bool(True)
    getattr(process, 'patElectrons'+pf).electronIDSources = cms.PSet(
	simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
	simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
	simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
      	simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
   	simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
       	simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    	simpleEleId95cIso  = cms.InputTag("simpleEleId95cIso"),
       	simpleEleId90cIso  = cms.InputTag("simpleEleId90cIso"),
    	simpleEleId85cIso  = cms.InputTag("simpleEleId85cIso"),
        simpleEleId80cIso  = cms.InputTag("simpleEleId80cIso"),
    	simpleEleId70cIso  = cms.InputTag("simpleEleId70cIso"),
        simpleEleId60cIso  = cms.InputTag("simpleEleId60cIso"),
    )
    # muon vertex 
    getattr(process, 'pfMuonsFromVertex'+pf).vertices=cms.InputTag("goodVertices") # require muon to come from the good vertices as defined above
    getattr(process, 'pfMuonsFromVertex'+pf).d0Cut   =cms.double(0.02) # transverse IP w.r.t. PV
    getattr(process, 'pfMuonsFromVertex'+pf).dzCut   =cms.double(1.)   # longitudinal IP w.r.t. PV
    getattr(process, 'pfSelectedMuons'+pf).cut = cms.string(
 		"abs( eta ) < 2.4 && pt > 5 && muonRef().isNonnull()" 
    )
    
    # electron vertex
    getattr(process, 'pfElectronsFromVertex'+pf).vertices=cms.InputTag("goodVertices") # require eles to come from the good vertices as defined above
    getattr(process, 'pfElectronsFromVertex'+pf).d0Cut   =cms.double(0.04) # transverse IP w.r.t. PV
    getattr(process, 'pfElectronsFromVertex'+pf).dzCut   =cms.double(1.)   # longitudinal IP w.r.t. PV
    getattr(process, 'pfSelectedElectrons'+pf).cut = cms.string(
 		"abs( eta ) < 2.4 && pt > 5 && gsfTrackRef().isNonnull()" 
		+"&& gsfTrackRef().trackerExpectedHitsInner().numberOfHits() <= 1"       # conv rejection
    )
    # ISOLATION
    getattr(process, 'pfIsolatedMuons'+pf)    .combinedIsolationCut = cms.double(0.20)
    getattr(process, 'pfIsolatedElectrons'+pf).combinedIsolationCut = cms.double(0.20)

    # Jet corrections 
    getattr(process,'patJetCorrFactors'+pf).levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
    getattr(process,'patJetCorrFactors'+pf).rho    = cms.InputTag('kt6PFJets','rho')
    getattr(process,'pfJets'+pf).doAreaFastjet     = True
    getattr(process,'pfJets'+pf).Rho_EtaMax        = process.kt6PFJets.Rho_EtaMax

    # Disable pileup removal on PF
    getattr(process,'pfNoPileUp'+pf).enable   = False

    # set to false to disable jet to be cleaned from Taus
    getattr(process,"pfNoTau"+pf).enable      = False

### Specific to first PF collection: MEDIUM
# ID Cuts
process.pfMuonsIDPF = cms.EDFilter("MuonIDPFCandidateSelector", 
		src = cms.InputTag("pfMuonsFromVertexPF"), 
		cut = cms.string( 
			"muonID('GlobalMuonPromptTight') && " + 
			"isTrackerMuon &&" + 
			"track.numberOfValidHits > 10 && "+ 
			"track.hitPattern.numberOfValidPixelHits > 0" 
			), 
		) 
process.pfSelectedMuonsPF.src = 'pfMuonsIDPF' 
process.pfMuonSequencePF.replace( process.pfSelectedMuonsPF, 
		process.pfMuonsIDPF * 
		process.pfSelectedMuonsPF 
		) 
process.pfElectronsIDPF = cms.EDFilter("ElectronIDPFCandidateSelector", 
		src = cms.InputTag("pfElectronsFromVertexPF"), 
		recoGsfElectrons = cms.InputTag('gsfElectrons'), 
		electronIdMap = cms.InputTag('simpleEleId95relIso'), 
		bitsToCheck = cms.vstring('id'), 
		) 
process.pfSelectedElectronsPF.src = 'pfElectronsIDPF' 
process.pfElectronSequencePF.replace( process.pfSelectedElectronsPF, 
		process.pfElectronsIDPF * 
		process.pfSelectedElectronsPF 
		) 

### Specific to second PF collection: TIGHT
# ID cuts
process.pfMuonsIDPF2 = cms.EDFilter("MuonIDPFCandidateSelector", 
		 src = cms.InputTag("pfMuonsFromVertexPF2"), 
		 cut = cms.string( 
			"muonID('GlobalMuonPromptTight') && " + 
			"isTrackerMuon && numberOfMatches > 1 && " + 
			"track.numberOfValidHits > 10 && "+ 
			"track.hitPattern.numberOfValidPixelHits > 0 &&" + 
			"globalTrack.ptError/pt < 0.1"
			), 
		 ) 
process.pfSelectedMuonsPF2.src = 'pfMuonsIDPF2' 
process.pfMuonSequencePF2.replace( process.pfSelectedMuonsPF2, 
		process.pfMuonsIDPF2 * 
		process.pfSelectedMuonsPF2 
		) 
process.pfElectronsIDPF2 = cms.EDFilter("ElectronIDPFCandidateSelector", 
		src = cms.InputTag("pfElectronsFromVertexPF2"), 
		recoGsfElectrons = cms.InputTag('gsfElectrons'), 
		electronIdMap = cms.InputTag('simpleEleId90relIso'), 
		bitsToCheck = cms.vstring('id'), 
		) 
process.pfSelectedElectronsPF2.src = 'pfElectronsIDPF2' 
process.pfElectronSequencePF2.replace( process.pfSelectedElectronsPF2, 
		process.pfElectronsIDPF2 * 
		process.pfSelectedElectronsPF2 
		) 

### Specific to second PF collection: LOOSE
# ID cuts
process.pfMuonsIDPF3 = cms.EDFilter("MuonIDPFCandidateSelector", 
		 src = cms.InputTag("pfMuonsFromVertexPF3"), 
		 cut = cms.string( 
			"muonID('GlobalMuonPromptTight') && " + 
			"isTrackerMuon && " + 
			"track.numberOfValidHits > 10 && "+ 
			"track.hitPattern.numberOfValidPixelHits > 0 " 
			), 
		 ) 
process.pfSelectedMuonsPF3.src = 'pfMuonsIDPF3' 
process.pfMuonSequencePF3.replace( process.pfSelectedMuonsPF3, 
		process.pfMuonsIDPF3 * 
		process.pfSelectedMuonsPF3 
		) 

### GenJets ####################################################################
# produce ak5GenJets (collection missing in case of some Spring10 samples)
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.mygenjets = cms.Sequence( process.genParticlesForJets * process.ak5GenJets )

### Analysis configuration #####################################################
process.load("DiLeptonAnalysis.NTupleProducer.ntupleproducer_cfi")
process.analyze.isRealData = cms.untracked.bool(options.runon=='data')

# Add some jet collections
process.analyze.jets = (
   # Calo jets
     cms.PSet( prefix = cms.untracked.string('CA'),
               tag = cms.untracked.InputTag('ak5CaloJets'),
               isPat = cms.untracked.bool(False),
               tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
               jet_id = cms.untracked.InputTag('ak5JetID'),
               sel_minpt  = process.analyze.sel_mincorjpt,
               sel_maxeta = process.analyze.sel_maxjeta,
               corrections = cms.string('ak5CaloL2L3'),
               btag_matchdeltaR = cms.double(0.25),
               ),
    # PF jets from PF2PAT
    cms.PSet( prefix = cms.untracked.string('PF2PAT'),
              tag = cms.untracked.InputTag('patJetsPF'),
              isPat = cms.untracked.bool(True),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              sel_minpt  = cms.double(15.0),
              sel_maxeta = process.analyze.sel_maxjeta,
              # The corrections are irrelevant for PF2PAT
              corrections = cms.string(''), 
              btag_matchdeltaR = cms.double(0.25),
              ),
    cms.PSet( prefix = cms.untracked.string('PF2PAT2'),
              tag = cms.untracked.InputTag('patJetsPF2'),
              isPat = cms.untracked.bool(True),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              sel_minpt  = cms.double(15.0),
              sel_maxeta = process.analyze.sel_maxjeta,
              # The corrections are irrelevant for PF2PAT
              corrections = cms.string(''), 
              btag_matchdeltaR = cms.double(0.25),
              ),
    cms.PSet( prefix = cms.untracked.string('PF2PAT3'),
              tag = cms.untracked.InputTag('patJetsPF3'),
              isPat = cms.untracked.bool(True),
              tag_jetTracks  = cms.untracked.InputTag('ak5JetTracksAssociatorAtVertex'),
              sel_minpt  = cms.double(15.0),
              sel_maxeta = process.analyze.sel_maxjeta,
              # The corrections are irrelevant for PF2PAT
              corrections = cms.string(''), 
              btag_matchdeltaR = cms.double(0.25),
              ),
    )
# # Add residual correction for running on data
# FIXME: NOT EVALUTED YET => NOT PRESENT IN 4_2 GLOBAL TAG FOR NOW
# if options.runon == 'data':
#         process.analyze.jetCorrs = process.analyze.jetCorrs.value() + 'Residual'
#         for extJet in process.analyze.jets:
#             extJet.corrections = extJet.corrections.value() + 'Residual'
#         for pf in pfPostfixes:            
#             getattr(process,'patJetCorrFactors'+pf).levels.extend( ['L2L3Residual'] )

# Add some PF lepton collections
process.analyze.leptons = (
    # PF Electrons
    cms.PSet( type = cms.untracked.string('electron'),
              prefix = cms.untracked.string('PfEl'),
              tag = cms.untracked.InputTag('patElectronsPF'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Muons
    cms.PSet( type = cms.untracked.string('muon'),
              prefix = cms.untracked.string('PfMu'),
              tag = cms.untracked.InputTag('patMuonsPF'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Taus
    cms.PSet( type = cms.untracked.string('tau'),
              prefix = cms.untracked.string('PfTau'),
              tag = cms.untracked.InputTag('patTausPF'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Electrons
    cms.PSet( type = cms.untracked.string('electron'),
              prefix = cms.untracked.string('PfEl2'),
              tag = cms.untracked.InputTag('patElectronsPF2'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Muons
    cms.PSet( type = cms.untracked.string('muon'),
              prefix = cms.untracked.string('PfMu2'),
              tag = cms.untracked.InputTag('patMuonsPF2'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Taus
    cms.PSet( type = cms.untracked.string('tau'),
              prefix = cms.untracked.string('PfTau2'),
              tag = cms.untracked.InputTag('patTausPF2'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Electrons
    cms.PSet( type = cms.untracked.string('electron'),
              prefix = cms.untracked.string('PfEl3'),
              tag = cms.untracked.InputTag('patElectronsPF3'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Muons
    cms.PSet( type = cms.untracked.string('muon'),
              prefix = cms.untracked.string('PfMu3'),
              tag = cms.untracked.InputTag('patMuonsPF3'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    # PF Taus
    cms.PSet( type = cms.untracked.string('tau'),
              prefix = cms.untracked.string('PfTau3'),
              tag = cms.untracked.InputTag('patTausPF3'),
              sel_minpt = process.analyze.sel_minelpt,
              sel_maxeta = process.analyze.sel_maxeleta,
              maxnobjs = cms.untracked.uint32(20)
              ),
    )
              
process.readAK5PF = cms.EDAnalyzer('JetCorrectorDBReader', 
	payloadName    = cms.untracked.string('AK5Calo'),
	printScreen    = cms.untracked.bool(False),
	createTextFile = cms.untracked.bool(True),
	globalTag      = cms.untracked.string('Fall10')
)

#### DEBUG #####################################################################
# process.dump = cms.EDAnalyzer("EventContentAnalyzer")
# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1) # number of events to ignore at start (default is one)
# )
# process.ProfilerService = cms.Service("ProfilerService",
#                                      firstEvent = cms.untracked.int32(2),
#                                      lastEvent = cms.untracked.int32(51),
#                                      paths = cms.untracked.vstring(['p'])
#                                      )
# process.Tracer = cms.Service("Tracer")

# to enable pileUpsubtraction for MET
# process.pfMETPF.src=cms.InputTag("pfNoPileUpPF")
# process.pfMETPF2.src=cms.InputTag("pfNoPileUpPF2")
#### Path ######################################################################
process.p = cms.Path(
    process.scrapingVeto * (
	process.goodVertices
        + process.HBHENoiseFilterResultProducer
        # + process.EcalAnomalousEventFilter
	+ process.kt6PFJets
	+ process.ak5PFJets
        + process.mygenjets
        + process.simpleEleIdSequence
        + process.metCorSequence
        + process.patPF2PATSequencePF
        + process.patPF2PATSequencePF2
        + process.patPF2PATSequencePF3
#	+ process.dump
        + process.analyze
        )
    )

process.p += process.readAK5PF


# remove ak5GenJets from the path if it will run on data
if options.runon=='data':
    process.p.remove(process.mygenjets)
if options.runon!='data':
    process.p.remove(process.scrapingVeto)

