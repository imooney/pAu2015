#include "StPicoAnalyzer.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

#include "StEmcCollection.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"

#include "StRefMultCorr/StRefMultCorr.h"

#include "StMessMgr.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TVector2.h"

#include "StThreeVectorF.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

//#include "fastjet/contrib/SoftDrop.hh"

#include <cmath>
#include <vector>

ClassImp(StPicoAnalyzer)

namespace {

  //cuts
  static const Double_t kJetR = 0.6;
  static const Double_t kTrackPtMin = 0.2;
  static const Double_t kTrackPtMax = 30.0;
  static const Double_t kTowerEtMin = kTrackPtMin;
  static const Double_t kTowerEtMax = kTrackPtMax;
  static const Double_t kTrackEtaMax = 1.0;
  static const Double_t kJetPtMin = 1.0; // for initial small statistics tests
  static const Double_t kJetPtMax = 50.0;
  static const Double_t kVertexZMax = 30.0;
  static const Double_t kVertexDiff = 3.0;
  static const Int_t    kMinNHits = 20;
  static const Double_t kMinHitRatio = 0.52;
  static const Double_t kMaxDca = 1.0;
  static const Double_t kEAmax = 64000.0;
  static const Double_t kNEFmax = 0.9;


  //basic quantities
  //const double kBemcRadius = 225.405;   // cm, BEMC tower face radius
  const int    kNTowers    = 4800;
  //const double kTwoPi      = 2.0 * M_PI;
  const Double_t mpi = 0.13957039; //mass of charged pion

  //tower info
  struct TowerHit {
    int    id;        // soft ID, 1-4800
    double eta, phi;  // relative to the primary vertex
    double Eraw;      // energy pre-hadronic-correction
  };

}

const Double_t StPicoAnalyzer::kEECBinEdges[StPicoAnalyzer::kNEECBins + 1] = {
    0.02,   0.0364, 0.0543, 0.0809, 0.0988, 0.1207, 0.1474, 0.1799,
    0.2197, 0.2683, 0.3276, 0.4,    0.4884, 0.5964, 0.7282
};

using namespace fastjet;

//these helper functions at the top rely on fastjet objects, so they're outside of the class for reasons
static void fillEEC(const fastjet::PseudoJet& jet, TH1F* hEEC, TH1F* hWeight);//defined later.

static Double_t neutralEnergyFraction(const fastjet::PseudoJet& jet)
{
  const std::vector<fastjet::PseudoJet> cons = jet.constituents();
  Double_t ptTower = 0., ptAll = 0.;
  for (size_t i = 0; i < cons.size(); ++i) {
    ptAll += cons[i].pt();
    if (cons[i].user_index() < 0) ptTower += cons[i].pt();   // towers: user_index = -softId
  }
  return (ptAll > 0.) ? ptTower / ptAll : -1.;
}

//_____________________________________________________________________________
StPicoAnalyzer::StPicoAnalyzer(const char* name, StPicoDstMaker* picoMaker, const char* outputFileName)
: StMaker(name), mPicoDstMaker(picoMaker), mEventCounter(0), mOutputFileName(outputFileName), mOutputFile(0), mTestHist(0), mEEC_low(0), mEEC_mid(0), mEEC_high(0), mEECWeight_low(0), mEECWeight_mid(0), mEECWeight_high(0), mNJets(0)
{
}

//_____________________________________________________________________________
StPicoAnalyzer::~StPicoAnalyzer()
{
}

bool StPicoAnalyzer::LoadBadTowers()
{
  mBadTowers.clear();
  std::ifstream in(mBadTowerFile.c_str());
  if (!in.is_open()) {
    Error("LoadBadTowers", "cannot open bad tower list '%s'", mBadTowerFile.c_str());
    return false;
  }
  std::string line;
  while (std::getline(in, line)) {
    const size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    std::replace(line.begin(), line.end(), ',', ' ');   // comma-separated -> whitespace
    std::istringstream ss(line);
    int id;
    while (ss >> id) if (id >= 1 && id <= kNTowers) mBadTowers.insert(id);
  }
  Info("LoadBadTowers", "loaded %d bad towers from %s",
       (int)mBadTowers.size(), mBadTowerFile.c_str());
  return true;
}

bool StPicoAnalyzer::IsGoodTower(int id) const
{
  //wants id to be softid (1 to 4800)
  return id >= 1 && id <= kNTowers && mBadTowers.find(id) == mBadTowers.end();
}

// towerE: raw energy per BTowHit index on input (0 for bad/empty towers), corrected on output
void StPicoAnalyzer::ApplyHadronicCorrection(StPicoDst* picoDst, StPicoEvent* picoEvent,
                                             std::vector<double>& towerE, double frac) const {

  for (size_t k = 0; k < picoDst->numberOfTracks(); ++k) {
    const StPicoTrack* trk = picoDst->track(k);
    if (!isTrackGood(trk,picoEvent)) continue;                          // looping over all tracks, but specifying only the good ones
    const int idx = trk->bemcTowerIndex();                    // mBEmcMatchedTowerIndex = 0 if no matched tower, < 0 if not exact match, and bemcTowerIndex() returns abs(mBEmcMatchedTowerIndex) - 1, so -1 means no match, 0 to 4799 means some kind of match.
    if (idx < 0 || idx >= (int)towerE.size()) continue;       // if idx >= 4800 something's wrong
    if (towerE[idx] <= 0) continue;                           // bad, empty, or already zeroed
    towerE[idx] -= frac * trk->pMom().Mag();
    if (towerE[idx] < 0) towerE[idx] = 0;
  }//end loop over tracks

}




//_____________________________________________________________________________
Int_t StPicoAnalyzer::Init()
{
  //IMPORTANT FOR UNCERTAINTIES FOR EEC: ensures correct error: w*sqrt(N) rather than default: sqrt(w*N).
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  LOG_INFO << "StPicoAnalyzer::Init() - starting up, output -> " << mOutputFileName.Data() << endm;
  mEventCounter = 0;
  
  mOutputFile = new TFile(mOutputFileName.Data(), "RECREATE");
  if (!mOutputFile || mOutputFile->IsZombie()) {
    LOG_ERROR << "StPicoAnalyzer::Init() - could not open output file "
	      << mOutputFileName.Data() << endm;
    return kStErr;
  }

  //mEmcGeom = StEmcGeom::instance("bemc");
  if (!LoadBadTowers()) return kStFatal;

  refmultCorr = new StRefMultCorr("refmult");

  mTestHist = new TH1F("h","test histogram",10,0,10);

  //~~~~~~~~~~~~~~~~~~~~~~QA~~~~~~~~~~~~~~~~~~~~//
  //event-level
  BookTH1F("hVz", "Primary vertex z;v_{z} (cm);Events",
	   120, -60., 60.);
  BookTH1F("hNTracks", "Number of tracks;N tracks;Events",
	   300, -0.5, 2999.5);
  BookTH1F("hNPrimTracks", "Number of primary tracks;N tracks;Events",
	   300, -0.5, 299.5);
  BookTH1F("hNGoodTracks", "Number of accepted tracks;N tracks;Events",
	   300, -0.5, 299.5);
  BookTH1F("hNTowersFired","Number of towers fired;N towers;Events",
	   500,-0.5,4999.5);
  BookTH1F("hBBCchannels", "Activity per channel, BBC;i;Total activity/channel",
	   51,-0.5, 50.5);
  BookTH1F("hiBBCEsum", "#sum ADC in inner east BBC;#sum EA^{BBCeast};Events",
	   1000,0., 100000.);
  BookTH1F("hVzDiff","| v_{z,VPD} - v_{z,TPC} |;#Delta v_{z}; Events",
	   100,-0.5,99.5);
  BookTH2F("hZDCx","ZDC coincidence rate vs. run ID;ZDC&; runID; Events",
	   100,0,100000,100,16000000,17000000); //fix this weird plot later
  //track-level
  //quality
  BookTH1F("hTrackNHitsFit", "N_{hits};N_{hits};Tracks",
	   100,0.,100.);
  BookTH1F("hTrackNHitsRatio", "N_{hits} ratio; N_{hits}/N_{hits-max};Tracks",
	   100,0.,2.);
  BookTH1F("hTrackDCA", "DCA; DCA [cm];Tracks",
	   100,0.,5.);
  //kinematics
  BookTH1F("hTrackPt","Track p_{T};p_{T} [GeV/c];Tracks",
	   100, 0, 32);
  BookTH1F("hTrackEta","Track #eta;#eta;Tracks",
	   100, -2, 2);
  BookTH1F("hTrackPhi","Track #phi;#phi;Tracks",
	   100, -7, 7);
  BookTH2F("hTrackEtaPhi","Track #eta vs. #phi;#eta;#phi;Tracks",
	   100, -2, 2, 100, -7, 7);
  BookTH1F("hTrackCharge","Track charge;Q;Tracks",
	   5, -2.5, 2.5);
  
  //tower-level
  BookTH1F("hTowerId","Tower hits;Tower ID;Total hits",
	   kNTowers+2, -0.5, kNTowers+1.5);
  BookTH1F("hTowerET","Tower ET;E_{T} [GeV];Total E_{T}",
	   100, 0, 32);
  BookTH1F("hTowerETcorr","Tower ET post-HC;E_{T}^{corr} [GeV];Total E_{T}^{corr}",
	   100, 0, 32);

  BookTH2F("hETrawVsCorr", ";E_{T}^{raw} [GeV];E_{T}^{corr} [GeV]",
			  100, 0, 20, 100, 0, 20);
  BookTH1F("hFracRemoved", ";(E_{T}^{raw}-E_{T}^{corr})/E_{T}^{raw};towers",
			  110, 0, 1.1);
  
  BookTProfile("pTowerETraw",  ";tower ID;<E_{T}^{raw}> [GeV]",  4801, -0.5, 4800.5);
  BookTProfile("pTowerETcorr", ";tower ID;<E_{T}^{corr}> [GeV]", 4801, -0.5, 4800.5);
  BookTProfile("pTowerFrac",   ";tower ID;<E_{T}^{corr}/E_{T}^{raw}>", 4801, -0.5, 4800.5);

  BookTH2F("hTowerEtaPhi","Tower #eta vs. #phi;#eta;#phi;Towers",
	   100, -2, 2, 100, -7, 7);
  
  //jet-level
  BookTH1F("hJetPt", "anti-k_{T} jets;p_{T,jet} (GeV/c);Jets",
	   100, 0., 60.);
  BookTH1F("hJetEta", "Accepted jets;#eta_{jet};Jets",
	   100, -1., 1.);
  BookTH1F("hJetNcons", "Jet constituent multiplicity;N_{cons};Jets", 
	   50, 0., 50.);
  BookTH1F("hJetNEF", "Jet neutral energy fraction;NEF;Jets", 55, 0., 1.1);

  //EEC
  /*
  const Double_t rMax = 2*kJetR;
  const Int_t nBins = rMax*100;
  */
  
  mEEC_low = new TH1F("hEEC_low",
                       "Raw EEC numerator, 15<p_{T,jet}<20;R_{L};#Sigma E_{i}E_{j}/p_{T,jet}^{2}",
		      kNEECBins, kEECBinEdges);//nBins, 0., rMax);
  mEEC_mid = new TH1F("hEEC_mid",
                       "Raw EEC numerator, 20<p_{T,jet}<30;R_{L};#Sigma E_{i}E_{j}/p_{T,jet}^{2}",
		      kNEECBins, kEECBinEdges);//nBins, 0., rMax);
  mEEC_high = new TH1F("hEEC_high",
                       "Raw EEC numerator, 30<p_{T,jet}<50;R_{L};#Sigma E_{i}E_{j}/p_{T,jet}^{2}",
                       kNEECBins, kEECBinEdges);//nBins, 0., rMax);

  mEECWeight_low = new TH1F("hEECWeight_low",
                             "Total EEC weight, 15<p_{T,jet}<20;dummy;weight",
			    1, 0., 1.);
  mEECWeight_mid = new TH1F("hEECWeight_mid",
                             "Total EEC weight, 20<p_{T,jet}<30;dummy;weight",
                             1, 0., 1.);
  mEECWeight_high = new TH1F("hEECWeight_high",
                             "Total EEC weight, 30<p_{T,jet}<50;dummy;weight",
                             1, 0., 1.);

  mNJets = new TH1F("hNJets",
                    "Accepted jets by p_{T} interval;bin;jets",
                    3, 0., 3.);
  mNJets->GetXaxis()->SetBinLabel(1, "low");
  mNJets->GetXaxis()->SetBinLabel(2, "mid");
  mNJets->GetXaxis()->SetBinLabel(3, "high");
  
  return kStOK;
}

//_____________________________________________________________________________
Int_t StPicoAnalyzer::Make()
{
  if (!mPicoDstMaker) {
    LOG_ERROR << "StPicoAnalyzer::Make() - no StPicoDstMaker!" << endm;
    return kStErr;
  }

  StPicoDst* picoDst = mPicoDstMaker->picoDst();
  if (!picoDst) {
    LOG_ERROR << "StPicoAnalyzer::Make() - no PicoDst!" << endm;
    return kStErr;
  }

  StPicoEvent* picoEvent = picoDst->event();
  if (!picoEvent) {
    LOG_ERROR << "StPicoAnalyzer::Make() - no PicoEvent!" << endm;
    return kStErr;
  }

  //bad runs!
  RunId=picoEvent->runId();
  //refmultcorr doesn't actually work for pAu, so keeping this commented out to prevent annoying error messages although it shouldn't hurt anything.
  /*
  refmultCorr->init(RunId);
  if(refmultCorr->isBadRun(RunId)) {
    LOG_INFO << "RUN " << RunId << " is officially bad" << endm;
    return kStOK; //shouldn't see this for 2015 pAu since bad runs aren't defined in the StRefMultCorr BadRun.h file
  }
  */
  
  //LOG_INFO << RunId << " not in 16142059 -  16149001, and neq 16135031 || 16135032, right?" << endm;
  //doing it by hand for now -- can do it more thoroughly later
  if (RunId >= 16142059 && RunId <= 16149001) { LOG_INFO << "RUN " << RunId << " failed by-hand bad-run check" << endm; return kStOK;}
  //something weird happened to the towers in run 16135032 (and it looks like it started at the end of run 16135031), so excluding both 
  if (RunId == 16135031 || RunId == 16135032) { LOG_INFO << "RUN " << RunId << " failed by-hand bad-run check" << endm; return kStOK;}

  //trigger!
  if ((picoEvent->isTrigger(500401))==0 && (picoEvent->isTrigger(500411))==0) {//JP2 for pAu2015 200 GeV
    return kStOK;
  }
  
  ++mEventCounter;

  //Eventually this could become its own function representing all event-level selections
  Double_t iBBCSum = 0.0;
  for (int i = 0; i < 30; ++i) {//looping over all (24) tiles (+a few to check the mapping)
    if (i < 16) {iBBCSum += picoEvent->bbcAdcEast(i);}//specifically inner tiles
    static_cast<TH1F*>(mQAHists["hBBCchannels"])->Fill(i,picoEvent->bbcAdcEast(i));
  }
  static_cast<TH1F*>(mQAHists["hiBBCEsum"])->Fill(iBBCSum);//pre-cut for now, to check
  
  if (iBBCSum > kEAmax)
    return kStOk;
  
  const TVector3 vtx = picoEvent->primaryVertex();
  const Float_t vpdvz = picoEvent->vzVpd();
  
  //event-level QA
  static_cast<TH1F*>(mQAHists["hVz"])->Fill(vtx.Z());
  static_cast<TH1F*>(mQAHists["hVzDiff"])->Fill(vtx.Z() - vpdvz);
  static_cast<TH1F*>(mQAHists["hNTracks"])->Fill(picoDst->numberOfTracks());
  static_cast<TH1F*>(mQAHists["hNPrimTracks"])->Fill(picoEvent->numberOfPrimaryTracks());
  static_cast<TH1F*>(mQAHists["hZDCx"])->Fill(picoEvent->ZDCx(), RunId);
  
  if (std::fabs(vtx.Z()) > kVertexZMax)
    return kStOK;

  if (std::fabs(vtx.Z() - vpdvz) > kVertexDiff)
    return kStOk;
  
  processJets(picoDst, picoEvent);
  
  // This is the "do nothing but print" placeholder. Replace this block with
  // real analysis: loop over picoDst->numberOfTracks(), apply cuts, fill
  // histograms/trees, etc.
  /*LOG_INFO << "StPicoAnalyzer::Make() - event# " << mEventCounter
            << " | run " << picoEvent->runId()
            << " | eventID " << picoEvent->eventId()
            << " | refMult " << picoEvent->refMult()
            << " | nTracks " << picoDst->numberOfTracks()
            << endm;
  
  mTestHist->Fill(5);
  */
  return kStOK;
}

//_____________________________________________________________________________
void StPicoAnalyzer::Clear(Option_t* opt)
{
  StMaker::Clear(opt);
}

//_____________________________________________________________________________
Int_t StPicoAnalyzer::Finish()
{
  LOG_INFO << "StPicoAnalyzer::Finish() - processed " << mEventCounter
            << " events total" << endm;

  if (mOutputFile) {
    mOutputFile->cd();
    TDirectory* qaDir = mOutputFile->mkdir("QA");
    qaDir->cd();
    for (const auto& kv : mQAHists) {
      kv.second->Write();
    }
    mOutputFile->cd();
    mOutputFile->Write();
    mOutputFile->Close();
    delete mOutputFile;
    mOutputFile = 0;
  }
  return kStOK;
}

//~~~~~~~~~~~~~~~~~~~~HELPERS~~~~~~~~~~~~~~~~~~~~~~//

int StPicoAnalyzer::isTrackGood (const StPicoTrack* trk, const StPicoEvent* picoEvent) const {//returns -1 if event needs to be skipped
  if (!trk) return 0;
  TVector3 p;
  if(!(trk->isPrimary())) return 0;
  p = trk->pMom();  
  /*
  if(doUseTofMatch){
    if(!(trk->isTofTrack()))return kFALSE;
  }
  */

  double pt = p.Perp();
  //double phi = p.Phi(); //not needed as no selections applied on it given the 2pi acceptance in phi
  
  //track selections
  //    const Double_t pt = trk->gPt();
  if (pt < kTrackPtMin) return 0;
  
  if (pt > kTrackPtMax) return -1;//this will effectively skip this event as nothing is filled until after this loop, which is aborted. Will do similarly for towers once implemented

    //const StThreeVectorF& p = trk->gMom();
    const Double_t pMag = p.Mag();//p.mag();
    if (pMag <= 0) return 0;

    const Double_t eta = p.PseudoRapidity();//0.5 * std::log((pMag + p.z()) / (pMag - p.z() + 1.e-12));
    if (std::fabs(eta) >= kTrackEtaMax) return 0;

    if (trk->nHitsFit() < kMinNHits) return 0;

    if (trk->nHitsMax() <= 0) return 0;

    const Double_t nHitsRatio = static_cast<Double_t>(trk->nHitsFit()) /
      static_cast<Double_t>(trk->nHitsMax());
    if (nHitsRatio < kMinHitRatio) return 0;

    double dca = (trk->gDCA (picoEvent->primaryVertex())).Mag();

    if (dca >= kMaxDca) return 0;

    return 1;
    
}

//_____________________________________________________________________________
void StPicoAnalyzer::processJets(StPicoDst* picoDst, StPicoEvent* picoEvent)
{
  std::vector<fastjet::PseudoJet> particles;
  particles.reserve(picoDst->numberOfTracks());
  int goodtrk_counter = 0;
  for (UInt_t i = 0; i < picoDst->numberOfTracks(); ++i) {
    StPicoTrack* trk = picoDst->track(i);
    if (isTrackGood(trk, picoEvent) < 0) return; //skipping this event because track pT > 30 GeV.
    if (isTrackGood(trk, picoEvent) == 0) continue; //if it returns 0, track is bad, so we skip to the next one

    //passed all selections
    goodtrk_counter ++;
        
    //track QA (post-selections)
    //quality
    static_cast<TH1F*>(mQAHists["hTrackNHitsFit"])->Fill(trk->nHitsFit());
    static_cast<TH1F*>(mQAHists["hTrackNHitsRatio"])->Fill(static_cast<Double_t>(trk->nHitsFit()) /
							   static_cast<Double_t>(trk->nHitsMax()));
    static_cast<TH1F*>(mQAHists["hTrackDCA"])->Fill((trk->gDCA (picoEvent->primaryVertex())).Mag());
    
    const Double_t energy = std::sqrt((trk->pMom()).Mag()*(trk->pMom()).Mag() + mpi*mpi);

    //kinematics
    static_cast<TH1F*>(mQAHists["hTrackPt"])->Fill((trk->pMom()).Perp());
    static_cast<TH1F*>(mQAHists["hTrackEta"])->Fill((trk->pMom()).PseudoRapidity());
    static_cast<TH1F*>(mQAHists["hTrackPhi"])->Fill((trk->pMom()).Phi());
    static_cast<TH2F*>(mQAHists["hTrackEtaPhi"])->Fill((trk->pMom()).PseudoRapidity(),(trk->pMom()).Phi());
    static_cast<TH1F*>(mQAHists["hTrackCharge"])->Fill(trk->charge());

    fastjet::PseudoJet pj((trk->pMom()).x(), (trk->pMom()).y(), (trk->pMom()).z(), energy);
    pj.set_user_index(static_cast<int>(i));
    particles.push_back(pj);
  }//end for loop over particles, conversion to pseudojets
  static_cast<TH1F*>(mQAHists["hNGoodTracks"])->Fill(goodtrk_counter);
  
  if (particles.empty()) {
    LOG_INFO << "No good tracks in this event" << endm;
    return;
  }
  
  StEmcGeom* bemcGeom = StEmcGeom::getEmcGeom("bemc");

  //for now, putting a tower loop here -- later, will need to associate towers with tracks, and add hadronic-corrected towers to the jets. Currently just doing some QA here
  Int_t nTowersFired = 0;
  std::vector<double> towerE(kNTowers, 0.0);
  std::vector<TowerHit> hits;                        // only surviving towers
  hits.reserve(picoDst->numberOfBTowHits());
  for (UInt_t itow = 0; itow < picoDst->numberOfBTowHits(); ++itow) {
    if (itow >= (UInt_t)kNTowers) { LOG_ERROR << "nBTowHits = " << picoDst->numberOfBTowHits() << endm; break; }
    StPicoBTowHit* tow = picoDst->btowHit(itow);
    if (!tow || tow->energy() <= kTowerEtMin ) continue;//avoids computation to check if E < max here first because then ET < max without having to calculate ET.
    //need: towerID, bemcGeom, ...
    float xT, yT, zT;
    bemcGeom->getXYZ(tow->numericIndex2SoftId(itow), xT, yT, zT);
    TVector3 towLoc {xT, yT, zT};
    TVector3 relPos { towLoc - picoEvent->primaryVertex()   };
    float towerEta  { static_cast<float>(relPos.PseudoRapidity())};
    if (fabs(towerEta) > 1.0) continue;
    float towerPhi  { static_cast<float>(TVector2::Phi_mpi_pi(relPos.Phi()))};//__0to2pi(relPos.Phi())) };
    float towerEt = (float) tow->energy() / (float) TMath::CosH(towerEta);
    if (towerEt < kTowerEtMin) continue;
    //if (towerEt >= kTowerEtMax) return; //quits the event, to be consistent with tracking. Doing this after the hadronic correction now, so track energy isn't double-counted with events thrown away unnecessarily.
    if (!IsGoodTower(tow->numericIndex2SoftId(itow))) continue; //if it's a dead/hot/bad status tower, skip it.
    towerE[itow] = tow->energy();//check all indexing later. Prevent off-by-ones. //think I was getting off-by-one before with "tow->numericIndex2SoftId(itow)"
    nTowersFired++;
    static_cast<TH1F*>(mQAHists["hTowerId"])->Fill(tow->numericIndex2SoftId(itow));
    static_cast<TH1F*>(mQAHists["hTowerET"])->Fill(towerEt);
    static_cast<TH2F*>(mQAHists["hTowerEtaPhi"])->Fill(towerEta,towerPhi);
    //static_cast<TH2F*>(mQAHists["hTowerIdEnergy"])->Fill(tow->id(), tow->energy());

    TowerHit th = { tow->numericIndex2SoftId(itow), towerEta, towerPhi, tow->energy() };
    hits.push_back(th);
  }
  static_cast<TH1F*>(mQAHists["hNTowersFired"])->Fill(nTowersFired);
  
  StPicoAnalyzer::ApplyHadronicCorrection(picoDst, picoEvent, towerE, 1.0);//100% hadronic correction
  
  // ---- second pass over the surviving towers only ----
  for (size_t i = 0; i < hits.size(); ++i) {
    const double Ecorr = towerE[hits[i].id - 1];//indexed to 0 instead of 1
    
    const double coshEta = std::cosh(hits[i].eta);
    const double ETraw   = hits[i].Eraw / coshEta;
    const double ETcorr  = Ecorr / coshEta;                // 0 for fully subtracted towers

    // raw vs corrected comparison, filled for every surviving tower
    static_cast<TH2F*>(mQAHists["hETrawVsCorr"])->Fill(ETraw, ETcorr);
    static_cast<TH1F*>(mQAHists["hFracRemoved"])->Fill((ETraw - ETcorr) /(double) ETraw);
    static_cast<TProfile*>(mQAHists["pTowerETraw"])->Fill(hits[i].id, ETraw); //softid here
    static_cast<TProfile*>(mQAHists["pTowerETcorr"])->Fill(hits[i].id, ETcorr);
    static_cast<TProfile*>(mQAHists["pTowerFrac"])->Fill(hits[i].id, ETcorr / (double) ETraw);

    if (ETcorr < kTowerEtMin) continue; //towers with ET < 0.2 GeV won't be clustered into jets.
    if (ETcorr > kTowerEtMax) return; //events with towers > 30 GeV 

    static_cast<TH1F*>(mQAHists["hTowerETcorr"])->Fill(ETcorr);
    fastjet::PseudoJet pj(ETcorr * std::cos(hits[i].phi), ETcorr * std::sin(hits[i].phi),
                          ETcorr * std::sinh(hits[i].eta), Ecorr);   // massless
    pj.set_user_index(-hits[i].id);
    particles.push_back(pj);
  }
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, kJetR);
  fastjet::ClusterSequence cs(particles, jetDef);

  std::vector<fastjet::PseudoJet> jets =
      fastjet::sorted_by_pt(cs.inclusive_jets(kJetPtMin));

  for (std::vector<fastjet::PseudoJet>::const_iterator it = jets.begin();
       it != jets.end(); ++it) {

    const fastjet::PseudoJet& jet = *it;

    if (jet.pt() <= kJetPtMin || jet.pt() >= kJetPtMax)
      continue;

    //fiducial acceptance cut
    if (std::fabs(jet.eta()) > kTrackEtaMax - kJetR)
      continue;

    const Double_t nef = neutralEnergyFraction(jet);
    static_cast<TH1F*>(mQAHists["hJetNEF"])->Fill(nef);
    if (nef < 0. || nef > kNEFmax) continue; //keeps only jets with 0 < NEF < NEFmax (0.9)

    static_cast<TH1F*>(mQAHists["hJetPt"])->Fill(jet.pt());
    static_cast<TH1F*>(mQAHists["hJetEta"])->Fill(jet.eta());
    static_cast<TH1F*>(mQAHists["hJetNcons"])->Fill(jet.constituents().size());

    const Double_t pt = jet.pt();
    if (pt > 15. && pt < 20.)
      mNJets->Fill(0.5);
    else if (pt >= 20. && pt < 30.)
      mNJets->Fill(1.5);
    else if (pt >= 30. && pt < 50.)
      mNJets->Fill(2.5);
    
    TH1F* hEEC = 0;
    TH1F* hWeight = 0;
    
    if (jet.pt() > 15. && jet.pt() < 20.) {
      hEEC = mEEC_low;
      hWeight = mEECWeight_low;
    } else if (jet.pt() >= 20. && jet.pt() < 30.) {
      hEEC = mEEC_mid;
      hWeight = mEECWeight_mid;
    } else if (jet.pt() >= 30. && jet.pt() < 50.) {
      hEEC = mEEC_high;
      hWeight = mEECWeight_high;
    } else {
      continue;
    }

    fillEEC(jet, hEEC, hWeight);
  }
}

//_____________________________________________________________________________
static void /*StPicoAnalyzer::*/fillEEC(const fastjet::PseudoJet& jet, TH1F* hEEC, TH1F* hWeight)
{
  const std::vector<fastjet::PseudoJet> constituents = jet.constituents();
  if (constituents.size() < 2)
    return;

  const Double_t jetPt2 = jet.pt() * jet.pt();

  //looping over distinct pairs rather than all pairs (i < j, rather than i != j). Either way is fine physics/math-wise.
  for (size_t i = 0; i < constituents.size(); ++i) {
    for (size_t j = i + 1; j < constituents.size(); ++j) {
      if (constituents[i].user_index() < 0 || constituents[j].user_index() < 0 ) { continue; }//only charged tracks in the correlator by constrruction
      const Double_t dEta = constituents[i].eta() - constituents[j].eta();
      const Double_t dPhi =
          TVector2::Phi_mpi_pi(constituents[i].phi() - constituents[j].phi());

      const Double_t RL = std::sqrt(dEta*dEta + dPhi*dPhi);
      const Double_t weight =
          constituents[i].E() * constituents[j].E() / jetPt2;

      if (weight <= 0.)
        continue;

      hEEC->Fill(RL, weight);
      hWeight->Fill(0.5, weight);
    }
  }

}

TH1* StPicoAnalyzer::BookTH1F(const char* name, const char* title,
                               Int_t nBins, Double_t lo, Double_t hi)
{
    TH1F* h = new TH1F(name, title, nBins, lo, hi);
    h->SetDirectory(0);   // detach from automatic directory ownership so later can go in its own TDirectory
    mQAHists[name] = h;
    return h;
}

TH1* StPicoAnalyzer::BookTH2F(const char* name, const char* title,
                               Int_t nBinsX, Double_t xlo, Double_t xhi,
                               Int_t nBinsY, Double_t ylo, Double_t yhi)
{
    TH2F* h = new TH2F(name, title, nBinsX, xlo, xhi, nBinsY, ylo, yhi);
    h->SetDirectory(0);   // detach from automatic directory ownership so later can go in its own TDirectory
    mQAHists[name] = h;
    return h;
}

TProfile* StPicoAnalyzer::BookTProfile(const char* name, const char* title,
                                       Int_t nBinsX, Double_t xlo, Double_t xhi,
                                       Option_t* option)
{
  TProfile* p = new TProfile(name, title, nBinsX, xlo, xhi, option);
  p->SetDirectory(0);   // detach from automatic directory ownership so it can go in its own TDirectory later
  mQAHists[name] = p;
  return p;
}
