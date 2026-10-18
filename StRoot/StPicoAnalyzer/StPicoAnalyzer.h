#ifndef StPicoAnalyzer_h
#define StPicoAnalyzer_h

// -----------------------------------------------------------------------------
// StPicoAnalyzer
// -----------------------------------------------------------------------------

#include <set>

#include "StMaker.h"

#include "StEmcUtil/geometry/StEmcGeom.h"

#include "TString.h"
#include "TProfile.h"

/*
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
*/
//#include "fastjet/contrib/SoftDrop.hh"

class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;
class TFile;
class TH1F;
class TH2F;

class StPicoAnalyzer : public StMaker {
 public:
  StPicoAnalyzer(const char* name, StPicoDstMaker* picoMaker, const char* outputFileName = "output.root");
  virtual ~StPicoAnalyzer();

  void SetBadTowerList(const char* f)    { mBadTowerFile = f; }
  void ApplyHadronicCorrection(StPicoDst* picoDst, StPicoEvent* picoEvent, std::vector<double>& towerE, double frac) const;

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t* opt = "");
  virtual Int_t Finish();

 protected:
  int RunId;

 private:
  StPicoDstMaker* mPicoDstMaker;  //! not owned, just a pointer we read from
  Long64_t mEventCounter;         //! events processed so far
  TString mOutputFileName;        //! e.g. "output.root" or "$JOBID.root"
  TFile* mOutputFile;             //! opened in Init(), written+closed in Finish()

  std::string   mBadTowerFile;   // comma-separated hot + dead tower IDs (1-4800)
  std::set<int> mBadTowers;

  StRefMultCorr  *refmultCorr; //currently only trying to use for bad run list access, assuming it's there

  bool   LoadBadTowers();
  bool   IsGoodTower(int id) const;

  std::map<TString, TH1*> mQAHists;

  TH1* BookTH1F(const char* name, const char* title, Int_t nBins, Double_t lo, Double_t hi);
  TH1* BookTH2F(const char* name, const char* title, Int_t nBinsX, Double_t xlo, Double_t xhi, Int_t nBinsY, Double_t ylo, Double_t yhi);
  TProfile* BookTProfile(const char* name, const char* title, Int_t nBinsX, Double_t xlo, Double_t xhi, Option_t* option = "");

  //EEC hists
  
  static const Int_t kNEECBins = 14;
  static const Double_t kEECBinEdges[kNEECBins + 1];  

  // Raw, unnormalized EEC numerator for several jet-pT intervals.
  // Normalize after the event loop using the corresponding integral/weight sum.
  TH1F* mEEC_low;
  TH1F* mEEC_mid;
  TH1F* mEEC_high;

  // Sum of E_i E_j / pT,jet^2 over all accepted track pairs in each jet-pT bin.
  TH1F* mEECWeight_low;
  TH1F* mEECWeight_mid;
  TH1F* mEECWeight_high;

  // Number of accepted jets in each pT interval.
  TH1F* mNJets;

  //cuts
  /*
  static const Double_t kJetR = 0.6;
  static const Double_t kTrackPtMin = 0.2;
  static const Double_t kTrackPtMax = 30.0;
  static const Double_t kTrackEtaMax = 1.0;
  static const Double_t kJetPtMin = 1.0; // for initial small statistics tests
  static const Double_t kJetPtMax = 50.0;
  static const Double_t kVertexZMax = 30.0;
  static const Int_t    kMinNHits = 20;
  static const Double_t kMinHitRatio = 0.52;
  static const Double_t kMaxDca = 1.0;
*/
  
  //analysis-specific functions
  int isTrackGood (const StPicoTrack* trk, const StPicoEvent* picoEvent) const;
  void processJets(StPicoDst* picoDst, StPicoEvent* picoEvent);
  //  void fillEEC(const class fastjet::PseudoJet& jet);//couldn't have this as member of class because of compile error with fastjet on RCF, so just defined separately from class in the cxx.
  
  TH1F* mTestHist;                //! test histogram -- replace with real histograms later
  
  ClassDef(StPicoAnalyzer, 2)
};

#endif
