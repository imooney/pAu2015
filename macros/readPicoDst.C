// -----------------------------------------------------------------------------
// readPicoDst.C
//
// Steering macro for a minimal PicoDst skeleton analysis.
//
// Usage (from within a starver'd ROOT/root4star session, after `cons` has
// built StPicoAnalyzer):
//
//   root4star -b -q 'readPicoDst.C("input.picoDst.root", -1)'
//   root4star -b -q 'readPicoDst.C("list_of_files.list", -1)'
//   root4star -b -q 'readPicoDst.C("input.picoDst.root", 500)'   // first 500 events only
//
// `input` can be:
//   - a single *.picoDst.root file
//   - a *.list / *.lis text file with one PicoDst path per line
// StPicoDstMaker handles both transparently (it builds a TChain internally).
// -----------------------------------------------------------------------------

class StChain;
StChain* chain = 0;

void readPicoDst(const Char_t* input = "picodst.list", Int_t nEvents = -1, const Char_t* outFile = "output.root", const char* badTowerList = "lists/bad_towers_pAu2015.list")
{
  // ---- load STAR shared libraries -----------------------------------------
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("libfastjet");
  gSystem->Load("libfastjettools");

  gSystem->Load("StRefMultCorr");

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StPicoAnalyzer");   // built from StRoot/StPicoAnalyzer by `cons`


  // ---- build the chain ------------------------------------------------------
  chain = new StChain();

  StPicoDstMaker* picoMaker =
      new StPicoDstMaker(StPicoDstMaker::IoRead, input, "picoDst");

  StPicoAnalyzer* analyzer =
    new StPicoAnalyzer("PicoAnalyzer", picoMaker, outFile);
  analyzer->SetBadTowerList(badTowerList);
  
  chain->Init();

  // ---- figure out how many events to run ------------------------------------
  Long64_t nEntries = picoMaker->chain()->GetEntries();
  if (nEvents < 0 || nEvents > nEntries) nEvents = nEntries;
  cout << "readPicoDst: " << nEntries << " events available, processing "
       << nEvents << endl;

  // ---- event loop -------------------------------------------------------------
  for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {
    chain->Clear();
    Int_t iret = chain->Make();
    if (iret) {
      cout << "Bad return code (" << iret << ") at event " << iEvent
           << " -- stopping." << endl;
      break;
    }
  }

  chain->Finish();
  delete chain;
}
