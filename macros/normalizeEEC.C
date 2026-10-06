
// Run this ONCE on the fully hadd'd output file -- after every SUMS job's
// raw histograms (numerator + weight accumulator) have been summed together
// by hadd, not before. Do not normalize per-job and then hadd the normalized
// histograms; hadd only correctly combines raw (unnormalized) sums.

#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include <iostream>

void DivideByBinWidth(TH1F* h)
{
    if (!h) return;
    for (Int_t b = 1; b <= h->GetNbinsX(); ++b) {
        const Double_t bw = h->GetXaxis()->GetBinWidth(b);
        if (bw <= 0.) continue;
        h->SetBinContent(b, h->GetBinContent(b) / bw);
        h->SetBinError(b, h->GetBinError(b) / bw);
    }
}

TH1F* NormalizeOne(TFile* inFile, TFile* outFile,
                    const TString& rawName, const TString& weightName)
{
    TH1F* hRaw = dynamic_cast<TH1F*>(inFile->Get(rawName));
    TH1F* hWeight = dynamic_cast<TH1F*>(inFile->Get(weightName));
 
    if (!hRaw) {
        std::cerr << "WARNING: could not find raw histogram \""
                   << rawName << "\" -- skipping." << std::endl;
        return nullptr;
    }
    if (!hWeight) {
        std::cerr << "WARNING: could not find weight histogram \""
                   << weightName << "\" -- skipping normalization of \""
                   << rawName << "\"." << std::endl;
        return nullptr;
    }
 
    const Double_t w = hWeight->GetBinContent(1);
    if (w <= 0.) {
        std::cerr << "WARNING: total weight for \"" << rawName
                   << "\" is non-positive (w = " << w
                   << ") -- skipping normalization." << std::endl;
        return nullptr;
    }
 
    outFile->cd();
 
    TH1F* hNorm = (TH1F*)hRaw->Clone(rawName + "Norm");
    hNorm->Scale(1. / w);
    DivideByBinWidth(hNorm);
 
    // Keep the raw title's selection text but relabel as normalized.
    TString title = hRaw->GetTitle();
    title.ReplaceAll("Raw EEC numerator", "Normalized EEC");
    hNorm->SetTitle(title);
 
    std::cout << "Normalized \"" << rawName << "\" (w = " << w
               << ", Integral(\"width\") = " << hNorm->Integral("width")
               << ")" << std::endl;
 
    return hNorm;
}
 
void normalizeEEC(const char* inputFile  = "merged.root",
                   const char* outputFile = "merged_normalized.root")
{
    TFile* inFile = TFile::Open(inputFile, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "ERROR: could not open input file \"" << inputFile
                   << "\"" << std::endl;
        return;
    }
 
    TFile* outFile = TFile::Open(outputFile, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "ERROR: could not open output file \"" << outputFile
                   << "\"" << std::endl;
        inFile->Close();
        return;
    }
 
    // One entry per pT,jet selection: {raw numerator name, weight accumulator name}
    /*const std::pair<TString, TString> selections[] = {
        {"hEEC_low",  "hEECWeight_low"},
        {"hEEC_mid",  "hEECWeight_mid"},
        {"hEEC_high", "hEECWeight_high"},
	};*/
    const TString rawNames[]    = {"hEEC_low",  "hEEC_mid",  "hEEC_high"};
    const TString weightNames[] = {"hEECWeight_low", "hEECWeight_mid", "hEECWeight_high"};
    const Int_t nSelections = 3;
 
    for (Int_t i = 0; i < nSelections; ++i) {
      TH1F* hNorm = NormalizeOne(inFile, outFile, rawNames[i], weightNames[i]);
      if (hNorm) hNorm->Write();
    }
    /*
    for (const auto& sel : selections) {
        TH1F* hNorm = NormalizeOne(inFile, outFile, sel.first, sel.second);
        if (hNorm) hNorm->Write();
    }
    */
    // Also carry the raw (merged) histograms into the output file for
    // convenience/cross-checking, unmodified.
    /*for (const auto& sel : selections) {
        TH1F* hRaw = dynamic_cast<TH1F*>(inFile->Get(sel.first));
        TH1F* hWeight = dynamic_cast<TH1F*>(inFile->Get(sel.second));
        if (hRaw) {
            outFile->cd();
            hRaw->Clone(sel.first)->Write();
        }
        if (hWeight) {
            outFile->cd();
            hWeight->Clone(sel.second)->Write();
        }
    }
    */

    for (Int_t i = 0; i < nSelections; ++i) {
      TH1F* hRaw = dynamic_cast<TH1F*>(inFile->Get(rawNames[i]));
      TH1F* hWeight = dynamic_cast<TH1F*>(inFile->Get(weightNames[i]));
      if (hRaw) {
	outFile->cd();
	hRaw->Clone(rawNames[i])->Write();
      }
      if (hWeight) {
	outFile->cd();
	hWeight->Clone(weightNames[i])->Write();
      }
    }
    
    outFile->Write();
    outFile->Close();
    inFile->Close();
 
    std::cout << "Done. Normalized output written to \"" << outputFile
               << "\"" << std::endl;
}
