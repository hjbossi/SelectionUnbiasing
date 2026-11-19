// driveBins.C
// Run startBasis.C once per pT bin, each writing its own output ROOT file,
// shared TEnv config file, and EEC terms CSV -- so the subjet radii and the
// EEC term grid can be set/edited independently per bin (see
// basis_functions.h / startBasis.C).
//
// For now both bins below use the same subjet radii ("0.10,0.20") and point
// eecTermsFile at the same shared "eec_terms_example.csv" -- a hand-written
// CSV (checked in alongside this macro) that reproduces the old compiled
// default_eec_terms() grid exactly. Unlike leaving eecTermsFile empty (which
// would make startBasis.C auto-generate a fresh, separate default file per
// bin), both bins here read the SAME file, so hand-editing
// eec_terms_example.csv later changes both bins' term grids at once. Point
// a bin at its own copy instead (and edit that copy) once you want the two
// bins to diverge.
//
// pT-DEPENDENT EEC ANGULAR CUTOFFS: the E column of the EEC terms CSV (and
// unbias_weights.cc's --dR-min) are REFERENCE cutoffs at pT = 100 GeV.
// unbias_weights.cc applies, per bin,
//     E_used      = E      * 100 / BinCenter
//     dR_min_used = dR_min * 100 / BinCenter
// where BinCenter is the value this macro writes for that bin
// ("Unbiasing.BinCenter": the pT-window midpoint when binCenter <= 0, or
// binCenter itself when > 0). Subjet radii are NOT scaled. With the
// midpoint defaults below: bin 100-140 (center 120) -> scale 0.8333,
// bin 140-180 (center 160) -> scale 0.625.
//
// Usage (run from the AnalysisCode/ directory -- all paths are relative to it):
//   root -l -b -q 'driveBins.C("/path/to/MCOutput/")'
// or set the input directory once via the environment:
//   export UNBIAS_INPUT_DIR=/path/to/MCOutput/
//   root -l -b -q driveBins.C
// The input directory is scanned recursively for the ROOT files holding the
// jet tree (every .root file in it is read).
//
// Produces, per bin (in ../Output/files/, plots in ../Output/plots/):
//   startBasis_<bin>.root                    -- biased/reference trees
//   unbiasing_config_<bin>.env                -- BinCenter/PtLow/PtHigh/
//                                                 SubjetRadii/EECTermsFile
//
// Reads (shared by both bins):
//   eec_terms_example.csv                     -- the EEC (A,B,m,E) term grid
//                                                (E = reference cut at 100 GeV)

// Pull in startBasis.C's definition directly (rather than loading it at
// runtime via gROOT->ProcessLine(".L startBasis.C")) so the call below
// resolves at compile time. Cling JIT-compiles the whole driveBins() body
// before executing it, so a ProcessLine(".L ...") inside that body loads
// too late for the compiler to see -- it reports startBasis as an
// undeclared identifier even though the load would happen fine at runtime.
#include "startBasis.C"

void driveBins(const char* inputDirArg = "") {
  // Input directory: macro argument, else $UNBIAS_INPUT_DIR, else abort
  // (a wrong/empty directory would otherwise just give 0 entries).
  const char* envDir = gSystem->Getenv("UNBIAS_INPUT_DIR");
  const char* inputDir = (inputDirArg && inputDirArg[0]) ? inputDirArg
                         : (envDir && envDir[0])         ? envDir
                         : "";
  if (!inputDir[0]) {
    std::cerr << "driveBins: no input directory. Pass it as an argument, e.g.\n"
              << "  root -l -b -q 'driveBins.C(\"/path/to/MCOutput/\")'\n"
              << "or set UNBIAS_INPUT_DIR." << std::endl;
    return;
  }
  const double p = 8.0;
  const double d = 0.1;
  const unsigned rngSeed = 12345;
  const double binCenter = -1.0;              // derive from the pT window (midpoint); also sets the EEC angular-cutoff scale 100/BinCenter
  const char* subjetRadii = "0.10,0.20";              // same subjet radii for both bins, for now
  const char* eecTermsFile = "eec_terms_example.csv"; // shared EEC term grid for both bins, for now

  // Generated files go under ../Output/ (next to AnalysisCode/, git-ignored):
  // ROOT/config files in Output/files/, plots in Output/plots/.
  const TString filesDir = "../Output/files";
  const TString plotsDir = "../Output/plots";
  gSystem->mkdir(filesDir, kTRUE);
  gSystem->mkdir(plotsDir, kTRUE);

  // Sample tag (required), e.g. UNBIAS_TAG=10M -> startBasis_100_140_10M.root,
  // so outputs made from different input samples never overwrite each other.
  const char* envTag = gSystem->Getenv("UNBIAS_TAG");
  if (!(envTag && envTag[0])) {
    std::cerr << "driveBins: UNBIAS_TAG is not set. Name the input sample, e.g.\n"
              << "  export UNBIAS_TAG=10M" << std::endl;
    return;
  }
  const TString tag = TString("_") + envTag;

  // ---- Bin 1: 100-140 GeV ----
  startBasis(inputDir,
             filesDir + "/startBasis_100_140" + tag + ".root",
             100, 140,
             p, d, rngSeed,
             binCenter,
             filesDir + "/unbiasing_config_100_140" + tag + ".env",
             subjetRadii,
             eecTermsFile,
             unbiasRootPlotDir(plotsDir));
  unbiasMplStartBasisPlot(filesDir + "/startBasis_100_140" + tag + ".root", plotsDir);

  // ---- Bin 2: 140-180 GeV ----
  startBasis(inputDir,
             filesDir + "/startBasis_140_180" + tag + ".root",
             140, 180,
             p, d, rngSeed,
             binCenter,
             filesDir + "/unbiasing_config_140_180" + tag + ".env",
             subjetRadii,
             eecTermsFile,
             unbiasRootPlotDir(plotsDir));
  unbiasMplStartBasisPlot(filesDir + "/startBasis_140_180" + tag + ".root", plotsDir);
}