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
// Usage:
//   root -l -b -q driveBins.C
//
// Produces, per bin:
//   startBasis_<bin>.root                    -- biased/reference trees
//   unbiasing_config_<bin>.env                -- BinCenter/PtLow/PtHigh/
//                                                 SubjetRadii/EECTermsFile
//
// Reads (shared by both bins):
//   eec_terms_example.csv                     -- the EEC (A,B,m,E) term grid

// Pull in startBasis.C's definition directly (rather than loading it at
// runtime via gROOT->ProcessLine(".L startBasis.C")) so the call below
// resolves at compile time. Cling JIT-compiles the whole driveBins() body
// before executing it, so a ProcessLine(".L ...") inside that body loads
// too late for the compiler to see -- it reports startBasis as an
// undeclared identifier even though the load would happen fine at runtime.
#include "startBasis.C"

void driveBins() {
  const char* inputDir = "/home/hbossi/SelectionUnbiasing/MCOutput/v3/";
  const double p = 8.0;
  const double d = 0.1;
  const unsigned rngSeed = 12345;
  const double binCenter = -1.0;              // derive from the pT window (midpoint)
  const char* subjetRadii = "0.10,0.20";              // same subjet radii for both bins, for now
  const char* eecTermsFile = "eec_terms_example.csv"; // shared EEC term grid for both bins, for now

  // ---- Bin 1: 100-140 GeV ----
  startBasis(inputDir,
             "startBasis_100_140.root",
             100, 140,
             p, d, rngSeed,
             binCenter,
             "unbiasing_config_100_140.env",
             subjetRadii,
             eecTermsFile);

  // ---- Bin 2: 140-180 GeV ----
  startBasis(inputDir,
             "startBasis_140_180.root",
             140, 180,
             p, d, rngSeed,
             binCenter,
             "unbiasing_config_140_180.env",
             subjetRadii,
             eecTermsFile);
}