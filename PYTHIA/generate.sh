#!/bin/bash
# Generate a PYTHIA pp-jet sample with ppjets_root.cc, split over parallel jobs.
#
# Usage:
#   PYTHIA/generate.sh <nEvents> <nJobs> <outDir> [firstSeed]
#
#   nEvents    total number of events (split evenly over the jobs)
#   nJobs      number of parallel jobs (= cores used)
#   outDir     output directory; created if missing, must not already contain .root files
#   firstSeed  random seed of job 1 (default 1); job i uses firstSeed + i - 1.
#              Use a different firstSeed when adding more statistics to an
#              existing sample, otherwise the new events repeat the old ones.
#
# Example (4M events on 16 cores):
#   PYTHIA/generate.sh 4000000 16 /path/to/MCOutput_4M
#
# Writes outDir/pythia_<i>.root for each job (logs in outDir/logs/). Point
# UNBIAS_INPUT_DIR at outDir: startBasis.C reads every .root file in it, so the
# per-job files need no merging.
#
# Needs (see SETUP.md): the conda env activated (ROOT, PYTHIA, $CONDA_PREFIX)
# and FJI set to the FastJet install prefix (or fastjet-config on PATH).
set -e

if [ $# -lt 3 ]; then
  sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
fi
NEVENTS=$1
NJOBS=$2
OUTDIR=$3
FIRSTSEED=${4:-1}

if [ -z "$CONDA_PREFIX" ]; then
  echo "generate.sh: CONDA_PREFIX is not set; activate the conda env first (see SETUP.md)." >&2
  exit 1
fi
if [ -n "$FJI" ]; then
  FASTJET_CONFIG=$FJI/bin/fastjet-config
else
  FASTJET_CONFIG=$(command -v fastjet-config || true)
fi
if [ ! -x "$FASTJET_CONFIG" ]; then
  echo "generate.sh: fastjet-config not found; set FJI to the FastJet install prefix." >&2
  exit 1
fi
FJLIB=$("$FASTJET_CONFIG" --prefix)/lib

mkdir -p "$OUTDIR/logs"
OUTDIR=$(cd "$OUTDIR" && pwd)
if ls "$OUTDIR"/*.root >/dev/null 2>&1; then
  echo "generate.sh: $OUTDIR already contains .root files; use an empty directory." >&2
  exit 1
fi

# Build once, next to the outputs (no .root extension, so startBasis.C ignores it).
SRC="$(cd "$(dirname "$0")" && pwd)/ppjets_root.cc"
EXE="$OUTDIR/logs/ppjets"
echo "Compiling $SRC"
c++ -std=c++17 -O2 "$SRC" -o "$EXE" \
    $("$FASTJET_CONFIG" --cxxflags --libs) -Wl,-rpath,"$FJLIB" \
    $(root-config --cflags --libs) \
    -I"$CONDA_PREFIX/include" -L"$CONDA_PREFIX/lib" -lpythia8 -Wl,-rpath,"$CONDA_PREFIX/lib"

export PYTHIA8DATA=${PYTHIA8DATA:-$CONDA_PREFIX/share/Pythia8/xmldoc}

# Split the events; the first (NEVENTS % NJOBS) jobs take one extra event.
PER=$((NEVENTS / NJOBS))
REST=$((NEVENTS % NJOBS))
echo "Running $NJOBS jobs, $NEVENTS events in total, seeds $FIRSTSEED..$((FIRSTSEED + NJOBS - 1))"
PIDS=()
for ((i = 1; i <= NJOBS; i++)); do
  N=$PER
  [ $i -le $REST ] && N=$((N + 1))
  SEED=$((FIRSTSEED + i - 1))
  "$EXE" "$N" "$OUTDIR/pythia_$i.root" "$SEED" > "$OUTDIR/logs/job_$i.log" 2>&1 &
  PIDS+=($!)
done

FAILED=0
for ((i = 1; i <= NJOBS; i++)); do
  if ! wait "${PIDS[$((i - 1))]}"; then
    echo "job $i failed, see $OUTDIR/logs/job_$i.log" >&2
    FAILED=1
  fi
done
[ $FAILED -ne 0 ] && exit 1

echo "Done. Output in $OUTDIR. Use it with:"
echo "  export UNBIAS_INPUT_DIR=$OUTDIR/"
