#!/usr/bin/env python3
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAT = re.compile(r"iter\s+(\d+)\s+loss=([0-9.eE+-]+)")


def main():
    if len(sys.argv) < 2:
        print("usage: plot_loss_from_nohup.py nohup.out [out.pdf]")
        return 2
    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("plot_loss_vs_iter.pdf")

    iters = []
    losses = []
    with infile.open("r", errors="ignore") as f:
        for line in f:
            m = PAT.search(line)
            if not m:
                continue
            iters.append(int(m.group(1)))
            losses.append(float(m.group(2)))

    if not iters:
        print("no loss values found (pattern: 'iter N loss=...')")
        return 1

    # sort by iteration in case nohup has interleaving
    data = sorted(zip(iters, losses))
    iters, losses = zip(*data)

    plt.figure(figsize=(7,5))
    plt.plot(iters, losses, marker='o', markersize=3, linewidth=1)
    plt.yscale('log')
    plt.xlabel('iteration')
    plt.ylabel('loss')
    plt.title('Loss vs Iteration')
    plt.grid(True, which='both', ls='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(outfile)
    print(f"wrote {outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
