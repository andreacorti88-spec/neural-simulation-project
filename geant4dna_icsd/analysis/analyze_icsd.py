#!/usr/bin/env python3
import sys, csv, math
from pathlib import Path

PAPER_M1 = {
    2: {20:0.77,50:2.21,100:4.30,300:5.35,600:2.90,1000:1.83,5000:7.10,10000:4.02},
    4: {20:0.44,50:1.77,100:3.65,300:4.55,600:2.54,1000:1.67,5000:6.56,10000:3.97},
    6: {20:0.83,50:2.43,100:4.83,300:6.65,600:3.59,1000:2.26,5000:8.09,10000:4.53},
}
PAPER_F2 = {
    2: {50:0.89,100:0.96,300:0.81,600:0.59,1000:0.42,5000:0.90,10000:0.69},
    4: {50:0.69,100:0.95,300:0.78,600:0.56,1000:0.39,5000:0.90,10000:0.71},
    6: {50:0.93,100:0.99,300:0.89,600:0.68,1000:0.50,5000:0.95,10000:0.77},
}

def load_counts(path):
    return [int(l.strip()) for l in open(path) if l.strip() and not l.startswith("#")]

def m1_f2_f3(counts):
    n = len(counts)
    mean = sum(counts)/n
    f2 = sum(1 for c in counts if c>=2)/n
    return mean, f2, n

def main():
    if len(sys.argv) != 3:
        print("Uso: python3 analyze_icsd.py <cartella_risultati> <opzione: 2|4|6>")
        sys.exit(1)
    result_dir = Path(sys.argv[1])
    option = int(sys.argv[2])
    files = sorted(result_dir.glob("icsd_*eV.csv"), key=lambda p:int(p.stem.split("_")[1].replace("eV","")))
    print(f"{'E (eV)':>8} {'N':>8} {'M1 sim':>9} {'M1 paper':>9} {'diff%':>8} {'F2 sim':>9} {'F2 paper':>9} {'diff%':>8}")
    print("-"*80)
    for f in files:
        energy = int(f.stem.split("_")[1].replace("eV",""))
        counts = load_counts(f)
        m1, f2, n = m1_f2_f3(counts)
        pm1 = PAPER_M1[option].get(energy)
        pf2 = PAPER_F2[option].get(energy)
        dm1 = 100*(m1-pm1)/pm1 if pm1 else float("nan")
        df2 = 100*(f2-pf2)/pf2 if pf2 else float("nan")
        print(f"{energy:>8} {n:>8} {m1:>9.3f} {pm1 or float('nan'):>9.3f} {dm1:>7.1f}% {f2:>9.3f} {pf2 or float('nan'):>9.3f} {df2:>7.1f}%")

if __name__ == "__main__":
    main()
