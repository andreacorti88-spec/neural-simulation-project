"""
Sezione 5.46: coda di frammentazione oltre il picco di Bragg del
carbonio-12, confrontata con il protone alla stessa profondita' relativa.

Legge i profili di dose in profondita' prodotti da braggProfile (fascio
largo 5x5cm, fette da 1mm -- statistica enormemente migliore della
mappa per-neurone di cnaoRingImpact per questo effetto, che a profondita'
oltre il picco e' raro) per carbonio-12 120 MeV/u e protone 70 MeV, e
confronta la dose residua oltre ciascun picco (in rapporto al proprio
picco) tra le due particelle.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_profile(path):
    with open(path) as f:
        lines = [l for l in f if not l.startswith('#')]
    r = csv.DictReader(lines)
    depths, doses = [], []
    for row in r:
        depths.append(float(row['depth_mm']))
        doses.append(float(row['sum_edep_MeV']))
    return np.array(depths), np.array(doses)


d_c, dose_c = load_profile('../braggProfile/build/bragg_carbon120_peak_and_tail.csv')
d_p, dose_p = load_profile('../braggProfile/build/bragg_proton70_peak_and_tail.csv')

peak_c_idx = np.argmax(dose_c)
peak_p_idx = np.argmax(dose_p)
peak_c_depth, peak_c_dose = d_c[peak_c_idx], dose_c[peak_c_idx]
peak_p_depth, peak_p_dose = d_p[peak_p_idx], dose_p[peak_p_idx]

print(f"Picco carbonio-12 120 MeV/u: {peak_c_depth:.1f}mm, dose={peak_c_dose:.3e} MeV")
print(f"Picco protone 70 MeV: {peak_p_depth:.1f}mm, dose={peak_p_dose:.3e} MeV")

ratio_c = dose_c / peak_c_dose
ratio_p = dose_p / peak_p_dose

print()
print("Rapporto dose/picco a distanze fisse OLTRE il proprio picco:")
print(f"{'margine oltre il picco (mm)':>28} {'carbonio-12':>14} {'protone':>14} {'rapporto C/p':>14}")
for margin in [5, 10, 15, 20, 30, 50]:
    tc = peak_c_depth + margin
    tp = peak_p_depth + margin
    ic = np.argmin(np.abs(d_c - tc))
    ip = np.argmin(np.abs(d_p - tp))
    rc, rp = ratio_c[ic], ratio_p[ip]
    print(f"{margin:>28} {rc:>14.6f} {rp:>14.8f} {(rc/rp if rp>0 else float('inf')):>14.1f}")

# ------------------------------------------------------------
# Figura: dose relativa al picco (scala log) vs profondita', per
# entrambe le particelle, allineate al proprio picco (asse x = margine
# oltre il picco in mm)
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 5))
mask_c = d_c >= peak_c_depth
mask_p = d_p >= peak_p_depth
ax.semilogy(d_c[mask_c] - peak_c_depth, np.clip(ratio_c[mask_c], 1e-7, None),
            label='Carbonio-12, 120 MeV/u (picco a 34.5mm)', color='C0')
ax.semilogy(d_p[mask_p] - peak_p_depth, np.clip(ratio_p[mask_p], 1e-7, None),
            label='Protone, 70 MeV (picco a 39.5mm)', color='C1')
ax.set_xlabel('Profondita\' oltre il proprio picco di Bragg (mm)')
ax.set_ylabel('Dose relativa al picco (scala log)')
ax.set_title('Coda di frammentazione: dose oltre il picco, normalizzata al picco stesso')
ax.set_xlim(0, 100)
ax.legend()
ax.grid(alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('fragmentation_tail.png', dpi=150)
print("\nGrafico salvato in fragmentation_tail.png")
