"""
Diagnosi formale del bias sistematico di posizione del ring attractor a rate
(osservato ripetutamente nelle sezioni 4.4.2, 5.1, 5.2 senza una spiegazione
verificata alla radice -- la sezione 9 lo elenca esplicitamente come sviluppo
futuro: "determinare se sia un artefatto puramente numerico o rifletta una
proprieta' reale, per quanto non desiderata, del modello discretizzato").

Usa la stessa equazione del modello verificato in ring_attractor_rate.py
(kernel eccitatorio locale J_exc + inibizione globale proporzionale alla
media della rete + attivazione rettificata satura a r_max) -- non il
modello a "cappello messicano" (due kernel Gaussiani opposti) usato altrove
in questo file in una versione precedente: durante la preparazione di
questa diagnosi si e' scoperto che il file ring_attractor_rate.py in cima
al repository era stato accidentalmente sovrascritto con quella variante
diversa, che non riproduce piu' i risultati documentati nella sezione 4.4.2
(si veda la sezione 9.1 per i dettagli del bug e della correzione).

Due test indipendenti, entrambi falsificabili:

PARTE 1 -- scaling con la risoluzione della griglia N. Se il bias e' un
artefatto della somma di Riemann che approssima l'integrale di convoluzione
continuo (rec_input = J @ r / N), deve ridursi sistematicamente all'aumentare
di N. Se invece e' una proprieta' reale della dinamica del kernel a cappello
messicano, deve restare approssimativamente costante al variare di N.
Lo stimolo e' posizionato deliberatamente a una coordinata (0.2537) che non
coincide MAI esattamente con un punto della griglia per nessuno degli N
testati, cosi' da non confondere lo scaling reale con un banale "aggancio"
del picco al punto di griglia piu' vicino quando N e' piccolo.

PARTE 2 -- dipendenza dalla posizione sub-griglia a N fissa. Se il bias
dipende dalla posizione dello stimolo RELATIVA al punto di griglia piu'
vicino (un segno diretto di artefatto di discretizzazione), deve mostrare
una struttura periodica con periodo 1/N quando si fa scorrere lo stimolo
con continuita' attraverso alcune celle della griglia.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(0)


def circ_dist(d):
    return 0.5 - np.abs(np.abs(d) - 0.5)


def signed_circ_diff(a, b):
    """Differenza angolare con segno a-b, in (-0.5, 0.5]."""
    return ((a - b + 0.5) % 1.0) - 0.5


def run_ring(N, stim_center, T_total=400, dt=0.2, tau=10.0,
             sigma_exc=0.05, A_exc=20.0, global_inhib=5.0, r_max=5.0,
             stim_strength=3.0, stim_duration=80, stim_sigma=0.03):
    x = np.arange(N) / N
    d_matrix = circ_dist(x[:, None] - x[None, :])
    J_exc = A_exc*np.exp(-d_matrix**2/(2*sigma_exc**2))

    def F(u):
        return np.clip(u, 0, r_max)

    n_steps = int(T_total/dt)
    stim_duration_steps = int(stim_duration/dt)
    stim_profile = stim_strength*np.exp(-circ_dist(x-stim_center)**2/(2*stim_sigma**2))

    r = np.zeros(N)
    for step in range(n_steps):
        I_ext = stim_profile if step < stim_duration_steps else 0.0
        rec_input = J_exc @ r / N - global_inhib*r.mean()
        dr = (-r + F(rec_input + I_ext)) / tau
        r = r + dt*dr

    peak_idx = np.argmax(r)
    return x[peak_idx], r[peak_idx], x


# ============================================================
# PARTE 1: scaling del bias con N (risoluzione della griglia)
# ============================================================
print("=" * 70)
print("PARTE 1: bias di posizione in funzione della risoluzione N")
print("=" * 70)

STIM_CENTER = 0.2537  # non allineato alla griglia per nessun N testato
Ns = [50, 100, 200, 400, 800, 1600]
biases_N = []
peaks_N = []

for N in Ns:
    peak_pos, peak_val, _ = run_ring(N, STIM_CENTER)
    bias = signed_circ_diff(peak_pos, STIM_CENTER)
    biases_N.append(bias)
    peaks_N.append(peak_val)
    print(f"N={N:5d}  griglia=1/N={1/N:.5f}  picco a x={peak_pos:.5f}  "
          f"bias={bias:+.5f}  |bias|/griglia={abs(bias)/(1/N):.2f}  "
          f"ampiezza picco={peak_val:.3f}")

biases_N = np.array(biases_N)
peaks_N = np.array(peaks_N)

# Se il bias fosse un puro artefatto di quantizzazione della lettura
# (argmax sulla griglia), scalerebbe come 1/N. Verifichiamo il rapporto.
ratio_first_last = abs(biases_N[0]) / abs(biases_N[-1]) if biases_N[-1] != 0 else float('inf')
N_ratio = Ns[-1] / Ns[0]
print(f"\nRapporto |bias| tra N={Ns[0]} e N={Ns[-1]}: {ratio_first_last:.2f}x "
      f"(un artefatto di quantizzazione 1/N implicherebbe un fattore ~{N_ratio}x)")

# ============================================================
# PARTE 2: dipendenza dalla posizione sub-griglia, N fissa
# ============================================================
print()
print("=" * 70)
print("PARTE 2: bias in funzione della posizione sub-griglia (N=200 fisso)")
print("=" * 70)

N_FIXED = 200
grid_step = 1.0 / N_FIXED
# scorre lo stimolo con continuita' attraverso 3 celle di griglia consecutive
n_sweep = 31
base_center = 0.30  # punto di partenza (coincide con un nodo di griglia)
offsets = np.linspace(0, 3*grid_step, n_sweep)
biases_sub = []

for off in offsets:
    center = base_center + off
    peak_pos, peak_val, _ = run_ring(N_FIXED, center)
    bias = signed_circ_diff(peak_pos, center)
    biases_sub.append(bias)

biases_sub = np.array(biases_sub)
# offset relativo entro UNA cella di griglia (modulo grid_step), per
# verificare se il bias ha davvero periodo 1/N
offsets_mod = offsets % grid_step

print(f"Deviazione standard del bias lungo lo sweep: {biases_sub.std():.5f} "
      f"(su una cella di griglia larga {grid_step:.5f})")
print(f"Bias minimo/massimo nello sweep: {biases_sub.min():+.5f} / {biases_sub.max():+.5f}")

# correlazione tra bias e posizione entro la cella di griglia
corr = np.corrcoef(offsets_mod, biases_sub)[0, 1]
print(f"Correlazione tra bias e posizione entro la cella di griglia: r={corr:.3f}")

# ============================================================
# PARTE 3: dipendenza dal passo di integrazione temporale dt
# ============================================================
print()
print("=" * 70)
print("PARTE 3: bias in funzione del passo di integrazione dt (N=200 fisso)")
print("=" * 70)
print("Se il bias fosse un artefatto di troncamento di Eulero esplicito,")
print("dovrebbe ridursi al diminuire di dt (l'errore locale di Eulero scala")
print("come dt).")

dts = [0.2, 0.1, 0.05, 0.025]
biases_dt = []
for dt_test in dts:
    peak_pos, peak_val, _ = run_ring(200, STIM_CENTER, dt=dt_test)
    bias = signed_circ_diff(peak_pos, STIM_CENTER)
    biases_dt.append(bias)
    print(f"dt={dt_test:.3f}  picco a x={peak_pos:.5f}  bias={bias:+.5f}  "
          f"ampiezza picco={peak_val:.3f}")

biases_dt = np.array(biases_dt)
ratio_dt = abs(biases_dt[0]) / abs(biases_dt[-1]) if biases_dt[-1] != 0 else float('inf')
print(f"\nRapporto |bias| tra dt={dts[0]} e dt={dts[-1]}: {ratio_dt:.2f}x "
      f"(un artefatto di troncamento di Eulero, lineare in dt, implicherebbe "
      f"un fattore ~{dts[0]/dts[-1]:.1f}x)")

# ============================================================
# Figura riassuntiva
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.2))

axes[0].plot(Ns, np.abs(biases_N), 'o-', color='C0', label='|bias| osservato')
# linea di riferimento: scaling 1/N normalizzato al primo punto
ref = np.abs(biases_N[0]) * (Ns[0] / np.array(Ns))
axes[0].plot(Ns, ref, '--', color='gray', label='scaling atteso se 1/N (artefatto)')
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_xlabel('N (punti di griglia)')
axes[0].set_ylabel('|bias di posizione|')
axes[0].set_title('Parte 1: bias vs risoluzione della griglia')
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].plot(offsets, biases_sub, 'o-', color='C1')
for k in range(4):
    axes[1].axvline(k*grid_step, color='cyan', linestyle='--', linewidth=0.8, alpha=0.6)
axes[1].set_xlabel('Posizione dello stimolo (offset da x=0.30)')
axes[1].set_ylabel('Bias di posizione (con segno)')
axes[1].set_title(f'Parte 2: bias vs posizione sub-griglia (N={N_FIXED})\nlinee tratteggiate = nodi di griglia')
axes[1].grid(alpha=0.3)

axes[2].plot(dts, np.abs(biases_dt), 'o-', color='C2', label='|bias| osservato')
ref_dt = np.abs(biases_dt[0]) * (np.array(dts) / dts[0])
axes[2].plot(dts, ref_dt, '--', color='gray', label='scaling atteso se lineare in dt\n(artefatto di Eulero)')
axes[2].set_xscale('log')
axes[2].set_yscale('log')
axes[2].set_xlabel('dt (ms)')
axes[2].set_ylabel('|bias di posizione|')
axes[2].set_title('Parte 3: bias vs passo di integrazione dt')
axes[2].legend(fontsize=8)
axes[2].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('diagnosi_bias_posizione_scaling.png', dpi=150)
print("\nGrafico salvato in diagnosi_bias_posizione_scaling.png")
