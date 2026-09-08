# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
CATENA DI DUE ATTRATTORI: l'attivazione dell'uno determina l'altro
====================================================================
Estensione del ring attractor a due popolazioni collegate:
  RING A: riceve lo stimolo esterno diretto, si comporta esattamente
          come il ring attractor singolo gia' verificato

  RING B: NON riceve alcuno stimolo esterno diretto. L'unico input che
          riceve e' una proiezione della posizione del bump di A,
          spostata di un offset fisso (0.3 sull'anello) -- una mappa
          associativa fissa, tipo "se A si stabilizza qui, B tende a
          stabilizzarsi la'"

Questo e' il primo passo meccanico verso una "catena di attivazioni":
stimolare A -> A si stabilizza -> A "accende" B in una posizione
corrispondente, SENZA che nessuno dica esplicitamente a B dove andare.

Cosa questo NON e': non c'e' alcuna forma di ragionamento, scelta, o
valutazione di alternative. E' pura propagazione causale attraverso
uno stato intermedio stabile -- ma e' esattamente l'ingrediente di
base necessario prima di poter costruire qualcosa di piu' sofisticato
(es. piu' stati concatenati, competizione tra piu' possibili "B",
un segnale di rinforzo che sceglie quale mappa usare).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def circ_dist(d):
    return 0.5 - np.abs(np.abs(d) - 0.5)



N = 200
x = np.arange(N) / N
tau = 10.0
dt = 0.2
T_total = 500
n_steps = int(T_total/dt)

A_exc = 20.0
sigma_exc = 0.05
global_inhib = 5.0
r_max = 5.0
stim_strength = 3.0
stim_width = 0.03
stim_duration_ms = 80

d_matrix = circ_dist(x[:, None] - x[None, :])
J_exc = A_exc * np.exp(-d_matrix**2 / (2*sigma_exc**2))

OFFSET_AB = 0.3           # spostamento fisso della mappa associativa A->B
COUPLING_STRENGTH = 1.3    # intensita' dell'influenza di A su B



def F(u):
    return np.clip(u, 0, r_max)



def run_chained(stim_center_A, apply_stim=True):
    rA = np.zeros(N)
    rB = np.zeros(N)
    stim_duration_steps = int(stim_duration_ms/dt)
    stim_profile = (stim_strength * np.exp(-circ_dist(x-stim_center_A)**2/(2*stim_width**2))
                     if apply_stim else np.zeros(N))
   history_A = np.zeros((n_steps, N))
   history_B = np.zeros((n_steps, N))

   for step in range(n_steps):
       I_ext_A = stim_profile if (apply_stim and step < stim_duration_steps) else 0.0

        rec_A = J_exc @ rA / N - global_inhib*rA.mean()
        drA = (-rA + F(rec_A + I_ext_A)) / tau
        rA_new = rA + dt*drA

        # Ring B: nessuno stimolo esterno -- solo la proiezione spostata di A
        cross_input = np.zeros(N)
        if rA.max() > 0.5:
            peak_A = x[np.argmax(rA)]
            target_B = (peak_A + OFFSET_AB) % 1.0
            cross_input = COUPLING_STRENGTH * rA.max() * np.exp(
                -circ_dist(x-target_B)**2/(2*stim_width**2))

        rec_B = J_exc @ rB / N - global_inhib*rB.mean()
        drB = (-rB + F(rec_B + cross_input)) / tau
        rB_new = rB + dt*drB

        rA, rB = rA_new, rB_new
        history_A[step] = rA
        history_B[step] = rB

   return rA, rB, history_A, history_B



if __name__ == '__main__':
    print("=" * 65)
    print("TEST 1: A stimolato in posizioni diverse, B riceve SOLO l'influenza da A")
    print("=" * 65)
    test_centers = [0.1, 0.3, 0.5, 0.7]
    for center_A in test_centers:
        rA_final, rB_final, _, _ = run_chained(center_A)
        posA = x[np.argmax(rA_final)]
        posB = x[np.argmax(rB_final)] if rB_final.max() > 0.5 else None
        posB_str = f"{posB:.3f}" if posB is not None else "NESSUN BUMP"
        expected_B = (center_A + OFFSET_AB) % 1.0
        print(f" Stimolo A a x={center_A:.2f} -> A: x={posA:.3f}, "
