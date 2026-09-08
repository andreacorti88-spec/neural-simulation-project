# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
POLITICA CONTESTUALE: risposte diverse per contesti diversi
====================================================================
Finora il sistema imparava UN'UNICA associazione corretta (bandito a
3 braccia semplice): sempre la stessa risposta, indipendentemente da
dove si stimolasse A.

Qui il compito e' piu' difficile e piu' realistico: la posizione dello
stimolo su A definisce un CONTESTO (3 contesti possibili), e ciascun
contesto ha una propria risposta corretta diversa dagli altri:

    Contesto 1 (A stimolato vicino a x=0.1)             -> risposta corretta: candidato 1
    Contesto 2 (A stimolato vicino a x=0.5)             -> risposta corretta: candidato 2
    Contesto 3 (A stimolato vicino a x=0.8)             -> risposta corretta: candidato 3

Il sistema deve imparare TRE associazioni contemporaneamente, non una
sola -- e a fine training va verificato che abbia davvero imparato a
DISTINGUERE i contesti (non solo a preferire un candidato in generale).

Meccanismo: invece di un unico vettore di pesi, manteniamo una TABELLA
pesi[contesto][candidato] -- una forma elementare di "memoria associativa
condizionata al contesto", concettualmente il primo passo verso una
politica decisionale vera e propria (che sceglie l'azione in base allo
stato, non sempre la stessa azione).
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

CANDIDATE_OFFSETS = [0.15, 0.45, 0.75]
COUPLING_STRENGTH = 1.3
NOISE_STD = 0.15



def F(u):
    return np.clip(u, 0, r_max)



def run_trial(stim_center_A, weights, seed_noise):
    rA = np.zeros(N)
    rB = np.zeros(N)
    stim_duration_steps = int(stim_duration_ms/dt)
    stim_profile = stim_strength * np.exp(-circ_dist(x-stim_center_A)**2/(2*stim_width**2))
    rng = np.random.RandomState(seed_noise)

      for step in range(n_steps):
          I_ext_A = stim_profile if step < stim_duration_steps else 0.0
          rec_A = J_exc @ rA / N - global_inhib*rA.mean()
          rA = rA + dt*(-rA + F(rec_A + I_ext_A)) / tau

          cross_input = np.zeros(N)
          if rA.max() > 0.5:
              peak_A = x[np.argmax(rA)]
              for offset, w in zip(CANDIDATE_OFFSETS, weights):
                  target = (peak_A + offset) % 1.0
                  cross_input += w * COUPLING_STRENGTH * rA.max() * np.exp(
                      -circ_dist(x-target)**2/(2*stim_width**2))
              cross_input += rng.normal(0, NOISE_STD, N)

          rec_B = J_exc @ rB / N - global_inhib*rB.mean()
          rB = rB + dt*(-rB + F(rec_B + cross_input)) / tau
    if rB.max() < 0.5:
        return None
    winner_pos = x[np.argmax(rB)]
    dists = [circ_dist(winner_pos - ((stim_center_A+off) % 1.0)) for off in
CANDIDATE_OFFSETS]
    return int(np.argmin(dists))



# ---------------------------------------------------------------
# TRE CONTESTI, TRE RISPOSTE CORRETTE DIVERSE
# ---------------------------------------------------------------
CONTEXTS = [0.1, 0.5, 0.8]           # posizioni di stimolo per A
CORRECT_PER_CONTEXT = [0, 1, 2]       # risposta corretta per ciascun contesto
LEARNING_RATE = 0.15
N_TRIALS = 150   # piu' prove: il compito e' piu' difficile (3 associazioni)

# tabella pesi[contesto][candidato], tutti uguali all'inizio
weights_table = np.full((3, 3), 0.5)
history = np.zeros((N_TRIALS, 3, 3))
choices_per_context = {0: [], 1: [], 2: []}

rng_context = np.random.RandomState(0)

print("Contesti e risposte corrette (nascoste al sistema):")
for c, correct in enumerate(CORRECT_PER_CONTEXT):
    print(f" Contesto {c+1} (stimolo a x={CONTEXTS[c]}) -> risposta corretta: candidato
{correct+1}")
print()

for trial in range(N_TRIALS):
    context = rng_context.randint(0, 3)    # contesto scelto a caso ad ogni prova
    stim_pos = CONTEXTS[context]
    correct = CORRECT_PER_CONTEXT[context]

    winner = run_trial(stim_pos, weights_table[context], seed_noise=trial+1)
    if winner is not None:
         reward = 1.0 if winner == correct else 0.0
         weights_table[context, winner] += LEARNING_RATE * (reward - weights_table[context,
winner])
         weights_table = np.clip(weights_table, 0.05, 2.0)
         choices_per_context[context].append(1 if winner == correct else 0)

    history[trial] = weights_table.copy()

    if trial % 30 == 0 or trial == N_TRIALS-1:
        print(f"Prova {trial+1:3d} (contesto {context+1}): "
              f"pesi contesto {context+1} ora = {np.round(weights_table[context],3)}")

# ---------------------------------------------------------------
# VALUTAZIONE FINALE: il sistema ha imparato a distinguere i 3 contesti?
# ---------------------------------------------------------------
print("\n" + "=" * 65)
print("VALUTAZIONE FINALE (pesi appresi, test pulito senza aggiornamento)")
print("=" * 65)
for c in range(3):
    test_wins = []
