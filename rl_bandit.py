"""
APPRENDIMENTO PER RINFORZO: la rete impara QUALE alternativa scegliere
====================================================================
Estensione finale: i pesi delle tre mappe candidate A->B non sono piu'
fissi. Ad ogni prova:
  1. A viene stimolato, B compete tra i 3 candidati e sceglie un vincitore
  2. Viene dato un segnale di RICOMPENSA: +1 se il vincitore e' quello
     "corretto" (definito da noi, come una regola esterna fissa -- il
     sistema non sa a priori quale sia), 0 altrimenti
  3. Il peso del candidato SCELTO viene aggiornato verso la ricompensa
     ricevuta (regola tipo Rescorla-Wagner / bandito multi-braccio):
         peso_scelto += tasso_apprendimento * (ricompensa - peso_scelto)

Nessun altro peso viene toccato -- il sistema impara SOLO dalle proprie
scelte, non da un confronto esplicito con l'alternativa corretta.

Questo e' un vero, piccolo ciclo di reinforcement learning: prova,
osserva il risultato, aggiorna il comportamento futuro. E' il primo
esperimento del progetto in cui il sistema non si limita a reagire o
a propagare, ma MODIFICA IL PROPRIO COMPORTAMENTO in base a un
risultato.
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
NOISE_STD = 0.15   # rumore continuo -- necessario per "esplorare" le alternative


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
    dists = [circ_dist(winner_pos - ((stim_center_A+off) % 1.0)) for off in CANDIDATE_OFFSETS]
    return int(np.argmin(dists))


# ---------------------------------------------------------------
# CICLO DI TRAINING: il candidato "corretto" e' fissato a priori
# (regola esterna, il sistema non la conosce, deve scoprirla per
# tentativi)
# ---------------------------------------------------------------
CORRECT_CANDIDATE = 1  # candidato 2 (indice 1) e' quello "giusto"
LEARNING_RATE = 0.15
N_TRIALS = 60
STIM_CENTER = 0.3

weights = np.array([0.5, 0.5, 0.5])  # si parte da pesi UGUALI, nessuna preferenza
weight_history = [weights.copy()]
choices = []
rewards = []

print(f"Candidato corretto (nascosto al sistema): {CORRECT_CANDIDATE+1}")
print(f"Pesi iniziali: {weights}\n")

for trial in range(N_TRIALS):
    winner = run_trial(STIM_CENTER, weights, seed_noise=trial+1)
    if winner is None:
        choices.append(-1)
        rewards.append(0)
        continue
    reward = 1.0 if winner == CORRECT_CANDIDATE else 0.0
    weights[winner] += LEARNING_RATE * (reward - weights[winner])
    weights = np.clip(weights, 0.05, 2.0)  # evita pesi negativi o instabili

    choices.append(winner)
    rewards.append(reward)
    weight_history.append(weights.copy())

    if trial % 10 == 0 or trial == N_TRIALS-1:
        print(f"Prova {trial+1:3d}: scelto candidato {winner+1}, "
              f"ricompensa={reward:.0f}, pesi ora={np.round(weights,3)}")

weight_history = np.array(weight_history)
choices = np.array(choices)

# accuratezza in finestre di 10 prove
window = 10
accuracy_curve = [np.mean(np.array(choices[max(0,k-window):k+1]) == CORRECT_CANDIDATE)
                   for k in range(len(choices))]

print(f"\nAccuratezza nelle prime {window} prove: "
      f"{np.mean(choices[:window]==CORRECT_CANDIDATE)*100:.0f}%")
print(f"Accuratezza nelle ultime {window} prove: "
      f"{np.mean(choices[-window:]==CORRECT_CANDIDATE)*100:.0f}%")
print(f"Pesi finali: {np.round(weights,3)} "
      f"(il candidato corretto, {CORRECT_CANDIDATE+1}, dovrebbe avere il peso piu alto)")

# ---------------------------------------------------------------
# VISUALIZZAZIONE
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for i in range(3):
    axes[0].plot(weight_history[:, i], label=f'Peso candidato {i+1}'
                 + (' (corretto)' if i == CORRECT_CANDIDATE else ''),
                 linewidth=2 if i == CORRECT_CANDIDATE else 1.2)
axes[0].set_xlabel('Numero di prova')
axes[0].set_ylabel('Peso associativo')
axes[0].set_title('Evoluzione dei pesi durante il training')
axes[0].legend(loc='upper left', fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].plot(accuracy_curve, color='#2ecc71', linewidth=1.5)
axes[1].axhline(1/3, color='gray', linestyle='--', linewidth=1,
                 label='livello casuale (1/3)')
axes[1].set_xlabel('Numero di prova')
axes[1].set_ylabel(f'Accuratezza (media mobile su {window} prove)')
axes[1].set_title('Curva di apprendimento: la scelta corretta diventa piu\' frequente')
axes[1].legend(loc='lower right', fontsize=8)
axes[1].grid(alpha=0.3)
axes[1].set_ylim(-0.05, 1.05)

plt.tight_layout()
plt.savefig('rl_bandit.png', dpi=150)
print("\nGrafico salvato in rl_bandit.png")
