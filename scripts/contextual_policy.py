"""
Estensione del bandito a 3 braccia (rl_bandit.py) a 3 contesti indipendenti:
la posizione di stimolo su A (0.1/0.5/0.8) seleziona quale delle 3 risposte
e' corretta. Tabella pesi[contesto][candidato] invece di un vettore unico --
verifica a fine training che il sistema distingua i contesti, non solo
preferisca un candidato in generale.
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
    dists = [circ_dist(winner_pos - ((stim_center_A+off) % 1.0)) for off in CANDIDATE_OFFSETS]
    return int(np.argmin(dists))


CONTEXTS = [0.1, 0.5, 0.8]
CORRECT_PER_CONTEXT = [0, 1, 2]
LEARNING_RATE = 0.15
N_TRIALS = 150   # piu' prove del bandito singolo: 3 associazioni da imparare

weights_table = np.full((3, 3), 0.5)
history = np.zeros((N_TRIALS, 3, 3))
choices_per_context = {0: [], 1: [], 2: []}

rng_context = np.random.RandomState(0)

print("Contesti e risposte corrette (nascoste al sistema):")
for c, correct in enumerate(CORRECT_PER_CONTEXT):
    print(f"  Contesto {c+1} (stimolo a x={CONTEXTS[c]}) -> risposta corretta: candidato {correct+1}")
print()

for trial in range(N_TRIALS):
    context = rng_context.randint(0, 3)
    stim_pos = CONTEXTS[context]
    correct = CORRECT_PER_CONTEXT[context]

    winner = run_trial(stim_pos, weights_table[context], seed_noise=trial+1)
    if winner is not None:
        reward = 1.0 if winner == correct else 0.0
        weights_table[context, winner] += LEARNING_RATE * (reward - weights_table[context, winner])
        weights_table = np.clip(weights_table, 0.05, 2.0)
        choices_per_context[context].append(1 if winner == correct else 0)

    history[trial] = weights_table.copy()

    if trial % 30 == 0 or trial == N_TRIALS-1:
        print(f"Prova {trial+1:3d} (contesto {context+1}): "
              f"pesi contesto {context+1} ora = {np.round(weights_table[context],3)}")

print("\n" + "=" * 65)
print("VALUTAZIONE FINALE (pesi appresi, test pulito senza aggiornamento)")
print("=" * 65)
for c in range(3):
    test_wins = []
    for t in range(10):
        w = run_trial(CONTEXTS[c], weights_table[c], seed_noise=9000+c*10+t)
        test_wins.append(w)
    acc = np.mean(np.array(test_wins) == CORRECT_PER_CONTEXT[c])
    print(f"  Contesto {c+1} (x={CONTEXTS[c]}): accuratezza sul candidato corretto "
          f"({CORRECT_PER_CONTEXT[c]+1}) = {acc*100:.0f}%  |  pesi finali = {np.round(weights_table[c],3)}")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

window = 15
colors = ['#4a5aff', '#ff6b35', '#2ecc71']
for c in range(3):
    seq = choices_per_context[c]
    curve = [np.mean(seq[max(0, k-window):k+1]) for k in range(len(seq))]
    axes[0].plot(curve, color=colors[c], linewidth=1.5,
                 label=f'Contesto {c+1} (corretto: candidato {CORRECT_PER_CONTEXT[c]+1})')
axes[0].axhline(1/3, color='gray', linestyle='--', linewidth=1, label='livello casuale (1/3)')
axes[0].set_xlabel(f'Numero di volte che il contesto e\' stato incontrato')
axes[0].set_ylabel(f'Accuratezza (media mobile su {window} prove)')
axes[0].set_title('Curve di apprendimento, una per contesto')
axes[0].legend(loc='lower right', fontsize=7)
axes[0].grid(alpha=0.3)
axes[0].set_ylim(-0.05, 1.05)

im = axes[1].imshow(weights_table, cmap='viridis', vmin=0, vmax=1.2, aspect='auto')
axes[1].set_xticks(range(3)); axes[1].set_xticklabels(['Candidato 1', 'Candidato 2', 'Candidato 3'])
axes[1].set_yticks(range(3)); axes[1].set_yticklabels([f'Contesto {c+1}' for c in range(3)])
axes[1].set_title('Politica appresa: pesi[contesto][candidato]\n(diagonale = risposte corrette)')
for c in range(3):
    for cand in range(3):
        axes[1].text(cand, c, f'{weights_table[c,cand]:.2f}', ha='center', va='center',
                      color='white' if weights_table[c,cand] < 0.8 else 'black', fontsize=10)
plt.colorbar(im, ax=axes[1], label='peso appreso')

plt.tight_layout()
plt.savefig('contextual_policy.png', dpi=150)
print("\nGrafico salvato in contextual_policy.png")

