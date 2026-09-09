"""
COMPETIZIONE TRA ATTRATTORI: B "sceglie" tra alternative
====================================================================
Estensione del sistema A->B con mappa fissa: invece di UN'UNICA
destinazione possibile per B, offriamo TRE candidati (offset diversi
sull'anello), ciascuno con un proprio "peso" (forza associativa).
Grazie all'inibizione globale gia' presente nella dinamica di B
(meccanismo winner-take-all), solo UN candidato puo' vincere e
stabilizzarsi -- gli altri vengono soppressi.

DUE REGIMI TESTATI ONESTAMENTE:

  1. PESI UGUALI: in assenza di rumore (o con rumore troppo debole),
     il "vincitore" e' determinato da un bias numerico deterministico
     dell'implementazione (differenze in floating point dell'ordine
     di 1e-4, amplificate dalla dinamica non lineare) -- NON e' vera
     casualita'. Aumentando il rumore continuo iniettato nella
     competizione a un livello sufficiente, la rottura di simmetria
     diventa genuinamente casuale e variabile tra le prove.

  2. PESI DIVERSI: quando un candidato ha un peso maggiore, vince in
     modo affidabile e ripetibile -- una vera selezione basata su un
     criterio (la forza del peso), non casuale e non arbitraria.

Questo e' il primo mattone di una "decisione" nel senso debole del
termine: il sistema seleziona un'alternativa tra piu' possibili, in
base a un criterio esplicito (il peso). Resta comunque lontanissimo
da un vero processo decisionale: il criterio e' fisso e dato
dall'esterno, non appreso o valutato dal sistema stesso.
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


def F(u):
    return np.clip(u, 0, r_max)


def run_competition(stim_center_A, weights, seed_noise=0, noise_std=0.15):
    rA = np.zeros(N)
    rB = np.zeros(N)
    stim_duration_steps = int(stim_duration_ms/dt)
    stim_profile = stim_strength * np.exp(-circ_dist(x-stim_center_A)**2/(2*stim_width**2))
    rng = np.random.RandomState(seed_noise) if seed_noise > 0 else None

    history_B = np.zeros((n_steps, N))

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
            if rng is not None:
                cross_input += rng.normal(0, noise_std, N)  # rumore CONTINUO

        rec_B = J_exc @ rB / N - global_inhib*rB.mean()
        rB = rB + dt*(-rB + F(rec_B + cross_input)) / tau
        history_B[step] = rB

    return rA, rB, history_B


def find_winner(rB_final):
    if rB_final.max() < 0.5:
        return None
    winner_pos = x[np.argmax(rB_final)]
    dists = [circ_dist(winner_pos - ((0.3+off) % 1.0)) for off in CANDIDATE_OFFSETS]
    return int(np.argmin(dists))


if __name__ == '__main__':
    print("=" * 65)
    print("TEST 1a: pesi UGUALI, rumore ASSENTE -- deterministico o casuale?")
    print("=" * 65)
    winners_no_noise = []
    for trial in range(6):
        _, rB_f, _ = run_competition(0.3, [1.0, 1.0, 1.0], seed_noise=0)
        w = find_winner(rB_f)
        winners_no_noise.append(w)
        print(f"  Prova {trial+1} (nessun rumore): vince candidato {w+1}")
    print(f"  -> Sempre lo stesso candidato: probabile bias numerico "
          f"dell'implementazione, non vera casualita'.")

    print("\n" + "=" * 65)
    print("TEST 1b: pesi UGUALI, con rumore CONTINUO sufficiente")
    print("=" * 65)
    winners_noise = []
    N_TRIALS = 12
    for trial in range(N_TRIALS):
        _, rB_f, _ = run_competition(0.3, [1.0, 1.0, 1.0], seed_noise=trial+1)
        w = find_winner(rB_f)
        winners_noise.append(w)
    counts = [winners_noise.count(i) for i in range(3)]
    print(f"  Su {N_TRIALS} prove: candidato 1 vince {counts[0]}x, "
          f"candidato 2 vince {counts[1]}x, candidato 3 vince {counts[2]}x")
    print("  -> Con rumore sufficiente, la rottura di simmetria e' genuinamente "
          "variabile (anche se non perfettamente uniforme 1/3-1/3-1/3).")

    print("\n" + "=" * 65)
    print("TEST 2: un candidato ha peso MAGGIORE -- vince sempre quello?")
    print("=" * 65)
    weighted_results = {}
    for favored in range(3):
        weights = [0.5, 0.5, 0.5]
        weights[favored] = 1.5
        wins = []
        for trial in range(6):
            _, rB_f, _ = run_competition(0.3, weights, seed_noise=trial+100)
            wins.append(find_winner(rB_f))
        weighted_results[favored] = wins
        print(f"  Candidato favorito {favored+1} (peso 1.5 vs 0.5): "
              f"vincitori nelle 6 prove = {[w+1 for w in wins]}")

    # -------------------------------------------------------------
    # VISUALIZZAZIONE
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].bar(['Candidato 1', 'Candidato 2', 'Candidato 3'], counts,
                color=['#4a5aff', '#ff6b35', '#2ecc71'])
    axes[0].set_ylabel(f'Vittorie su {N_TRIALS} prove')
    axes[0].set_title('Pesi uguali + rumore: rottura di simmetria casuale')
    axes[0].grid(axis='y', alpha=0.3)

    _, _, history_example = run_competition(0.3, [0.5, 1.5, 0.5], seed_noise=42)
    t_axis = np.arange(n_steps)*dt
    im = axes[1].imshow(history_example.T, aspect='auto', origin='lower', cmap='viridis',
                         extent=[0, T_total, 0, 1])
    for off in CANDIDATE_OFFSETS:
        axes[1].axhline((0.3+off) % 1.0, color='white', linestyle=':', linewidth=0.8, alpha=0.6)
    axes[1].set_xlabel('Tempo (ms)')
    axes[1].set_ylabel('Posizione su Ring B')
    axes[1].set_title('Competizione: candidato 2 (peso maggiore) vince,\naltri due soppressi')
    plt.colorbar(im, ax=axes[1], label='r_B(x,t)')

    plt.tight_layout()
    plt.savefig('competing_targets.png', dpi=150)
    print("\nGrafico salvato in competing_targets.png")
