"""
Catena di due ring attractor: A riceve lo stimolo diretto, B riceve solo
la proiezione della posizione del bump di A (offset fisso 0.3 sull'anello,
mappa associativa hard-coded). Nessun meccanismo di scelta -- propagazione
causale attraverso uno stato intermedio stabile, propedeutico a
competizione/rinforzo tra piu' mappe candidate (vedi script successivi).
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

OFFSET_AB = 0.3
COUPLING_STRENGTH = 1.3


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

        # B: nessuno stimolo diretto, solo la proiezione spostata di A
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
    print("TEST 1: A stimolato in posizioni diverse, B guidato solo dall'influenza di A")
    print("=" * 65)
    test_centers = [0.1, 0.3, 0.5, 0.7]
    for center_A in test_centers:
        rA_final, rB_final, _, _ = run_chained(center_A)
        posA = x[np.argmax(rA_final)]
        posB = x[np.argmax(rB_final)] if rB_final.max() > 0.5 else None
        posB_str = f"{posB:.3f}" if posB is not None else "NESSUN BUMP"
        expected_B = (center_A + OFFSET_AB) % 1.0
        print(f"  Stimolo A a x={center_A:.2f} -> A: x={posA:.3f}, "
              f"B: x={posB_str} (atteso ~{expected_B:.3f})")

    print("\n" + "=" * 65)
    print("TEST 2: controllo, senza stimolo ad A ne' A ne' B dovrebbero attivarsi")
    print("=" * 65)
    rA0, rB0, _, _ = run_chained(0.3, apply_stim=False)
    print(f"  A picco finale: {rA0.max():.6f}")
    print(f"  B picco finale: {rB0.max():.6f}")
    if rA0.max() < 0.01 and rB0.max() < 0.01:
        print("  >>> Confermato: nessuna attivazione spontanea. <<<")

    # spazio-tempo di A e B, caso stimolo a x=0.3
    _, _, history_A, history_B = run_chained(0.3)
    t_axis = np.arange(n_steps)*dt

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    im0 = axes[0].imshow(history_A.T, aspect='auto', origin='lower', cmap='inferno',
                          extent=[0, T_total, 0, 1])
    axes[0].axvline(stim_duration_ms, color='cyan', linestyle='--', linewidth=1.2,
                     label='fine stimolo esterno (solo su A)')
    axes[0].set_ylabel('Posizione su Ring A')
    axes[0].set_title('Ring A: stimolo esterno diretto')
    axes[0].legend(loc='upper right', fontsize=8)
    plt.colorbar(im0, ax=axes[0], label='r_A(x,t)')

    im1 = axes[1].imshow(history_B.T, aspect='auto', origin='lower', cmap='viridis',
                          extent=[0, T_total, 0, 1])
    axes[1].axvline(stim_duration_ms, color='cyan', linestyle='--', linewidth=1.2)
    axes[1].set_xlabel('Tempo (ms)')
    axes[1].set_ylabel('Posizione su Ring B')
    axes[1].set_title('Ring B: nessuno stimolo diretto, si attiva via proiezione da A '
                       '(nota il ritardo)')
    plt.colorbar(im1, ax=axes[1], label='r_B(x,t)')

    plt.tight_layout()
    plt.savefig('chained_attractors.png', dpi=150)
    print("\nGrafico salvato in chained_attractors.png")
