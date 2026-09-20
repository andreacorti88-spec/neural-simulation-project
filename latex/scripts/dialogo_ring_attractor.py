"""
Estensione di chained_attractors.py: la' A->B era a senso unico (mappa a
offset fissa, B non ha mai influenza su A). Qui la connessione e'
bidirezionale con DUE mappe di offset distinte (A->B e B->A), e i due ring
si scambiano il turno esplicitamente invece di co-evolvere in continuo --
altrimenti la dinamica converge quasi subito a un unico stato fuso e non
c'e' nulla da osservare come "scambio".

Un'informazione iniziale (un valore scalare in [0,1), lo stesso formalismo
gia' usato per lo stimolo su A negli script precedenti) viene data solo ad
A al turno 0. Da li' in poi nessun altro ingresso esterno: ogni turno
successivo il ring che "ascolta" riceve solo la proiezione della posizione
corrente dell'altro ring. Si registra la posizione del bump a ogni turno --
la "trascrizione" dello scambio -- e si verifica empiricamente cosa succede
davvero (convergenza a un punto fisso, ciclo stabile, o deriva), senza
assumerlo a priori.
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

A_exc = 20.0
sigma_exc = 0.05
global_inhib = 5.0
r_max = 5.0
stim_strength = 3.0
stim_width = 0.03
stim_duration_ms = 80

d_matrix = circ_dist(x[:, None] - x[None, :])
J_exc = A_exc * np.exp(-d_matrix**2 / (2*sigma_exc**2))

OFFSET_AB = 0.3   # come B interpreta cio' che sente da A
OFFSET_BA = 0.2   # come A interpreta cio' che sente da B (diversa da AB apposta:
                   # con due mappe uguali lo scambio tornerebbe sempre allo stesso punto)
COUPLING_STRENGTH = 1.3

T_TURN = 200        # ms per turno
N_TURNS = 6          # turno 0 = messaggio iniziale su A, poi 5 turni di ascolto alternato
steps_per_turn = int(T_TURN/dt)


def F(u):
    return np.clip(u, 0, r_max)


def bump_pos(r):
    return x[np.argmax(r)] if r.max() > 0.5 else None


def run_dialogo(messaggio, bidirezionale=True):
    """Turno 0: messaggio esterno su A. Turni dispari: B ascolta A.
    Turni pari (>0): A ascolta B (solo se bidirezionale=True)."""
    rA = np.zeros(N)
    rB = np.zeros(N)
    stim_profile = stim_strength * np.exp(-circ_dist(x-messaggio)**2/(2*stim_width**2))
    stim_duration_steps = int(stim_duration_ms/dt)

    history_A = np.zeros((N_TURNS*steps_per_turn, N))
    history_B = np.zeros((N_TURNS*steps_per_turn, N))
    pos_A_turni = []
    pos_B_turni = []

    for turn in range(N_TURNS):
        listener = 'B' if turn % 2 == 1 else ('A' if turn > 0 else None)

        for step_in_turn in range(steps_per_turn):
            step = turn*steps_per_turn + step_in_turn

            I_ext_A = stim_profile if (turn == 0 and step_in_turn < stim_duration_steps) else 0.0

            cross_to_A = np.zeros(N)
            if listener == 'A' and bidirezionale:
                posB = bump_pos(rB)
                if posB is not None:
                    target = (posB + OFFSET_BA) % 1.0
                    cross_to_A = COUPLING_STRENGTH * rB.max() * np.exp(
                        -circ_dist(x-target)**2/(2*stim_width**2))

            cross_to_B = np.zeros(N)
            if listener == 'B':
                posA = bump_pos(rA)
                if posA is not None:
                    target = (posA + OFFSET_AB) % 1.0
                    cross_to_B = COUPLING_STRENGTH * rA.max() * np.exp(
                        -circ_dist(x-target)**2/(2*stim_width**2))

            rec_A = J_exc @ rA / N - global_inhib*rA.mean()
            rA = rA + dt*(-rA + F(rec_A + I_ext_A + cross_to_A)) / tau

            rec_B = J_exc @ rB / N - global_inhib*rB.mean()
            rB = rB + dt*(-rB + F(rec_B + cross_to_B)) / tau

            history_A[step] = rA
            history_B[step] = rB

        pos_A_turni.append(bump_pos(rA))
        pos_B_turni.append(bump_pos(rB))

    return pos_A_turni, pos_B_turni, history_A, history_B


if __name__ == '__main__':
    print("=" * 70)
    print("TEST 1: scambio bidirezionale per diversi messaggi iniziali su A")
    print("=" * 70)
    messaggi_test = [0.10, 0.40, 0.70]
    for msg in messaggi_test:
        pos_A, pos_B, _, _ = run_dialogo(msg)
        trascrizione = []
        for t in range(N_TURNS):
            chi = 'A (messaggio)' if t == 0 else ('B ascolta' if t % 2 == 1 else 'A ascolta')
            pa = f"{pos_A[t]:.3f}" if pos_A[t] is not None else "nessun bump"
            pb = f"{pos_B[t]:.3f}" if pos_B[t] is not None else "nessun bump"
            trascrizione.append(f"    turno {t} ({chi}): A={pa}  B={pb}")
        print(f"\n  Messaggio iniziale x={msg:.2f}:")
        print("\n".join(trascrizione))

    print("\n" + "=" * 70)
    print("TEST 2: controllo -- senza canale B->A, A cambia mai dopo il turno 0?")
    print("=" * 70)
    pos_A_uni, pos_B_uni, _, _ = run_dialogo(0.40, bidirezionale=False)
    # turno 0->1: il bump si assesta dallo stimolo grezzo alla sua posizione di riposo
    # (il bias di assestamento del ring attractor e' un artefatto gia' noto e documentato
    # nel resoconto del progetto, sezione 9 -- non e' quello che questo test verifica).
    # Qui verifichiamo che DOPO l'assestamento (turni 1..5) A resti davvero fissa,
    # dato che senza canale B->A non ha piu' alcun ingresso.
    pos_A_dopo_assestamento = [p for p in pos_A_uni[1:] if p is not None]
    stabile = all(abs(circ_dist(p - pos_A_uni[1])) < 0.01 for p in pos_A_dopo_assestamento)
    print(f"  Posizioni di A nei turni successivi (solo A->B attivo): "
          f"{[f'{p:.3f}' if p is not None else None for p in pos_A_uni]}")
    print(f"  (turno 0->1: assestamento naturale del bump, atteso; "
          f"da qui in poi nessun ingresso esterno raggiunge A)")
    if stabile:
        print("  >>> Confermato: dopo l'assestamento iniziale, senza il canale di ritorno "
              "A resta fissa (come in chained_attractors.py). <<<")
    else:
        print("  >>> ATTENZIONE: A si muove anche dopo l'assestamento e senza canale di "
              "ritorno -- controllare il codice, non dovrebbe succedere. <<<")

    print("\n" + "=" * 70)
    print("TEST 3: lo scambio bidirezionale si stabilizza o continua a variare?")
    print("=" * 70)
    pos_A_bi, pos_B_bi, history_A, history_B = run_dialogo(0.40, bidirezionale=True)
    ultimo_giro = circ_dist(pos_A_bi[-1] - pos_A_bi[-3]) if pos_A_bi[-1] is not None and pos_A_bi[-3] is not None else None
    if ultimo_giro is not None:
        print(f"  Spostamento di A tra il penultimo e l'ultimo suo turno di ascolto: "
              f"{abs(ultimo_giro):.4f}")
        if abs(ultimo_giro) < 0.01:
            print("  >>> Lo scambio converge a un punto fisso (un 'accordo'). <<<")
        else:
            print("  >>> Lo scambio NON converge in 6 turni: la posizione continua "
                  "a variare a ogni giro (deriva o ciclo). <<<")

    # figura: spazio-tempo dei due ring + trascrizione posizione-per-turno
    t_axis = np.arange(N_TURNS*steps_per_turn)*dt
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=False)

    im0 = axes[0].imshow(history_A.T, aspect='auto', origin='lower', cmap='inferno',
                          extent=[0, N_TURNS*T_TURN, 0, 1])
    for t in range(1, N_TURNS):
        axes[0].axvline(t*T_TURN, color='cyan', linestyle='--', linewidth=0.8, alpha=0.6)
    axes[0].set_ylabel('Posizione su Ring A')
    axes[0].set_title(f'Ring A (messaggio iniziale x=0.40) -- linee tratteggiate = confine tra turni')
    plt.colorbar(im0, ax=axes[0], label='r_A(x,t)')

    im1 = axes[1].imshow(history_B.T, aspect='auto', origin='lower', cmap='viridis',
                          extent=[0, N_TURNS*T_TURN, 0, 1])
    for t in range(1, N_TURNS):
        axes[1].axvline(t*T_TURN, color='cyan', linestyle='--', linewidth=0.8, alpha=0.6)
    axes[1].set_xlabel('Tempo (ms)')
    axes[1].set_ylabel('Posizione su Ring B')
    axes[1].set_title('Ring B')
    plt.colorbar(im1, ax=axes[1], label='r_B(x,t)')

    turni = list(range(N_TURNS))
    pa_plot = [p if p is not None else np.nan for p in pos_A_bi]
    pb_plot = [p if p is not None else np.nan for p in pos_B_bi]
    axes[2].plot(turni, pa_plot, 'o-', color='C3', label='A (posizione a fine turno)')
    axes[2].plot(turni, pb_plot, 's-', color='C0', label='B (posizione a fine turno)')
    axes[2].set_xlabel('Turno')
    axes[2].set_ylabel('Posizione sul ring [0,1)')
    axes[2].set_title('Trascrizione dello scambio: posizione di A e B turno per turno')
    axes[2].legend(loc='best', fontsize=8)
    axes[2].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('dialogo_ring_attractor.png', dpi=150)
    print("\nGrafico salvato in dialogo_ring_attractor.png")
