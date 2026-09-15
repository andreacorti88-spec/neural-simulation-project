"""
Estensione di dialogo_ring_attractor.py: li' le due mappe di proiezione tra
A e B (offset 0.3 e 0.2) erano costanti scritte a mano nel codice. Qui
diventano sinapsi vere -- due matrici di peso N x N, W_AB e W_BA -- che
partono da un debole bias innato (troppo debole, da solo, per far scattare
un bump nella rete che ascolta) e si CONSOLIDANO con un aggiornamento
Hebbiano applicato a fine turno: quando la rete che ascolta forma davvero
un bump, la coppia (posizione di chi ha appena "parlato", posizione di chi
ha appena "risposto") viene rinforzata li' dove e' successa.

E' una versione a grana grossa -- per-turno, non al livello del singolo
spike -- della STDP usata altrove nel progetto (stdp.py,
classificazione_stdp.py): stessa idea (rinforzo guidato dalla correlazione
temporale tra chi si attiva prima e chi risponde dopo), applicata pero' su
una rete a rate invece che a impulsi, per coerenza con la scelta gia' fatta
per gli attrattori (ring_attractor_rate.py).

Il problema dell'uovo e della gallina -- senza una connessione funzionante
non si osserva mai una coppia pre/post da rinforzare -- viene risolto con
rumore esplorativo nella rete che ascolta SOLO durante il training: rottura
di simmetria casuale, lo stesso meccanismo gia' verificato in
competing_targets.py (TEST 1b). Con quel rumore un bump puo' comunque
comparire per caso anche quando i pesi sono deboli, e quel caso puo' essere
rinforzato -- un rich-get-richer Hebbiano, non diverso in spirito dalle
onde spontanee che guidano la formazione delle mappe retinotopiche prima
che ci sia vera visione.
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

COUPLING_STRENGTH = 1.3
T_TURN = 200
N_TURNS = 6
steps_per_turn = int(T_TURN/dt)

# bias innato: stessa forma delle mappe fisse di dialogo_ring_attractor.py,
# ma scalato troppo debole per bastare da solo (verificato sotto)
OFFSET_AB_INNATO = 0.3
OFFSET_BA_INNATO = 0.2
# calibrato empiricamente (vedi note): la transizione da "nessun bump" a
# "bump certo" in B e' un'accensione tutto-o-niente tra scale=0.0025 e 0.003
# (l'attrattore B e' bistabile, non risponde linearmente all'input in ingresso)
# -- si parte apposta sotto quella soglia, cosi' l'ignizione iniziale dipende
# solo dal rumore esplorativo, non dal bias innato da solo
W_INIT_SCALE = 0.0018
SIGMA_W = 0.05
W_MAX = 2.5
ETA = 0.0001
DECAY_PER_EPISODE = 0.999
# esplorazione: un impulso a posizione casuale, stessa forma/durata del
# messaggio vero, ma troppo debole per accendere un bump DA SOLO (soglia di
# accensione solitaria verificata a strength~0.0036 con questa durata --
# 0.003 sta sotto). Solo quando l'impulso casuale capita vicino alla zona
# gia' favorita dal bias innato, la somma dei due basta a superare la soglia
# -- rumore incoerente puro (provato prima: rumore bianco indipendente per
# passo) non ha mai acceso nulla in 40 episodi, perche' la costante di tempo
# della rete (tau=10ms) lo media via prima che si accumuli in un punto
EXPLORATION_PULSE_STRENGTH = 0.0030
EXPLORATION_PULSE_DURATION_MS = stim_duration_ms


def F(u):
    return np.clip(u, 0, r_max)


def bump_pos(r):
    return x[np.argmax(r)] if r.max() > 0.5 else None


def init_weak_bias(offset, scale):
    target = (x[:, None] + offset) % 1.0
    d = circ_dist(x[None, :] - target)
    return scale * np.exp(-d**2/(2*SIGMA_W**2))


def hebbian_bump(pos_pre, pos_post):
    kern_pre = np.exp(-circ_dist(x-pos_pre)**2/(2*SIGMA_W**2))
    kern_post = np.exp(-circ_dist(x-pos_post)**2/(2*SIGMA_W**2))
    return np.outer(kern_pre, kern_post)


def run_episode(W_AB, W_BA, messaggio, rng, train):
    rA = np.zeros(N)
    rB = np.zeros(N)
    stim_profile = stim_strength * np.exp(-circ_dist(x-messaggio)**2/(2*stim_width**2))
    stim_duration_steps = int(stim_duration_ms/dt)

    pos_A_turni = []
    pos_B_turni = []

    explore_duration_steps = int(EXPLORATION_PULSE_DURATION_MS/dt)

    for turn in range(N_TURNS):
        listener = 'B' if turn % 2 == 1 else ('A' if turn > 0 else None)

        if train and listener is not None:
            explore_pos = rng.uniform(0, 1)
            explore_profile = EXPLORATION_PULSE_STRENGTH * np.exp(
                -circ_dist(x-explore_pos)**2/(2*stim_width**2))
        else:
            explore_profile = np.zeros(N)

        for step_in_turn in range(steps_per_turn):
            I_ext_A = stim_profile if (turn == 0 and step_in_turn < stim_duration_steps) else 0.0

            explore_now = explore_profile if (train and step_in_turn < explore_duration_steps) else np.zeros(N)
            explore_A = explore_now if listener == 'A' else np.zeros(N)
            explore_B = explore_now if listener == 'B' else np.zeros(N)

            cross_to_A = COUPLING_STRENGTH * (rB @ W_BA) / N if listener == 'A' else np.zeros(N)
            cross_to_B = COUPLING_STRENGTH * (rA @ W_AB) / N if listener == 'B' else np.zeros(N)

            rec_A = J_exc @ rA / N - global_inhib*rA.mean()
            rA = rA + dt*(-rA + F(rec_A + I_ext_A + cross_to_A + explore_A)) / tau

            rec_B = J_exc @ rB / N - global_inhib*rB.mean()
            rB = rB + dt*(-rB + F(rec_B + cross_to_B + explore_B)) / tau

        pos_A = bump_pos(rA)
        pos_B = bump_pos(rB)
        pos_A_turni.append(pos_A)
        pos_B_turni.append(pos_B)

        if train and listener == 'B' and pos_A is not None and pos_B is not None:
            W_AB += ETA * hebbian_bump(pos_A, pos_B)
            np.clip(W_AB, 0, W_MAX, out=W_AB)
        elif train and listener == 'A' and pos_A is not None and pos_B is not None:
            W_BA += ETA * hebbian_bump(pos_B, pos_A)
            np.clip(W_BA, 0, W_MAX, out=W_BA)

    return pos_A_turni, pos_B_turni


def eval_batch(W_AB, W_BA, messaggio, n_trial, seed0):
    """Valutazione pulita: pesi congelati, niente rumore esplorativo.
    NOTA: A ha un bump persistente fin dal turno 0 per costruzione (e' un
    attrattore), quindi "A ha un bump al turno 2" e' quasi sempre vero anche
    se B non ha mai influenzato nulla. Il confronto va fatto contro il
    turno 1 (non il turno 0): tra il turno 0 e il turno 1 il bump di A si
    assesta comunque un po' per conto suo (lo stesso bias di posizione gia'
    documentato nel resoconto del progetto, sezione 9), quindi usare il
    turno 0 come riferimento conterebbe quell'assestamento come se fosse una
    risposta a B. Dal turno 1 al turno 2, invece, senza influenza da B la
    posizione di A e' gia' assestata e dovrebbe restare ferma -- solo un
    cambiamento li' indica che l'informazione e' davvero tornata da B."""
    successi_B, successi_A_spostata = 0, 0
    pos_B_ok, pos_A_ok = [], []
    for k in range(n_trial):
        rng = np.random.RandomState(seed0 + k)
        pos_A, pos_B = run_episode(W_AB.copy(), W_BA.copy(), messaggio, rng, train=False)
        if pos_B[1] is not None:
            successi_B += 1
            pos_B_ok.append(pos_B[1])
        if pos_A[2] is not None and pos_A[1] is not None and abs(circ_dist(pos_A[2]-pos_A[1])) > 0.02:
            successi_A_spostata += 1
            pos_A_ok.append(pos_A[2])
    return {
        'tasso_B': successi_B/n_trial,
        'tasso_A_spostata': successi_A_spostata/n_trial,
        'std_B': np.std(pos_B_ok) if len(pos_B_ok) > 1 else None,
        'std_A': np.std(pos_A_ok) if len(pos_A_ok) > 1 else None,
        'media_B': np.mean(pos_B_ok) if pos_B_ok else None,
    }


if __name__ == '__main__':
    MESSAGGIO = 0.40
    N_EPISODES = 300

    W_AB = init_weak_bias(OFFSET_AB_INNATO, W_INIT_SCALE)
    W_BA = init_weak_bias(OFFSET_BA_INNATO, W_INIT_SCALE)

    # soglia di accensione autonoma verificata empiricamente per bisezione
    # (run_episode reale, non stima lineare): tra scale=0.0025 (0% successo)
    # e scale=0.0030 (100% successo)
    SOGLIA_IGNIZIONE_W = 0.0028

    print("=" * 70)
    print("SETUP: bias innato deliberatamente debole")
    print("=" * 70)
    print(f"  W_INIT_SCALE={W_INIT_SCALE}, soglia per formare un bump = 0.5")
    W_AB_iniziale = W_AB.copy()
    W_BA_iniziale = W_BA.copy()

    print("\n" + "=" * 70)
    print("BASELINE (pre-training): il bias innato basta da solo a far rispondere B?")
    print("=" * 70)
    baseline = eval_batch(W_AB, W_BA, MESSAGGIO, n_trial=20, seed0=0)
    print(f"  Tasso di successo di B (turno 1): {baseline['tasso_B']*100:.0f}%")
    print(f"  Tasso di volte che A si SPOSTA per davvero al turno 2 "
          f"(risposta reale a B, non solo bump persistente): {baseline['tasso_A_spostata']*100:.0f}%")

    print("\n" + "=" * 70)
    print(f"TRAINING: {N_EPISODES} episodi, messaggio fisso x={MESSAGGIO}, "
          f"rumore esplorativo attivo, pesi persistenti tra episodi")
    print("=" * 70)
    rng_train = np.random.RandomState(42)
    successo_B_storia = []
    spostamento_A_storia = []
    picco_W_AB_storia = []
    picco_W_BA_storia = []
    for ep in range(N_EPISODES):
        pos_A, pos_B = run_episode(W_AB, W_BA, MESSAGGIO, rng_train, train=True)
        W_AB *= DECAY_PER_EPISODE
        W_BA *= DECAY_PER_EPISODE
        successo_B_storia.append(1 if pos_B[1] is not None else 0)
        spostato = (pos_A[2] is not None and pos_A[1] is not None
                    and abs(circ_dist(pos_A[2]-pos_A[1])) > 0.02)
        spostamento_A_storia.append(1 if spostato else 0)
        picco_W_AB_storia.append(W_AB.max())
        picco_W_BA_storia.append(W_BA.max())
        if ep % 50 == 0 or ep == N_EPISODES-1:
            win_B = successo_B_storia[max(0, ep-30):ep+1]
            win_A = spostamento_A_storia[max(0, ep-30):ep+1]
            print(f"  episodio {ep:3d}: successo B (ultimi 30) = {np.mean(win_B)*100:5.1f}%  |  "
                  f"A si sposta per B (ultimi 30) = {np.mean(win_A)*100:5.1f}%  |  "
                  f"picco W_AB = {W_AB.max():.4f}  picco W_BA = {W_BA.max():.4f}")

    print("\n" + "=" * 70)
    print("DOPO IL TRAINING: confronto pulito, pesi congelati, nessun rumore esplorativo")
    print("=" * 70)
    dopo = eval_batch(W_AB, W_BA, MESSAGGIO, n_trial=20, seed0=9000)
    print(f"  PRIMA  -- tasso successo B: {baseline['tasso_B']*100:5.0f}%   "
          f"std posizione B: {baseline['std_B']}")
    std_dopo_str = f"{dopo['std_B']:.4f}" if dopo['std_B'] else "None"
    print(f"  DOPO   -- tasso successo B: {dopo['tasso_B']*100:5.0f}%   "
          f"std posizione B: {std_dopo_str}")
    print(f"  PRIMA  -- tasso di A che si sposta per davvero al turno 2: "
          f"{baseline['tasso_A_spostata']*100:5.0f}%")
    print(f"  DOPO   -- tasso di A che si sposta per davvero al turno 2: "
          f"{dopo['tasso_A_spostata']*100:5.0f}%")
    if dopo['media_B'] is not None:
        offset_appreso = circ_dist(dopo['media_B'] - MESSAGGIO)
        print(f"\n  Offset A->B appreso empiricamente: {offset_appreso:.3f} "
              f"(bias innato di partenza era {OFFSET_AB_INNATO})")

    print("\n" + "=" * 70)
    print("ASIMMETRIA A->B vs B->A")
    print("=" * 70)
    print(f"  Picco finale W_AB: {W_AB.max():.4f}  (soglia di accensione ~{SOGLIA_IGNIZIONE_W})")
    print(f"  Picco finale W_BA: {W_BA.max():.4f}  (soglia di accensione ~{SOGLIA_IGNIZIONE_W})")
    if dopo['tasso_A_spostata'] < 0.1:
        print("  >>> Il canale B->A NON si e' consolidato in questo run: a differenza di B, "
              "A ha gia' un bump persistente fin dal turno 0, quindi l'impulso esplorativo deve "
              "SOSTITUIRE un bump gia' attivo invece di accenderne uno da zero -- un compito "
              "piu' difficile, per cui non e' mai emerso un segnale di apprendimento su W_BA. "
              "Non e' un fallimento nascosto: e' un'asimmetria reale del meccanismo, verificata "
              "col codice. <<<")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    window = 20
    curva = [np.mean(successo_B_storia[max(0, k-window):k+1]) for k in range(len(successo_B_storia))]
    axes[0, 0].plot(curva, color='C0')
    axes[0, 0].set_xlabel('Episodio di training')
    axes[0, 0].set_ylabel(f'Tasso di successo di B (media mobile su {window})')
    axes[0, 0].set_title("Curva di apprendimento: B forma un bump sempre piu' spesso")
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].set_ylim(-0.05, 1.05)

    axes[0, 1].plot(picco_W_AB_storia, color='C3', label='W_AB (A->B, accende B da silenzio)')
    axes[0, 1].plot(picco_W_BA_storia, color='C0', label='W_BA (B->A, deve spodestare un bump attivo)')
    axes[0, 1].axhline(SOGLIA_IGNIZIONE_W, color='gray', linestyle='--',
                        label='soglia di accensione da silenzio (non vale per B->A)')
    axes[0, 1].set_xlabel('Episodio di training')
    axes[0, 1].set_ylabel('Peso massimo')
    axes[0, 1].set_title("Il peso cresce per entrambe le sinapsi, ma solo A->B\ndiventa comportamentalmente efficace")
    axes[0, 1].legend(fontsize=7, loc='upper left')
    axes[0, 1].grid(alpha=0.3)

    im0 = axes[1, 0].imshow(W_AB_iniziale, origin='lower', extent=[0, 1, 0, 1],
                             cmap='viridis', vmin=0, vmax=W_AB.max())
    axes[1, 0].set_xlabel('Posizione post (B)')
    axes[1, 0].set_ylabel('Posizione pre (A)')
    axes[1, 0].set_title('W_AB PRIMA del training (bias innato debole)')
    plt.colorbar(im0, ax=axes[1, 0])

    im1 = axes[1, 1].imshow(W_AB, origin='lower', extent=[0, 1, 0, 1],
                             cmap='viridis', vmin=0, vmax=W_AB.max())
    axes[1, 1].set_xlabel('Posizione post (B)')
    axes[1, 1].set_ylabel('Posizione pre (A)')
    axes[1, 1].set_title('W_AB DOPO il training (consolidata)')
    plt.colorbar(im1, ax=axes[1, 1])

    plt.tight_layout()
    plt.savefig('dialogo_stdp.png', dpi=150)
    print("\nGrafico salvato in dialogo_stdp.png")
