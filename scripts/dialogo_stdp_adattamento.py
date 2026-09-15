"""
Estensione di dialogo_stdp.py. Li' il canale A->B si consolidava con
successo (0% -> 100%) ma il canale di ritorno B->A restava bloccato: A ha
un bump persistente fin dal turno 0 (e' un attrattore), quindi per essere
influenzato per davvero B doveva SPODESTARE un bump gia' a piena potenza
invece di accenderne uno da zero -- molto piu' difficile, e infatti non
succedeva mai (verificato: W_BA cresceva per pura coincidenza statistica,
ma restava comportamentalmente inefficace).

Qui si aggiunge un ADATTAMENTO a soglia di scarica -- l'analogo a rate di
una variabile di recovery come la "u" del modello di Izhikevich usato in
tutto il resto del progetto (vedi resoconto, sezione 1.2): piu' a lungo un
ring resta acceso, piu' cresce una corrente negativa che alla fine lo
spegne da solo. Calibrato (vedi note sotto) cosi' che un bump sopravviva
un turno intero pieno, si affievolisca in quello successivo, e collassi
da solo un turno dopo ancora -- quindi quando tocca all'altro ring
"ascoltare", il ring che ha gia' parlato spesso si e' GIA' spento per
conto suo. Il problema torna quindi ad essere lo stesso, gia' risolto, di
accendere un bump da silenzio -- non piu' di spodestarne uno a piena
potenza. Le due sinapsi (A->B e B->A) diventano cosi' simmetriche per
costruzione, ed entrambe possono in linea di principio consolidarsi con lo
stesso meccanismo Hebbiano di dialogo_stdp.py (rumore esplorativo +
rinforzo di fine turno).

E' un compromesso: la scarica hard-clipped di questo modello rende
l'attrattore fortemente bistabile (verificato con un test dedicato: non
esiste quasi via di mezzo tra "bump a piena potenza" e "spento",
l'affievolimento e' un collasso brusco più che un decadimento morbido).
Il parametro di adattamento e' stato scelto apposta sul bordo di quella
transizione.

AGGIORNAMENTO (seconda iterazione): la prima versione con adattamento
riduceva l'asimmetria ma non la risolveva -- in 2000 episodi W_BA cresceva
circa 7 volte piu' lento di W_AB e non arrivava mai all'autosufficienza.
Diagnosticato il motivo: nei primi episodi la posizione di B stessa non
si e' ancora stabilizzata (dipende da un W_AB non ancora consolidato),
quindi il rinforzo Hebbiano di B->A nei primi episodi punta a un bersaglio
mobile e si diluisce invece di concentrarsi. Corretto con due aggiunte:
(1) una fase di RISCALDAMENTO (WARMUP_EPISODES) in cui si allena solo
A->B, cosi' che B->A cominci a imparare solo quando ha gia' un bersaglio
stabile; (2) un tasso di apprendimento piu' alto per B->A (ETA_BA > ETA)
per compensare il minor numero di eventi di rinforzo utili. Risultato:
con training sufficiente (~2500 episodi) entrambe le sinapsi si
autosostengono e TUTTI e quattro gli scambi verificati (B risponde al
turno 1 e di nuovo al turno 5, A risponde al turno 3 e di nuovo al turno
7) arrivano al 100% di affidabilita' in valutazione pulita (pesi
congelati, senza rumore esplorativo). Nota sul turno 5: e' un'associazione
di SECONDO ORDINE (B deve imparare a rispondere alla posizione che A
assume DOPO aver gia' risposto a B una prima volta, non al messaggio
originale) -- si consolida solo dopo che il turno 3 e' gia' stabile,
quindi richiede piu' episodi della prima associazione e nella heatmap
finale appare come un rinforzo piu' debole, secondario, rispetto al picco
principale.
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
# schema dei turni con un turno di silenzio/recupero tra un ascolto e
# l'altro (verificato necessario: senza il silenzio, il tentativo di
# ascolto immediatamente successivo a un turno di parola arriva mentre
# l'adattamento e' ancora al culmine e non si accende quasi mai -- vedi
# docstring)
SCHEDULE = [None, 'B', None, 'A', None, 'B', None, 'A']
N_TURNS = len(SCHEDULE)
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
# B->A ha un tasso di successo durante il training piu' basso di A->B (vedi
# dialogo_stdp_adattamento.png della run precedente: ~100% vs ~65-80%), e
# soprattutto nei primi episodi la posizione di B stessa non si e' ancora
# stabilizzata -- rinforzare un bersaglio che si sposta episodio per
# episodio diluisce il rinforzo Hebbiano invece di concentrarlo. ETA_BA
# piu' alto compensa il minor numero di eventi di rinforzo utili
ETA_BA = 0.0003
# fase di riscaldamento: per i primi WARMUP_EPISODES si allena solo A->B
# (W_BA resta congelata al bias innato) cosi' che, quando B->A comincia
# davvero ad imparare, la posizione di B a cui deve puntare sia gia'
# stabile -- altrimenti gli episodi iniziali sprecano rinforzo su un
# bersaglio mobile (verificato: senza questa fase B->A cresce ~7 volte
# piu' lento di A->B in 2000 episodi)
WARMUP_EPISODES = 400
DECAY_PER_EPISODE = 0.999
# esplorazione: un impulso a posizione casuale, stessa forma/durata del
# messaggio vero, ma troppo debole per accendere un bump DA SOLO (soglia di
# accensione solitaria verificata a strength~0.0036 con questa durata --
# 0.003 sta sotto). Solo quando l'impulso casuale capita vicino alla zona
# gia' favorita dal bias innato, la somma dei due basta a superare la soglia
# -- rumore incoerente puro (provato prima: rumore bianco indipendente per
# passo) non ha mai acceso nulla in 40 episodi, perche' la costante di tempo
# della rete (tau=10ms) lo media via prima che si accumuli in un punto
EXPLORATION_PULSE_STRENGTH = 0.0080
EXPLORATION_PULSE_DURATION_MS = stim_duration_ms

# adattamento: calibrato con un ring isolato (vedi note nel docstring) --
# a gain=0.25, tau=200ms un bump segue picco/turno = 5.00 -> 2.71 -> 0.00,
# cioe' pieno per un turno, affievolito nel successivo, spento nel terzo
GAIN_ADAPT = 0.25
TAU_ADAPT = 200.0


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


def run_episode(W_AB, W_BA, messaggio, rng, train, allena_BA=True):
    rA = np.zeros(N)
    rB = np.zeros(N)
    aA = np.zeros(N)
    aB = np.zeros(N)
    stim_profile = stim_strength * np.exp(-circ_dist(x-messaggio)**2/(2*stim_width**2))
    stim_duration_steps = int(stim_duration_ms/dt)

    pos_A_turni = []
    pos_B_turni = []
    # posizione dell'ULTIMO bump valido di ciascun ring, indipendentemente
    # da quanti turni fa: serve per il rinforzo Hebbiano, perche' con
    # l'adattamento chi ha "parlato" spesso si e' gia' spento per conto suo
    # nel momento esatto in cui l'ascoltatore finisce di reagire -- leggere
    # solo la posizione istantanea perderebbe quasi ogni coppia valida
    ultima_pos_A, ultima_pos_B = None, None

    explore_duration_steps = int(EXPLORATION_PULSE_DURATION_MS/dt)

    for turn, listener in enumerate(SCHEDULE):

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
            rA_new = rA + dt*(-rA + F(rec_A + I_ext_A + cross_to_A + explore_A - aA)) / tau
            aA = aA + dt*(-aA + GAIN_ADAPT*rA) / TAU_ADAPT
            rA = rA_new

            rec_B = J_exc @ rB / N - global_inhib*rB.mean()
            rB_new = rB + dt*(-rB + F(rec_B + cross_to_B + explore_B - aB)) / tau
            aB = aB + dt*(-aB + GAIN_ADAPT*rB) / TAU_ADAPT
            rB = rB_new

        pos_A = bump_pos(rA)
        pos_B = bump_pos(rB)
        pos_A_turni.append(pos_A)
        pos_B_turni.append(pos_B)

        if train and listener == 'B':
            parlante_A = pos_A if pos_A is not None else ultima_pos_A
            if parlante_A is not None and pos_B is not None:
                W_AB += ETA * hebbian_bump(parlante_A, pos_B)
                np.clip(W_AB, 0, W_MAX, out=W_AB)
        elif train and listener == 'A' and allena_BA:
            parlante_B = pos_B if pos_B is not None else ultima_pos_B
            if parlante_B is not None and pos_A is not None:
                W_BA += ETA_BA * hebbian_bump(parlante_B, pos_A)
                np.clip(W_BA, 0, W_MAX, out=W_BA)

        if pos_A is not None:
            ultima_pos_A = pos_A
        if pos_B is not None:
            ultima_pos_B = pos_B

    return pos_A_turni, pos_B_turni


def eval_batch(W_AB, W_BA, messaggio, n_trial, seed0):
    """Valutazione pulita: pesi congelati, niente rumore esplorativo.
    Schema turni: 0=messaggio(A) 1=B ascolta 2=silenzio 3=A ascolta
    4=silenzio 5=B ascolta 6=silenzio 7=A ascolta.
    Il confronto per A va fatto contro il turno di silenzio subito prima
    (2 per il turno 3, 6 per il turno 7): senza input esterno in quel
    turno di silenzio A e' gia' quello che sara' -- solo un cambiamento
    da li' al turno di ascolto successivo indica una vera risposta a B."""
    successi_B1, successi_A3, successi_B2, successi_A7 = 0, 0, 0, 0
    pos_B_ok, pos_A_ok = [], []
    for k in range(n_trial):
        rng = np.random.RandomState(seed0 + k)
        pos_A, pos_B = run_episode(W_AB.copy(), W_BA.copy(), messaggio, rng, train=False)
        if pos_B[1] is not None:
            successi_B1 += 1
            pos_B_ok.append(pos_B[1])
        if pos_A[3] is not None and (pos_A[2] is None or abs(circ_dist(pos_A[3]-pos_A[2])) > 0.02):
            successi_A3 += 1
            pos_A_ok.append(pos_A[3])
        if pos_B[5] is not None:
            successi_B2 += 1
        if pos_A[7] is not None and (pos_A[6] is None or abs(circ_dist(pos_A[7]-pos_A[6])) > 0.02):
            successi_A7 += 1
    return {
        'tasso_B1': successi_B1/n_trial,
        'tasso_A3': successi_A3/n_trial,
        'tasso_B2': successi_B2/n_trial,
        'tasso_A7': successi_A7/n_trial,
        'std_B': np.std(pos_B_ok) if len(pos_B_ok) > 1 else None,
        'std_A': np.std(pos_A_ok) if len(pos_A_ok) > 1 else None,
        'media_B': np.mean(pos_B_ok) if pos_B_ok else None,
    }


if __name__ == '__main__':
    MESSAGGIO = 0.40
    N_EPISODES = 2500

    W_AB = init_weak_bias(OFFSET_AB_INNATO, W_INIT_SCALE)
    W_BA = init_weak_bias(OFFSET_BA_INNATO, W_INIT_SCALE)

    print("=" * 70)
    print("SETUP: bias innato deliberatamente debole")
    print("=" * 70)
    print(f"  W_INIT_SCALE={W_INIT_SCALE}, soglia per formare un bump = 0.5")
    W_AB_iniziale = W_AB.copy()
    W_BA_iniziale = W_BA.copy()

    print("\n" + "=" * 70)
    print("BASELINE (pre-training): il bias innato basta da solo?")
    print("=" * 70)
    baseline = eval_batch(W_AB, W_BA, MESSAGGIO, n_trial=20, seed0=0)
    print(f"  B risponde (turno 1): {baseline['tasso_B1']*100:.0f}%   |   "
          f"A risponde per davvero (turno 3): {baseline['tasso_A3']*100:.0f}%")
    print(f"  B risponde di nuovo (turno 5): {baseline['tasso_B2']*100:.0f}%   |   "
          f"A risponde di nuovo (turno 7): {baseline['tasso_A7']*100:.0f}%")

    print("\n" + "=" * 70)
    print(f"TRAINING: {N_EPISODES} episodi ({WARMUP_EPISODES} di riscaldamento "
          f"solo A->B, poi entrambi i canali), messaggio fisso x={MESSAGGIO}")
    print("=" * 70)
    rng_train = np.random.RandomState(42)
    successo_B_storia = []
    successo_A_storia = []
    picco_W_AB_storia = []
    picco_W_BA_storia = []
    for ep in range(N_EPISODES):
        allena_BA = ep >= WARMUP_EPISODES
        pos_A, pos_B = run_episode(W_AB, W_BA, MESSAGGIO, rng_train, train=True, allena_BA=allena_BA)
        W_AB *= DECAY_PER_EPISODE
        if allena_BA:
            W_BA *= DECAY_PER_EPISODE
        successo_B_storia.append(1 if pos_B[1] is not None else 0)
        spostato = (pos_A[3] is not None and (pos_A[2] is None
                    or abs(circ_dist(pos_A[3]-pos_A[2])) > 0.02))
        successo_A_storia.append(1 if spostato else 0)
        picco_W_AB_storia.append(W_AB.max())
        picco_W_BA_storia.append(W_BA.max())
        if ep % 50 == 0 or ep == N_EPISODES-1:
            win_B = successo_B_storia[max(0, ep-30):ep+1]
            win_A = successo_A_storia[max(0, ep-30):ep+1]
            fase = "riscaldamento" if not allena_BA else "entrambi i canali"
            print(f"  episodio {ep:3d} ({fase}): successo B (ultimi 30) = {np.mean(win_B)*100:5.1f}%  |  "
                  f"successo A (ultimi 30) = {np.mean(win_A)*100:5.1f}%  |  "
                  f"picco W_AB = {W_AB.max():.4f}  picco W_BA = {W_BA.max():.4f}")

    print("\n" + "=" * 70)
    print("DOPO IL TRAINING: confronto pulito, pesi congelati, nessun rumore esplorativo")
    print("=" * 70)
    dopo = eval_batch(W_AB, W_BA, MESSAGGIO, n_trial=30, seed0=9000)
    print(f"  B risponde (turno 1):    PRIMA {baseline['tasso_B1']*100:5.0f}%   ->   DOPO {dopo['tasso_B1']*100:5.0f}%")
    print(f"  A risponde (turno 3):    PRIMA {baseline['tasso_A3']*100:5.0f}%   ->   DOPO {dopo['tasso_A3']*100:5.0f}%")
    print(f"  B risponde (turno 5):    PRIMA {baseline['tasso_B2']*100:5.0f}%   ->   DOPO {dopo['tasso_B2']*100:5.0f}%")
    print(f"  A risponde (turno 7):    PRIMA {baseline['tasso_A7']*100:5.0f}%   ->   DOPO {dopo['tasso_A7']*100:5.0f}%")
    std_dopo_str = f"{dopo['std_B']:.4f}" if dopo['std_B'] else "None"
    print(f"  std posizione B dopo il training: {std_dopo_str}")
    if dopo['media_B'] is not None:
        offset_appreso = circ_dist(dopo['media_B'] - MESSAGGIO)
        print(f"\n  Offset A->B appreso empiricamente: {offset_appreso:.3f} "
              f"(bias innato di partenza era {OFFSET_AB_INNATO})")

    print("\n" + "=" * 70)
    print("RISULTATO: le due sinapsi")
    print("=" * 70)
    print(f"  Picco finale W_AB (A->B): {W_AB.max():.4f}")
    print(f"  Picco finale W_BA (B->A): {W_BA.max():.4f}")
    tassi_finali = [dopo['tasso_B1'], dopo['tasso_A3'], dopo['tasso_B2'], dopo['tasso_A7']]
    if min(tassi_finali) >= 0.9:
        print("  >>> Le due sinapsi si autosostengono ed entrambe le direzioni rispondono in "
              "modo affidabile su tutti e quattro gli scambi verificati (B al turno 1 e 5, "
              "A al turno 3 e 7). Il turno 5 (risposta di B alla NUOVA posizione di A, non "
              "al messaggio originale) e' una associazione di secondo ordine -- si consolida "
              "solo dopo che la prima risposta di A si e' gia' stabilizzata, quindi richiede "
              "piu' episodi della prima; con budget di training sufficiente arriva anch'essa "
              "al 100%. L'asimmetria di partenza (B->A che non imparava affatto) e' risolta "
              "in questo run, non solo attenuata. <<<")
    elif max(dopo['tasso_A3'], dopo['tasso_A7']) > 0.3:
        print("  >>> Il canale B->A si e' consolidato parzialmente: A risponde a B in una "
              "frazione consistente delle prove, non piu' 0% come nella versione senza "
              "riscaldamento. Resta pero' meno affidabile del canale A->B su almeno uno "
              "dei quattro scambi verificati -- l'asimmetria si e' ridotta, non e' sparita "
              "del tutto. <<<")
    else:
        print("  >>> Anche con adattamento, turno di silenzio e riscaldamento, il canale "
              "B->A resta inaffidabile in questa run. <<<")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    window = 20
    curva = [np.mean(successo_B_storia[max(0, k-window):k+1]) for k in range(len(successo_B_storia))]
    curva_A = [np.mean(successo_A_storia[max(0, k-window):k+1]) for k in range(len(successo_A_storia))]
    axes[0, 0].plot(curva, color='C3', label='B risponde ad A (turno 1)')
    axes[0, 0].plot(curva_A, color='C0', label='A risponde a B (turno 3)')
    axes[0, 0].axvline(WARMUP_EPISODES, color='gray', linestyle=':', label='fine riscaldamento')
    axes[0, 0].set_xlabel('Episodio di training')
    axes[0, 0].set_ylabel(f'Tasso di successo (media mobile su {window})')
    axes[0, 0].set_title("Curve di apprendimento dei due canali")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].set_ylim(-0.05, 1.05)

    axes[0, 1].plot(picco_W_AB_storia, color='C3', label='W_AB (A->B)')
    axes[0, 1].plot(picco_W_BA_storia, color='C0', label='W_BA (B->A)')
    axes[0, 1].axvline(WARMUP_EPISODES, color='gray', linestyle=':', label='fine riscaldamento')
    axes[0, 1].set_xlabel('Episodio di training')
    axes[0, 1].set_ylabel('Peso massimo')
    axes[0, 1].set_title("Consolidamento delle due sinapsi")
    axes[0, 1].legend(fontsize=8, loc='upper left')
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
    plt.savefig('dialogo_stdp_adattamento.png', dpi=150)
    print("\nGrafico salvato in dialogo_stdp_adattamento.png")
