"""
Estensione di dialogo_stdp_adattamento.py: la' un solo messaggio fisso
(x=0.40) imparava uno scambio bidirezionale perfettamente affidabile in
entrambe le direzioni. Qui la domanda e' se la STESSA coppia di sinapsi
W_AB/W_BA puo' imparare a distinguere PIU' messaggi diversi senza che si
interferiscano a vicenda -- il vero test di generalizzazione, nello
stesso spirito di contextual_policy.py (che estende rl_bandit.py da un
singolo contesto a piu' contesti).

Messaggi usati: 0.15 e 0.65 (distanza circolare ~0.5, il piu' possibile
separati sull'anello). Ogni episodio ne sceglie uno a caso, cosi' come
contextual_policy.py sceglie un contesto a caso a ogni prova. Le sinapsi
restano matrici N x N condivise fra tutti i messaggi -- se i due messaggi
sono abbastanza separati rispetto alla larghezza del kernel Hebbiano
(SIGMA_W=0.05, quindi le due "zone" occupano regioni ben distinte della
matrice), l'aspettativa e' che NON interferiscano; il codice sotto lo
verifica invece di assumerlo, con una valutazione pulita per messaggio a
fine training e un controllo esplicito sull'ALTRO messaggio ogni volta che
se ne allena uno (per scoprire un eventuale degrado nel tempo).

Riassunto di dialogo_stdp_adattamento.py (di cui questo file eredita tutta
la meccanica): il canale A->B si consolida con successo (0% -> 100%) usando
adattamento a soglia di scarica + rumore esplorativo + rinforzo Hebbiano di
fine turno. Il canale di ritorno B->A era piu' lento (bersaglio mobile nei
primi episodi) ma si risolve con una fase di riscaldamento (solo A->B) e
un tasso di apprendimento piu' alto per B->A. Con un solo messaggio,
tutti e quattro gli scambi verificati (B al turno 1 e 5, A al turno 3 e 7)
arrivano al 100% di affidabilita' in valutazione pulita.

AGGIORNAMENTO: dialogo_multi_messaggio.py aveva gia' risposto alla domanda
di questo file, ma non come sperato. Con due messaggi (0.15 e 0.65) il
canale B->A restava pulito al 100% per entrambi su due run indipendenti
(anche con esposizione garantita identica, round-robin), ma il canale
A->B mostrava un'interferenza competitiva reale: un solo messaggio vinceva
e arrivava al 100%, l'altro crollava a 0% -- e CAMBIANDO run cambiava
quale dei due vinceva, escludendo che fosse solo esposizione ineguale.

Qui si prova una correzione mirata proprio a quel meccanismo: una
NORMALIZZAZIONE OMEOSTATICA del totale dei pesi in ENTRATA su ogni
neurone di B (il totale di ogni colonna di W_AB non puo' superare un
tetto fisso, W_AB_COL_CAP). L'idea, motivata biologicamente (i neuroni
veri regolano omeostaticamente il totale delle proprie sinapsi in
entrata): se i due messaggi competono per la STESSA zona di risposta in
B, la normalizzazione impedisce a uno dei due di saturarla del tutto a
spese dell'altro; se invece puntano a zone diverse, non li tocca. Un
test di regressione con un solo messaggio conferma che la normalizzazione
non interferisce con l'apprendimento quando non c'e' competizione (il
totale di colonna resta ben sotto il tetto per tutto il training).
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
# tentativo di correzione dell'interferenza osservata in dialogo_multi_messaggio.py
# (con esposizione uguale garantita in round-robin, un messaggio arrivava
# comunque al 100% e l'altro crollava a 0% su A->B, mentre B->A restava
# sempre pulito): normalizzazione omeostatica del totale dei pesi in entrata
# su ogni neurone di B. Cap stimato da un run di riferimento a un solo
# messaggio (dialogo_stdp_adattamento.py, picco W_AB=0.0476 -- somma di
# colonna stimata ~picco*sigma_w*sqrt(2*pi)*N =~ 1.2): il cap e' fissato un
# po' sopra quel valore, cosi' un singolo messaggio ha ancora spazio per
# diventare autosufficiente, ma due messaggi che puntano alla stessa zona
# di B non possono piu' saturarla entrambi senza competere per lo stesso
# budget
W_AB_COL_CAP = 1.5
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
# piu' lento di A->B in 2000 episodi). Con piu' messaggi il budget si
# divide tra loro (round-robin), quindi qui e' raddoppiato rispetto alla
# versione a un solo messaggio per dare a CIASCUN messaggio lo stesso
# numero di episodi di riscaldamento (~400) della versione originale
WARMUP_EPISODES = 800
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


def normalizza_colonne(W, cap):
    """Scaling omeostatico: il totale dei pesi in ENTRATA su ogni neurone
    di B (colonna di W_AB) non puo' superare 'cap'. Se due messaggi
    competono per la stessa zona di B, questo impedisce a uno dei due di
    monopolizzare quella colonna a spese dell'altro -- senza impedire a
    messaggi che puntano a zone DIVERSE di crescere ciascuno per conto suo."""
    somme = W.sum(axis=0)
    fattore = np.minimum(1.0, cap / np.maximum(somme, 1e-12))
    W *= fattore[None, :]


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
                normalizza_colonne(W_AB, W_AB_COL_CAP)
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
    MESSAGGI = [0.15, 0.65]
    N_EPISODES = 9000
    EVAL_EVERY = 500
    N_TRIAL_EVAL = 15

    W_AB = init_weak_bias(OFFSET_AB_INNATO, W_INIT_SCALE)
    W_BA = init_weak_bias(OFFSET_BA_INNATO, W_INIT_SCALE)

    print("=" * 70)
    print(f"SETUP: {len(MESSAGGI)} messaggi = {MESSAGGI}, stessa coppia di sinapsi condivisa")
    print("=" * 70)
    W_AB_iniziale = W_AB.copy()
    W_BA_iniziale = W_BA.copy()

    print("\n" + "=" * 70)
    print("BASELINE (pre-training), per ciascun messaggio")
    print("=" * 70)
    for m in MESSAGGI:
        b = eval_batch(W_AB, W_BA, m, n_trial=15, seed0=0)
        print(f"  x={m:.2f}:  B(t1)={b['tasso_B1']*100:.0f}%  A(t3)={b['tasso_A3']*100:.0f}%  "
              f"B(t5)={b['tasso_B2']*100:.0f}%  A(t7)={b['tasso_A7']*100:.0f}%")

    print("\n" + "=" * 70)
    print(f"TRAINING: {N_EPISODES} episodi ({WARMUP_EPISODES} di riscaldamento solo A->B), "
          f"messaggi alternati in round-robin (esposizione garantita identica) tra {MESSAGGI}")
    print("=" * 70)
    rng_train = np.random.RandomState(42)
    picco_W_AB_storia = []
    picco_W_BA_storia = []
    episodi_check = []
    curva_min_per_messaggio = {m: [] for m in MESSAGGI}

    for ep in range(N_EPISODES):
        allena_BA = ep >= WARMUP_EPISODES
        messaggio = MESSAGGI[ep % len(MESSAGGI)]
        run_episode(W_AB, W_BA, messaggio, rng_train, train=True, allena_BA=allena_BA)
        W_AB *= DECAY_PER_EPISODE
        if allena_BA:
            W_BA *= DECAY_PER_EPISODE
        picco_W_AB_storia.append(W_AB.max())
        picco_W_BA_storia.append(W_BA.max())

        if ep % EVAL_EVERY == 0 or ep == N_EPISODES - 1:
            episodi_check.append(ep)
            fase = "riscaldamento" if not allena_BA else "entrambi i canali"
            riga = f"  episodio {ep:4d} ({fase}):"
            for m in MESSAGGI:
                r = eval_batch(W_AB, W_BA, m, n_trial=N_TRIAL_EVAL, seed0=90000 + ep)
                minimo = min(r['tasso_B1'], r['tasso_A3'], r['tasso_B2'], r['tasso_A7'])
                curva_min_per_messaggio[m].append(minimo)
                riga += f"  x={m:.2f} peggiore-dei-4={minimo*100:5.1f}%"
            print(riga)

    print("\n" + "=" * 70)
    print("DOPO IL TRAINING: valutazione pulita per ciascun messaggio")
    print("=" * 70)
    dopo_per_messaggio = {}
    for m in MESSAGGI:
        dopo = eval_batch(W_AB, W_BA, m, n_trial=40, seed0=9000)
        dopo_per_messaggio[m] = dopo
        print(f"  x={m:.2f}:  B(t1)={dopo['tasso_B1']*100:5.0f}%  A(t3)={dopo['tasso_A3']*100:5.0f}%  "
              f"B(t5)={dopo['tasso_B2']*100:5.0f}%  A(t7)={dopo['tasso_A7']*100:5.0f}%")
        if dopo['media_B'] is not None:
            offset = circ_dist(dopo['media_B'] - m)
            print(f"           offset A->B appreso: {offset:.3f}  (bias innato: {OFFSET_AB_INNATO})")

    print("\n" + "=" * 70)
    print("RISULTATO: interferenza tra i messaggi?")
    print("=" * 70)
    print(f"  Picco finale W_AB: {W_AB.max():.4f}   Picco finale W_BA: {W_BA.max():.4f}")
    tutti_i_tassi = [v for d in dopo_per_messaggio.values()
                      for v in (d['tasso_B1'], d['tasso_A3'], d['tasso_B2'], d['tasso_A7'])]
    offsets = [circ_dist(dopo_per_messaggio[m]['media_B'] - m) for m in MESSAGGI
               if dopo_per_messaggio[m]['media_B'] is not None]
    offset_spread = (max(offsets) - min(offsets)) if len(offsets) > 1 else None
    if min(tutti_i_tassi) >= 0.9:
        print("  >>> Nessuna interferenza rilevabile: entrambi i messaggi, con la STESSA coppia "
              "di sinapsi condivisa, arrivano al 100% (o quasi) su tutti e quattro gli scambi. "
              f"Gli offset appresi per i due messaggi sono {[f'{o:.3f}' for o in offsets]} -- "
              f"differiscono di {offset_spread:.3f} tra loro, se fosse pura sovrapposizione "
              "casuale ci aspetteremmo che collassassero sullo stesso valore o si "
              "distruggessero a vicenda; non e' successo. <<<")
    elif min(tutti_i_tassi) < 0.5:
        print("  >>> Interferenza rilevabile: almeno un messaggio non ha raggiunto un livello "
              "affidabile -- puo' darsi che il budget di training condiviso non sia ancora "
              "sufficiente per consolidare piu' associazioni (ciascun messaggio riceve solo "
              "una frazione degli episodi totali), non necessariamente vera interferenza "
              "distruttiva. Verificare le curve per-messaggio nel grafico. <<<")
    else:
        print("  >>> Risultato misto: alcuni scambi sono affidabili, altri no -- vedi il "
              "dettaglio sopra per capire quali. <<<")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    colori = ['C3', 'C0', 'C2', 'C4']
    for i, m in enumerate(MESSAGGI):
        axes[0, 0].plot(episodi_check, curva_min_per_messaggio[m], 'o-',
                         color=colori[i % len(colori)], label=f'messaggio x={m:.2f}')
    axes[0, 0].axvline(WARMUP_EPISODES, color='gray', linestyle=':', label='fine riscaldamento')
    axes[0, 0].set_xlabel('Episodio di training')
    axes[0, 0].set_ylabel('Tasso di successo del PEGGIORE dei 4 scambi')
    axes[0, 0].set_title("Un messaggio degrada l'altro? (valutazione pulita, periodica)")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].set_ylim(-0.05, 1.05)

    axes[0, 1].plot(picco_W_AB_storia, color='C3', label='W_AB (A->B)')
    axes[0, 1].plot(picco_W_BA_storia, color='C0', label='W_BA (B->A)')
    axes[0, 1].axvline(WARMUP_EPISODES, color='gray', linestyle=':', label='fine riscaldamento')
    axes[0, 1].set_xlabel('Episodio di training')
    axes[0, 1].set_ylabel('Peso massimo')
    axes[0, 1].set_title("Consolidamento delle due sinapsi (condivise fra i messaggi)")
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
    axes[1, 1].set_title(f'W_AB DOPO il training -- righe attese vicino a x={MESSAGGI}')
    plt.colorbar(im1, ax=axes[1, 1])

    plt.tight_layout()
    plt.savefig('dialogo_multi_normalizzato.png', dpi=150)
    print("\nGrafico salvato in dialogo_multi_normalizzato.png")
