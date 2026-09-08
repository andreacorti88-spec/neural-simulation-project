# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
SINTESI FINALE: struttura spaziale + sinapsi realistiche + scala massima
============================================================================
Questo esperimento unisce i tre miglioramenti sviluppati separatamente nel
progetto in un'unica rete:

  1. STRUTTURA SPAZIALE (sezione 4.1-4.2): i neuroni sono disposti in una
     forma 3D a due lobi (emisferi), con connettivita' che dipende dalla
     distanza -- densa localmente, rara a lungo raggio -- calcolata
     efficientemente con un KD-tree

  2. SINAPSI REALISTICHE (esperimento precedente): non piu' semplici "salti
     di tensione", ma sinapsi a conduttanza con potenziale di inversione,
     il modello biofisico standard in letteratura

  3. DENSITA' SINAPTICA REALISTICA: ~2600 sinapsi per neurone, la stessa
    densita' media del modello Allen Institute (il piu' avanzato al mondo)

La scala (118.000 neuroni) resta la stessa del precedente esperimento con
sinapsi realistiche, perche' e' dettata dallo stesso vincolo di memoria
(circa 307 milioni di sinapsi e' il tetto pratico sui tuoi 24 GB di RAM
con questo tipo di sinapsi) -- qui la STRUTTURA delle connessioni cambia
(spaziale invece che casuale), non il numero totale.

AVVERTENZA IMPORTANTE (dall'esperimento precedente): questa densita' di
connessione produce un\'attivita\' di rete che NON si stabilizza mai del
tutto in uno stato quieto -- l\'ampiezza delle oscillazioni spontanee
cresce di nuovo dopo essersi inizialmente smorzata, un comportamento
verificato su piu' prove indipendenti. Non e\' un errore, e' una proprieta\'
osservata di reti con questa densita\' di connessione ricorrente.

TEMPO ATTESO: in base ai test di calibrazione, l'esperimento completo
(costruzione + burn-in + baseline + stimolo + recovery) dovrebbe richiedere
tra i 25 e i 50 minuti, in base a quanto interviene il throttling termico.
"""

from brian2 import *
import numpy as np
from scipy.spatial import cKDTree
import time as pytime

print("Fase 1/4: generazione della forma 3D del cervello (118.000 punti)...")



def sample_brain_shape(n_target):
    points = []
    total = 0
    while total < n_target:
        batch = np.random.uniform(-1.3, 1.3, size=(n_target*3, 3))
        x, y, z = batch[:, 0], batch[:, 1], batch[:, 2]
        left = ((x+0.35)/0.55)**2 + (y/0.45)**2 + (z/0.42)**2
        right = ((x-0.35)/0.55)**2 + (y/0.45)**2 + (z/0.42)**2
        inside = (left <= 1.0) | (right <= 1.0)
        pts = batch[inside]
        points.append(pts)
        total += len(pts)
    return np.vstack(points)[:n_target]



N = 118_000
N_E = int(N * 0.8)
K_LOCAL = 2200   # connessioni verso i vicini spaziali piu' prossimi
K_FAR = 400      # connessioni "scorciatoia" a lunga distanza (small-world)

np.random.seed(42)
pts = sample_brain_shape(N)
print(f"Punti generati.")

print("Fase 2/4: calcolo della connettivita' spaziale (KD-tree, a blocchi "
      "per contenere l'uso di memoria)...")
t0 = pytime.time()
tree = cKDTree(pts)
batch_size = 2000
all_neigh = np.zeros((N, K_LOCAL), dtype=np.int32)
for start in range(0, N, batch_size):
    end = min(start + batch_size, N)
    _, neigh = tree.query(pts[start:end], k=K_LOCAL+1, workers=-1)
    all_neigh[start:end] = neigh[:, 1:].astype(np.int32)
print(f"Connettivita' locale calcolata in {pytime.time()-t0:.1f}s.")
src_local = np.repeat(np.arange(N), K_LOCAL)
tgt_local = all_neigh.flatten()
src_far = np.repeat(np.arange(N), K_FAR)
tgt_far = np.random.randint(0, N, size=N*K_FAR)
src_all = np.concatenate([src_local, src_far])
tgt_all = np.concatenate([tgt_local, tgt_far])
valid = src_all != tgt_all
src_all, tgt_all = src_all[valid], tgt_all[valid]
del src_local, tgt_local, src_far, tgt_far, all_neigh     # libera memoria

print(f"Connessioni totali generate: {len(src_all):,} "
      f"({len(src_all)/N:.0f} per neurone)")

print("Fase 3/4: costruzione della rete Brian2...")
start_scope()

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms
E_exc = 0*mV
E_inh = -80*mV
tau_e = 5*ms
tau_i = 10*ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + ge*(E_exc-v) + gi*(E_inh-v) + I_ext : volt
du/dt = a*(b*v - u) : volt/second
dge/dt = -ge/tau_e : 1/second
dgi/dt = -gi/tau_i : 1/second
I_ext : volt/second
'''

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV + 10*mV*rand(N)
neurons.u = b * neurons.v[:]
neurons.I_ext = 0*mV/ms

background = PoissonGroup(N, rates=1200*Hz)
bg_syn = Synapses(background, neurons, on_pre='ge_post += 15/second')
bg_syn.connect(j='i')

w_e = 0.12/second   # ridotto rispetto al tentativo precedente (0.22): a piena
                      # scala e con clustering spaziale, i pesi calibrati sulla
                      # rete puramente casuale risultavano troppo forti e
                      # producevano oscillazioni crescenti invece di stabilizzarsi
w_i = 0.6/second     # ridotto in proporzione (da 0.9)

src_is_E = src_all < N_E
