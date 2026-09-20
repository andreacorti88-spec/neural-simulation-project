"""
Unione dei tre assi sviluppati separatamente: struttura spaziale 3D a due
lobi con connettivita' dipendente dalla distanza (KD-tree, densa in loco +
shortcut a lungo raggio), sinapsi a conduttanza (rete_realistica_118K), densita'
~2600 sinapsi/neurone livello Allen Institute. Scala invariata a 118K rispetto
all'esperimento precedente (stesso tetto di memoria, ~307M sinapsi su 24GB) --
qui cambia la struttura delle connessioni (spaziale vs casuale), non il numero.

Nota: a questa densita' l'attivita' non si stabilizza mai del tutto -- le
oscillazioni spontanee ricrescono dopo uno smorzamento iniziale, comportamento
riproducibile su piu' run indipendenti, non un errore di calibrazione.
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
del src_local, tgt_local, src_far, tgt_far, all_neigh

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

# ridotti rispetto a rete_realistica_118K (0.22/0.9): con clustering spaziale
# i pesi calibrati sulla rete casuale erano troppo forti, oscillazioni crescenti
w_e = 0.12/second
w_i = 0.6/second

src_is_E = src_all < N_E
syn_exc = Synapses(neurons, neurons, on_pre='ge_post += w_e')
syn_exc.connect(i=src_all[src_is_E], j=tgt_all[src_is_E])
syn_inh = Synapses(neurons, neurons, on_pre='gi_post += w_i')
syn_inh.connect(i=src_all[~src_is_E], j=tgt_all[~src_is_E])

n_syn_total = len(syn_exc) + len(syn_inh) + len(bg_syn)
print(f"Rete costruita. Sinapsi totali: {n_syn_total:,} "
      f"({n_syn_total/N:.0f}/neurone)")

net = Network(collect())

print("\nFase 4/4: simulazione (burn-in + baseline + stimolo + recovery)...")
t0 = pytime.time()
net.run(200*ms)
print(f"Burn-in completato in {pytime.time()-t0:.1f}s.")

rate_mon = PopulationRateMonitor(neurons)
raster_subset = SpikeMonitor(neurons[:1500], record=True)
net.add(rate_mon, raster_subset)

t0 = pytime.time()
net.run(200*ms)
print(f"Baseline completata in {pytime.time()-t0:.1f}s.")

# stimolo localizzato vicino alla punta dell'emisfero sinistro
target = np.array([-0.85, 0.0, 0.0])
d_target = np.linalg.norm(pts[:N_E] - target, axis=1)
stim_idx = np.argsort(d_target)[:200]
stim_time = 400*ms
neurons.I_ext[stim_idx] = 40*mV/ms
net.run(15*ms)
neurons.I_ext[stim_idx] = 0*mV/ms
print(f"Stimolo applicato a {len(stim_idx)} neuroni vicino alla punta "
      f"dell'emisfero sinistro.")

t0 = pytime.time()
net.run(200*ms)
print(f"Recovery completata in {pytime.time()-t0:.1f}s.")

r = np.array(rate_mon.smooth_rate(window='flat', width=5*ms)/Hz)
t = np.array(rate_mon.t/ms)

baseline_rate = r[(t > 200) & (t < 400)].mean()
peak_rate = r[(t >= 400) & (t < 450)].max()
recovery_rate = r[t > 590].mean()
amplitude_early = r[(t > 200) & (t < 300)]
amplitude_late = r[(t > 500) & (t < 615)]

print("\n" + "=" * 60)
print("RISULTATI FINALI")
print("=" * 60)
print(f"Scala: {N:,} neuroni, {n_syn_total:,} sinapsi ({n_syn_total/N:.0f}/neurone)")
print(f"Struttura: spaziale (small-world 3D, forma a cervello)")
print(f"Sinapsi: a conduttanza (biofisicamente realistiche)")
print(f"Frequenza baseline: {baseline_rate:.2f} Hz")
print(f"Picco dopo lo stimolo: {peak_rate:.2f} Hz")
print(f"Frequenza dopo il recovery: {recovery_rate:.2f} Hz")
print(f"Ampiezza oscillazione (200-300ms): "
      f"{amplitude_early.max()-amplitude_early.min():.2f} Hz")
print(f"Ampiezza oscillazione (500-615ms): "
      f"{amplitude_late.max()-amplitude_late.min():.2f} Hz")

figure(figsize=(11, 7))

subplot(2, 1, 1)
plot(raster_subset.t/ms, raster_subset.i, '.', color='C0', markersize=1.5)
axvspan(stim_time/ms, (stim_time+15*ms)/ms, color='red', alpha=0.15)
ylabel('Indice neurone\n(sottoinsieme di 1500)')
title(f'Rete a {N:,} neuroni -- spaziale + sinapsi realistiche + '
      f'{n_syn_total/N:.0f} sinapsi/neurone')

subplot(2, 1, 2)
plot(t, r, color='black', linewidth=1)
axvspan(stim_time/ms, (stim_time+15*ms)/ms, color='red', alpha=0.15, label='stimolo')
axhline(baseline_rate, color='gray', linestyle='--', linewidth=1, label='baseline')
xlabel('Tempo (ms)')
ylabel('Frequenza popolazione (Hz)')
title('Risposta collettiva della rete allo stimolo')
legend(loc='upper right', fontsize=8)
grid(alpha=0.3)

tight_layout()
savefig('rete_sintesi_finale.png', dpi=150)
print("\nGrafico salvato in rete_sintesi_finale.png")
