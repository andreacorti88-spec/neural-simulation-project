# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
ESPERIMENTO FINALE: 40 MILIONI DI NEURONI (57.1% di un cervello di topo)
================================================================================
Questo e' il traguardo del percorso di scaling: la stessa struttura di
rete bilanciata E/I usata per l'esperimento a 300 neuroni (baseline,
stimolo mirato, recupero), ma alla scala massima verificata come
praticabile sul tuo Mac -- 40 milioni di neuroni, 57.1% della scala
di un cervello di topo (70 milioni).

ACCORGIMENTI PER LA MEMORIA (importanti a questa scala):
  - NON registriamo i tempi di spike di tutti i 40 milioni di neuroni
    (sarebbe troppa memoria) -- registriamo solo:
      1. La frequenza di scarica dell'INTERA popolazione nel tempo
         (PopulationRateMonitor -- efficiente, non salva ogni singolo
         spike, solo l'andamento aggregato)
      2. Il raster dettagliato di un piccolo SOTTOINSIEME di 2000
         neuroni (giusto per avere un'immagine visiva, come le reti
         precedenti)
  - Il numero di neuroni stimolati e' lo 0.5% degli eccitatori
    (~160.000 su 32 milioni) -- una frazione piccola ma sufficiente
    a produrre un effetto visibile sulla frequenza dell'intera rete

TEMPO ATTESO: in base ai dati di scaling raccolti, костruire la rete
(i 40M neuroni + ~1.2 miliardi di sinapsi) richiede circa 9-10 minuti.
Le fasi di baseline/stimolo/recovery aggiungono probabilmente altri
10-20 minuti. Preventiva un totale di 20-30 minuti, con le altre app
chiuse per lasciare piu' RAM libera possibile.
"""

from brian2 import *
import time as pytime
import numpy as np

print("Costruzione della rete (40 milioni di neuroni)...")
print("Questo passaggio da solo richiede diversi minuti. Non interrompere.")

start_scope()

N = 40_000_000
N_E = int(N * 0.8)
N_I = N - N_E

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
du/dt = a*(b*v - u) : volt/second
I : volt/second
'''

t_build_start = pytime.time()

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV + 15*mV*rand(N) - 5*mV
neurons.u = b * neurons.v[:] + 2*mV/ms*randn(N)
neurons.I = 0*mV/ms

E = neurons[:N_E]
I_pop = neurons[N_E:]

background = PoissonGroup(N, rates=1500*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
bg_syn.connect(j='i')

K_local = 15
w_exc = 1.0*mV
w_inh = 4.0*mV
K_ei = min(K_local, N_I)
K_ii = min(K_local, N_I)

syn_ee = Synapses(E, E, on_pre='v_post += w_exc')
syn_ee.connect(j=f'k for k in sample(N_E, size={K_local})', skip_if_invalid=True)
syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
syn_ei.connect(j=f'k for k in sample(N_I, size={K_ei})', skip_if_invalid=True)
syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
syn_ie.connect(j=f'k for k in sample(N_E, size={K_local})', skip_if_invalid=True)
syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
syn_ii.connect(j=f'k for k in sample(N_I, size={K_ii})', skip_if_invalid=True)

print(f"Rete costruita in {pytime.time()-t_build_start:.1f}s. Inizio burn-in...")

# ---------------------------------------------------------------
# BURN-IN
# ---------------------------------------------------------------
run(200*ms)
print("Burn-in completato.")

# ---------------------------------------------------------------
# MONITOR (leggeri, per non esaurire la memoria)
# ---------------------------------------------------------------
rate_mon = PopulationRateMonitor(neurons)
raster_subset = SpikeMonitor(neurons[:2000], record=True) # solo 2000 su 40M

# ---------------------------------------------------------------
# BASELINE
# ---------------------------------------------------------------
run(200*ms)
print("Baseline registrata.")
# ---------------------------------------------------------------
# STIMOLO: 0.5% dei neuroni eccitatori (~160.000 su 32 milioni)
# ---------------------------------------------------------------
N_STIM = int(N_E * 0.005)
stim_group = E[:N_STIM]
stim_time = 400*ms
stim_group.I = 40*mV/ms
run(20*ms)
stim_group.I = 0*mV/ms
print(f"Stimolo applicato a {N_STIM:,} neuroni ({N_STIM/N*100:.3f}% della rete totale).")

# ---------------------------------------------------------------
# RECOVERY
# ---------------------------------------------------------------
run(280*ms)
print("Recovery completato.")

# ---------------------------------------------------------------
# ANALISI
# ---------------------------------------------------------------
r = np.array(rate_mon.smooth_rate(window='flat', width=5*ms)/Hz)
t = np.array(rate_mon.t/ms)

baseline_rate = r[(t > 200) & (t < 400)].mean()
peak_rate = r[(t >= 400) & (t < 450)].max()
recovery_rate = r[t > 650].mean()

