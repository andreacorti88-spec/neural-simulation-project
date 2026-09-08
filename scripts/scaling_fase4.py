# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
SCALING FASE 4 (FINALE): l'ultimo tratto verso la scala di un cervello di topo
================================================================================
Finora, passo dopo passo, hai raggiunto 40 milioni di neuroni (57.1%
della scala di un cervello di topo — 70 milioni), usando solo 15.7 GB
di RAM. La crescita della memoria sta rallentando rispetto alle
stime iniziali (i dati reali continuano a essere piu' bassi del
previsto), quindi ricalibro ancora una volta la stima usando TUTTI e
9 i punti raccolti finora, con un margine di sicurezza ridotto al 10%.

Questo e' probabilmente l'ultimo passo di scaling puro: da qui in poi,
se la memoria lo permette, l'obiettivo e' arrivare il piu' vicino
possibile a 70.000.000 (un cervello di topo completo), sapendo che
potremmo non arrivarci per intero.
"""

from brian2 import *
import time as pytime
import numpy as np

try:
    import resource
    HAS_RESOURCE = True
except ImportError:
    HAS_RESOURCE = False



def run_network(N, K_local=15, sim_time=200*ms):
    start_scope()
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

   spikes = SpikeMonitor(neurons, record=False)

   t0 = pytime.time()
   run(sim_time)
   elapsed = pytime.time() - t0

   n_synapses = len(syn_ee) + len(syn_ei) + len(syn_ie) + len(syn_ii) + len(bg_syn)
   n_spikes = int(sum(spikes.count))

   mem_mb = None
   if HAS_RESOURCE:
       mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024*1024)

   return elapsed, n_synapses, n_spikes, mem_mb



# ---------------------------------------------------------------
# TUTTI I DATI RACCOLTI FINORA (fasi 1, 2, 3 — 9 punti reali)
# ---------------------------------------------------------------
N_observed = np.array([1, 3, 5, 10, 20, 25, 30, 35, 40]) * 1_000_000
mem_observed = np.array([1720, 4152, 6254, 9838, 11147, 11495, 12434, 15104, 15744])

slope_mb_per_neuron, intercept = np.polyfit(N_observed, mem_observed, 1)
SAFETY_FACTOR = 1.10

def predict_memory_mb(N):
    return SAFETY_FACTOR * (intercept + slope_mb_per_neuron * N)

MOUSE_BRAIN_NEURONS = 70_000_000
