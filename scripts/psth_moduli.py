# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
PROVE RIPETUTE E MEDIA (PSTH): la propagazione e' reale o e' rumore?
====================================================================
Nello script precedente, una singola prova sembrava mostrare una
chiara risposta del Modulo 2 allo stimolo del Modulo 1. Ma guardando
il grafico con attenzione, quel rialzo era paragonabile alle normali
fluttuazioni spontanee della rete — non potevamo essere sicuri che
fosse un vero effetto o solo una coincidenza casuale.

LA SOLUZIONE STANDARD in neuroscienza: ripetere l'esperimento molte
volte, allineare ogni prova rispetto al momento dello stimolo, e
fare la MEDIA. Se c'e' un vero effetto, nella media emerge sopra il
rumore (che tende a cancellarsi mediando su piu' prove). Se non c'e'
nessun vero effetto, la media resta piatta. Questa tecnica si chiama
PSTH (Peristimulus Time Histogram) ed e' lo standard per analizzare
dati neurali reali.

In questo script:
  1. Ripetiamo lo stimolo al Modulo 1 per 25 volte, spaziate nel tempo
  2. Per ogni ripetizione, registriamo l'attivita' del Modulo 2 nella
     finestra [-50, +150] ms attorno allo stimolo
  3. Facciamo la media (e calcoliamo l'errore standard) su tutte le
     25 ripetizioni
  4. Confrontiamo il risultato mediato con la variabilita' naturale
     di una singola prova, per capire se l'effetto e' reale
"""

from brian2 import *
import numpy as np

start_scope()

# ---------------------------------------------------------------
# STRUTTURA DELLA RETE (identica allo script dei due moduli)
# ---------------------------------------------------------------
N_E_mod = 120
N_I_mod = 30
N_mod = N_E_mod + N_I_mod
N = N_mod * 2

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

M1_E = neurons[0:120]; M1_I = neurons[120:150]
M2_E = neurons[150:270]; M2_I = neurons[270:300]

background = PoissonGroup(N, rates=1500*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
bg_syn.connect(j='i')

p_local = 0.1
w_exc = 1.0*mV
w_inh = 4.0*mV

def connect_module(E, I):
    see = Synapses(E, E, on_pre='v_post += w_exc')
    see.connect(condition='i!=j', p=p_local)
    sei = Synapses(E, I, on_pre='v_post += w_exc')
    sei.connect(p=p_local)
    sie = Synapses(I, E, on_pre='v_post -= w_inh')
    sie.connect(p=p_local)
    sii = Synapses(I, I, on_pre='v_post -= w_inh')
    sii.connect(condition='i!=j', p=p_local)
    return see, sei, sie, sii

syn1 = connect_module(M1_E, M1_I)
syn2 = connect_module(M2_E, M2_I)

p_inter = 0.02
w_inter = 1.0*mV
syn_12 = Synapses(M1_E, M2_E, on_pre='v_post += w_inter')
syn_12.connect(p=p_inter)
syn_21 = Synapses(M2_E, M1_E, on_pre='v_post += w_inter')
syn_21.connect(p=p_inter)

# ---------------------------------------------------------------
# BURN-IN
# ---------------------------------------------------------------
run(200*ms)
print("Burn-in completato.")

spikes = SpikeMonitor(neurons)

# ---------------------------------------------------------------
# 25 RIPETIZIONI DELLO STIMOLO, spaziate di 300ms l'una dall'altra
# (abbastanza tempo perche' la rete torni alla normalita' tra una
# ripetizione e la successiva)
# ---------------------------------------------------------------
N_TRIALS = 25
ISI = 300*ms           # tempo tra uno stimolo e il successivo
stim_dur = 20*ms
stim_group = M1_E[:20]

stim_times = []
run(50*ms)
for trial in range(N_TRIALS):
    stim_times.append(float(defaultclock.t/ms))
    stim_group.I = 40*mV/ms
    run(stim_dur)
    stim_group.I = 0*mV/ms
    run(ISI - stim_dur)

print(f"Completate {N_TRIALS} ripetizioni dello stimolo.")

# ---------------------------------------------------------------
# COSTRUZIONE DEL PSTH (media sulle prove)
# ---------------------------------------------------------------
t_spikes = np.array(spikes.t/ms)
i_spikes = np.array(spikes.i)
mod2_mask = i_spikes >= 150
mod1_mask = i_spikes < 150

window_pre, window_post, bin_size = 50, 150, 5
bins = np.arange(-window_pre, window_post, bin_size)
bin_centers = bins[:-1] + bin_size/2

