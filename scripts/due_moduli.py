# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
DUE MODULI COLLEGATI: verso un'architettura con "regioni"
====================================================================
Finora avevamo una singola rete di 300 neuroni, tutti potenzialmente
collegabili tra loro con la stessa probabilita'. Qui invece costruiamo
due "regioni" distinte:
  MODULO 1: 150 neuroni (120 eccitatori + 30 inibitori)
  MODULO 2: 150 neuroni (120 eccitatori + 30 inibitori)

Dentro ogni modulo, la connettivita' e' densa (10%, come nella rete
singola di prima). TRA i due moduli, la connettivita' e' molto piu'
debole (2%) e SOLO eccitatoria (solo E->E) — questo rispecchia come
funzionano davvero le connessioni a lungo raggio nella corteccia
cerebrale reale: dense localmente, sparse tra regioni distanti.

ESPERIMENTO: stimoliamo SOLO il Modulo 1 e osserviamo:
  - Come risponde il Modulo 1 (dovrebbe attivarsi molto, essendo
    stimolato direttamente)
  - Se e quanto il Modulo 2 (che NON riceve stimolo diretto) si
    attiva comunque, tramite le poche connessioni che lo legano
    al Modulo 1 — questo e' il concetto di PROPAGAZIONE TRA REGIONI
"""

from brian2 import *

start_scope()

# ---------------------------------------------------------------
# STRUTTURA DEI DUE MODULI
# ---------------------------------------------------------------
N_E_mod = 120
N_I_mod = 30
N_mod = N_E_mod + N_I_mod   # 150 per modulo
N = N_mod * 2               # 300 totali

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

# indici: Modulo 1 = [0:150), Modulo 2 = [150:300)
# dentro ogni modulo: prima gli eccitatori, poi gli inibitori
M1_E = neurons[0:120]
M1_I = neurons[120:150]
M2_E = neurons[150:270]
M2_I = neurons[270:300]

# ---------------------------------------------------------------
# RUMORE DI FONDO (uguale per tutti, come prima)
# ---------------------------------------------------------------
background = PoissonGroup(N, rates=1500*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
bg_syn.connect(j='i')

# ---------------------------------------------------------------
# CONNETTIVITA' DENTRO OGNI MODULO (densa, 10%, E e I come prima)
# ---------------------------------------------------------------
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

# ---------------------------------------------------------------
# CONNETTIVITA' TRA I DUE MODULI: debole (2%) e solo eccitatoria
# ---------------------------------------------------------------
p_inter = 0.02
w_inter = 1.0*mV

syn_12 = Synapses(M1_E, M2_E, on_pre='v_post += w_inter')
syn_12.connect(p=p_inter)
syn_21 = Synapses(M2_E, M1_E, on_pre='v_post += w_inter')
syn_21.connect(p=p_inter)

# ---------------------------------------------------------------
# BURN-IN + BASELINE (come nello script precedente)
# ---------------------------------------------------------------
run(200*ms)
print("Burn-in completato.")

rate_mon1 = PopulationRateMonitor(neurons[0:150])     # Modulo 1
rate_mon2 = PopulationRateMonitor(neurons[150:300])   # Modulo 2
spikes = SpikeMonitor(neurons)

run(200*ms)
