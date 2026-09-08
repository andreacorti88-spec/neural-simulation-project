# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
CONNETTIVITA' STRUTTURATA (SMALL-WORLD): verso un realismo spaziale
====================================================================
Tutte le reti costruite finora (anche quella a 40 milioni di neuroni)
usavano connettivita' CASUALE: ogni neurone si collega a K altri
scelti a caso in tutta la rete, senza alcuna nozione di "spazio" o
"vicinanza". Il cervello reale non funziona cosi': un neurone si
connette soprattutto ai suoi vicini fisici, con solo poche
connessioni rare a lunga distanza (le fibre della sostanza bianca).
Questo tipo di organizzazione si chiama "small-world" (Watts-Strogatz,
1998) ed e' uno dei pattern piu' universali nei sistemi biologici.

MODELLO USATO QUI:
  - I 1000 neuroni sono disposti su un ANELLO (posizione x da 0 a 1,
    con avvolgimento circolare)
  - La probabilita' di connessione tra due neuroni DECRESCE
    esponenzialmente con la distanza sull'anello (connessioni dense
    e locali), piu' una piccola probabilita' costante di connessioni
    "a scorciatoia" a lunga distanza (le shortcut del modello
    small-world, che tengono la rete comunque ben connessa nel suo
    insieme)

COSA OSSERVIAMO CHE NON SI VEDEVA CON LA CONNETTIVITA' CASUALE:
Stimolando un piccolo gruppo di neuroni ADIACENTI nello spazio (non
sparsi a caso come prima), il segnale si propaga come una VERA ONDA
che si allontana dal punto di stimolo -- misurabile: la distanza
media degli spike dal punto di stimolo cresce in modo monotono nel
tempo, poi si stabilizza quando l'onda ha investito l'intera rete.
Con connettivita' casuale questo fenomeno non puo' esistere, perche'
non c'e' alcuna nozione di "vicino" o "lontano" nella topologia.
"""

from brian2 import *
import numpy as np

start_scope()

N = 1000
N_E = 800
N_I = 200

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
du/dt = a*(b*v - u) : volt/second
I : volt/second
x : 1
'''

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV + 15*mV*rand(N) - 5*mV
neurons.u = b * neurons.v[:] + 2*mV/ms*randn(N)
neurons.I = 0*mV/ms
neurons.x = arange(N)/N   # posizione sull'anello, distribuita uniformemente

E = neurons[:N_E]
I_pop = neurons[N_E:]

background = PoissonGroup(N, rates=1500*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
bg_syn.connect(j='i')
w_exc = 1.0*mV
w_inh = 4.0*mV

# ---------------------------------------------------------------
# CONNETTIVITA' DIPENDENTE DALLA DISTANZA (il cuore del modello)
# ---------------------------------------------------------------
xi = 0.03           # lunghezza caratteristica: quanto "locale" e' la connettivita'
p_local = 0.5        # probabilita' massima per neuroni immediatamente vicini
p_shortcut = 0.005    # probabilita' costante di connessioni a lunga distanza

# distanza circolare tra due punti su un anello (0..1, con avvolgimento)
dist_expr = '(0.5 - abs(abs(x_pre-x_post) - 0.5))'
prob_expr = f'{p_local}*exp(-({dist_expr})/{xi}) + {p_shortcut}'

syn_ee = Synapses(E, E, on_pre='v_post += w_exc')
syn_ee.connect(condition='i!=j', p=prob_expr)
syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
syn_ei.connect(p=prob_expr)
syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
syn_ie.connect(p=prob_expr)
syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
syn_ii.connect(condition='i!=j', p=prob_expr)

print(f"Sinapsi create -- EE: {len(syn_ee)}, EI: {len(syn_ei)}, "
      f"IE: {len(syn_ie)}, II: {len(syn_ii)}")

# ---------------------------------------------------------------
# BURN-IN + BASELINE
# ---------------------------------------------------------------
run(200*ms)
print("Burn-in completato.")

spikes = SpikeMonitor(neurons)
