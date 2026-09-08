# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
RETE A 300 NEURONI: assestamento, attivita' spontanea e risposta a uno stimolo
================================================================================
Tre miglioramenti rispetto allo script precedente:

1. BURN-IN: facciamo "girare a vuoto" la rete per 200ms prima di iniziare
   a registrare, cosi' il forte transitorio sincronizzato iniziale (dovuto
   al fatto che tutti i neuroni partono da condizioni simili) si esaurisce
   e vediamo solo l'attivita' spontanea "vera" a regime.

2. SIMULAZIONE PIU' LUNGA: registriamo un periodo di attivita' spontanea
   (baseline) abbastanza lungo da poter dire se ci sono oscillazioni
   periodiche reali o solo fluttuazioni casuali.

3. STIMOLO ESTERNO MIRATO: dopo la fase di baseline, diamo un impulso
   forte a un piccolo sottoinsieme di 20 neuroni eccitatori (su 300) e
   osserviamo cosa succede al resto della rete: si propaga? si smorza?
   la rete torna al livello di attivita' di partenza?

Questo terzo punto e' concettualmente il piu' vicino a un "piccolo
cervello che reagisce a uno stimolo": input esterno -> propagazione
nella rete -> (si spera) ritorno a uno stato di equilibrio, invece di
esplodere o spegnersi.
"""

from brian2 import *

start_scope()

# ---------------------------------------------------------------
# DIMENSIONI DELLA RETE (identiche allo script precedente)
# ---------------------------------------------------------------
N_E = 240
N_I = 60
N = N_E + N_I

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
# variabilita' iniziale AMPIA (piu' di prima) per rompere meglio la simmetria
neurons.v = -65*mV + 15*mV*rand(N) - 5*mV
neurons.u = b * neurons.v[:] + 2*mV/ms*randn(N)
neurons.I = 0*mV/ms

E = neurons[:N_E]
I_pop = neurons[N_E:]

background = PoissonGroup(N, rates=1500*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
bg_syn.connect(j='i')

p_conn = 0.1
w_exc = 1.0*mV
w_inh = 4.0*mV
syn_ee = Synapses(E, E, on_pre='v_post += w_exc')
syn_ee.connect(condition='i!=j', p=p_conn)
syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
syn_ei.connect(p=p_conn)
syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
syn_ie.connect(p=p_conn)
syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
syn_ii.connect(condition='i!=j', p=p_conn)

# ---------------------------------------------------------------
# FASE 1: BURN-IN (200ms, nessun monitor collegato ancora)
# ---------------------------------------------------------------
run(200*ms)
print("Burn-in completato (200ms, non registrato).")

# ---------------------------------------------------------------
# FASE 2: BASELINE — colleghiamo i monitor e registriamo 200ms
# di attivita' spontanea, prima di dare qualsiasi stimolo
# ---------------------------------------------------------------
rate_mon = PopulationRateMonitor(neurons)
spikes = SpikeMonitor(neurons)

run(200*ms)
print("Baseline registrata (200ms).")

# ---------------------------------------------------------------
# FASE 3: STIMOLO — 20 neuroni eccitatori (su 300) ricevono un
# impulso di corrente forte per 20ms, poi torna tutto normale
# ---------------------------------------------------------------
N_STIM = 20
stim_group = E[:N_STIM]
stim_group.I = 40*mV/ms
stim_time = 400*ms   # tempo assoluto a cui parte lo stimolo (200 burn-in + 200 baseline)

run(20*ms)
stim_group.I = 0*mV/ms
print(f"Stimolo applicato a {N_STIM} neuroni eccitatori (20ms di impulso).")

# ---------------------------------------------------------------
# FASE 4: RECOVERY — osserviamo se e come la rete torna alla
# normalita' dopo lo stimolo
# ---------------------------------------------------------------
run(280*ms)
print("Fase di recovery completata.")

