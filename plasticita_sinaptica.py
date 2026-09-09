"""
PLASTICITA' SINAPTICA: facilitazione a breve termine
========================================================
Finora la sinapsi tra A e B aveva sempre la stessa intensita'.
Qui la facciamo CAMBIARE in base all'uso: e' il primo passo,
concettualmente, verso l'apprendimento biologico.

Il meccanismo che implementiamo si chiama "facilitazione a breve
termine" (short-term facilitation, STF): ogni volta che A spara,
l'efficacia della sinapsi aumenta un po' (variabile 'g'), ma poi
decade esponenzialmente nel tempo se A non spara piu'.
Risultato: se A spara ripetutamente e velocemente, ogni spike
successivo ha un effetto PIU' FORTE su B rispetto al primo.

Questo e' il meccanismo reale dietro fenomeni come la sommazione
temporale nei neuroni corticali.
"""

from brian2 import *

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
du/dt = a*(b*v - u) : volt/second
I : volt/second
'''

start_scope()

neurons = NeuronGroup(2, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV
neurons.u = b * neurons.v

# Stimoliamo A con una corrente PIU' FORTE del solito, cosi' spara
# a raffica e possiamo vedere bene l'effetto di facilitazione.
neurons.I = [25*mV/ms, 0*mV/ms]

# ---------------------------------------------------------------
# SINAPSI CON PLASTICITA'
# ---------------------------------------------------------------
# g = "peso" attuale della sinapsi (parte da un valore base)
# tau_facilitation = quanto velocemente il rafforzamento decade nel tempo
# facilitation_increment = di quanto aumenta 'g' ad ogni spike di A
synapse_eqs = '''
dg/dt = (g_base - g) / tau_facilitation : volt (clock-driven)
g_base : volt
tau_facilitation : second
facilitation_increment : volt
'''

syn = Synapses(neurons, neurons, model=synapse_eqs,
                on_pre='''
                v_post += g
                g += facilitation_increment
                ''',
                method='euler')
syn.connect(i=0, j=1)
syn.delay = 1.5*ms
syn.g = 8*mV                      # peso iniziale (sotto soglia da solo)
syn.g_base = 8*mV                 # valore a cui 'g' decade se A non spara
syn.tau_facilitation = 50*ms      # tempo caratteristico di decadimento
syn.facilitation_increment = 6*mV # quanto si rafforza ad ogni spike

# ---------------------------------------------------------------
# MONITORAGGIO E SIMULAZIONE
# ---------------------------------------------------------------
mon = StateMonitor(neurons, 'v', record=True)
g_mon = StateMonitor(syn, 'g', record=0)
spikes = SpikeMonitor(neurons)

run(200*ms)

# ---------------------------------------------------------------
# VISUALIZZAZIONE
# ---------------------------------------------------------------
figure(figsize=(10, 7))

subplot(3, 1, 1)
plot(mon.t/ms, mon.v[0]/mV, label='Neurone A', color='C0')
plot(mon.t/ms, mon.v[1]/mV, label='Neurone B', color='C1')
ylabel('Potenziale (mV)')
legend(loc='upper right', fontsize=8)
title('Potenziale di membrana')

subplot(3, 1, 2)
plot(g_mon.t/ms, g_mon.g[0]/mV, color='green')
ylabel('Peso sinaptico g (mV)')
title('Come cambia l\'intensita\' della sinapsi nel tempo (facilitazione)')

subplot(3, 1, 3)
plot(spikes.t[spikes.i == 0]/ms, [0]*sum(spikes.i == 0), 'o', color='C0', label='A')
plot(spikes.t[spikes.i == 1]/ms, [1]*sum(spikes.i == 1), 'o', color='C1', label='B')
yticks([0, 1], ['A', 'B'])
xlabel('Tempo (ms)')
title('Raster degli spike')
legend(loc='upper right', fontsize=8)

tight_layout()
savefig('plasticita_sinaptica.png', dpi=150)

print("Simulazione completata. Grafico salvato in plasticita_sinaptica.png")
print(f"Spike di A: {sum(spikes.i == 0)}")
print(f"Spike di B: {sum(spikes.i == 1)}")
print(f"Peso sinaptico iniziale: {syn.g[0]/mV:.1f} mV (osserva come sale nel grafico centrale)")
print("\nOsservazione chiave: i primi spike di A potrebbero non bastare a far")
print("sparare B (sinapsi ancora 'debole'), ma dopo alcuni spike ravvicinati")
print("la sinapsi si e' rafforzata abbastanza da far scattare B.")
