# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
PICCOLA RETE: propagazione del segnale su una catena di neuroni
====================================================================
Finora avevamo solo A -> B. Qui costruiamo una catena di 8 neuroni:
  N0 -> N1 -> N2 -> N3 -> N4 -> N5 -> N6 -> N7

Solo N0 riceve stimolo esterno. Ogni neurone successivo si attiva
SOLO se riceve abbastanza spinta da quello precedente.

Cosa possiamo osservare che con 2 neuroni non si vedeva:
  1. Il ritardo si ACCUMULA lungo la catena (ogni sinapsi aggiunge
     il suo delay)
  2. Se l'intensita' sinaptica e' troppo debole, il segnale puo'
     "morire" a meta' catena (non arriva fino in fondo) — un
     fenomeno reale nei circuiti neurali biologici
  3. Con connessioni piu' complesse (non solo a catena) iniziano
     a comparire fenomeni come sincronizzazione o feedback
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

N = 8   # numero di neuroni nella catena

start_scope()

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV
neurons.u = b * neurons.v

# Solo il primo neurone (indice 0) riceve stimolo esterno.
neurons.I = 0*mV/ms
neurons.I[0] = 20*mV/ms

# ---------------------------------------------------------------
# CONNESSIONI A CATENA: neurone i -> neurone i+1
# ---------------------------------------------------------------
syn = Synapses(neurons, neurons, on_pre='v_post += 22*mV')
# i, i+1 per i che va da 0 a N-2 (crea la catena 0->1->2->...->7)
syn.connect(i=arange(N-1), j=arange(1, N))
syn.delay = 1.5*ms

# ---------------------------------------------------------------
# MONITORAGGIO E SIMULAZIONE
# ---------------------------------------------------------------
spikes = SpikeMonitor(neurons)
run(200*ms)

# ---------------------------------------------------------------
# VISUALIZZAZIONE
# ---------------------------------------------------------------
figure(figsize=(10, 5))
plot(spikes.t/ms, spikes.i, 'o', markersize=8)
yticks(range(N), [f'N{i}' for i in range(N)])
xlabel('Tempo (ms)')
ylabel('Neurone nella catena')
title(f'Propagazione del segnale lungo una catena di {N} neuroni')
grid(alpha=0.3)
tight_layout()
savefig('rete_catena.png', dpi=150)

print(f"Simulazione completata con {N} neuroni in catena.")
print(f"Grafico salvato in rete_catena.png\n")

print(f"{'Neurone':>10} | {'N. spike':>10} | {'Primo spike (ms)':>18}")
print("-" * 45)
for idx in range(N):
    spike_times = spikes.t[spikes.i == idx]/ms
    n_spikes = len(spike_times)
