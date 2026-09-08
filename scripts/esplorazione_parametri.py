# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
ESPLORAZIONE PARAMETRI: delay sinaptico e intensita' della sinapsi
====================================================================
Eseguiamo la stessa simulazione (A -> B) piu' volte, cambiando:
  - il delay sinaptico (quanto tempo impiega il segnale a "viaggiare")
  - l'intensita' della sinapsi (quanto forte e' lo stimolo che B riceve)

Cosi' vediamo concretamente:
  1. Come il delay sposta nel tempo la risposta di B
  2. Come, sotto una certa intensita', B smette del tutto di sparare
     (concetto di SOGLIA SINAPTICA: se lo stimolo non basta a far
       superare la soglia di attivazione, il neurone non risponde)
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

def run_simulation(delay_ms, synapse_strength_mV):
    """Esegue la simulazione con un dato delay e una data intensita' sinaptica.
    Ritorna i tempi e gli indici degli spike (per A e B)."""
    start_scope() # resetta l'ambiente Brian2 tra una run e l'altra

      neurons = NeuronGroup(2, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                             method='euler')
      neurons.v = -65*mV
      neurons.u = b * neurons.v
      neurons.I = [15*mV/ms, 0*mV/ms]

      syn = Synapses(neurons, neurons, on_pre=f'v_post += {synapse_strength_mV}*mV')
      syn.connect(i=0, j=1)
      syn.delay = delay_ms*ms

      spikes = SpikeMonitor(neurons)
      run(200*ms)

      return spikes.t/ms, spikes.i, sum(spikes.i == 0), sum(spikes.i == 1)



# ---------------------------------------------------------------
# ESPERIMENTO 1: variare il delay sinaptico (intensita' fissa a 20 mV)
# ---------------------------------------------------------------
delays_to_test = [0.5, 1.5, 5.0, 10.0]

figure(figsize=(10, 6))
for idx, delay in enumerate(delays_to_test):
    t, i, n_a, n_b = run_simulation(delay, 20)
    subplot(len(delays_to_test), 1, idx+1)
    plot(t[i == 0], i[i == 0], 'o', color='C0', label='A' if idx == 0 else None)
    plot(t[i == 1], i[i == 1]+0.1, 'o', color='C1', label='B' if idx == 0 else None)
    yticks([0, 1.1], ['A', 'B'])
    ylabel(f'delay={delay}ms', fontsize=9)
    if idx == 0:
        legend(loc='upper right', fontsize=8)
xlabel('Tempo (ms)')
suptitle('Effetto del delay sinaptico sulla risposta di B')
tight_layout()
savefig('esperimento_delay.png', dpi=150)
print("Grafico 1 salvato: esperimento_delay.png")

# ---------------------------------------------------------------
