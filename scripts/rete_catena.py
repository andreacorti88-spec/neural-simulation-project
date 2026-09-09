"""
Catena N0->N1->...->N7, solo N0 stimolato esternamente.
Verifica: accumulo del delay lungo la catena, e condizioni sotto cui
il segnale si estingue prima di raggiungere N7.
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

N = 8

start_scope()

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV
neurons.u = b * neurons.v

neurons.I = 0*mV/ms
neurons.I[0] = 20*mV/ms

# catena i -> i+1
syn = Synapses(neurons, neurons, on_pre='v_post += 22*mV')
syn.connect(i=arange(N-1), j=arange(1, N))
syn.delay = 1.5*ms

spikes = SpikeMonitor(neurons)
run(200*ms)

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
    first = f"{spike_times[0]:.1f}" if n_spikes > 0 else "mai"
    print(f"{'N'+str(idx):>10} | {n_spikes:>10} | {first:>18}")

if len(spikes.t[spikes.i == N-1]) == 0:
    print(f"\nIl segnale NON e' arrivato fino all'ultimo neurone (N{N-1}).")
    print("Da alzare l'intensita' sinaptica (22*mV -> es. 26*mV) se serve arrivi in fondo.")
else:
    delay_totale = spikes.t[spikes.i == N-1][0]/ms - spikes.t[spikes.i == 0][0]/ms
    print(f"\nIl segnale e' arrivato fino a N{N-1}.")
    print(f"Ritardo totale accumulato dal primo all'ultimo neurone: {delay_totale:.1f} ms")
