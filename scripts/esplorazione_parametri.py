"""
Sweep su delay sinaptico e intensita' (A -> B): effetto sul timing di B
e curva soglia (intensita' minima perche' B risponda).
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
    """Ritorna tempi/indici degli spike (A e B) per un dato delay/intensita'."""
    start_scope()

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


# sweep 1: delay sinaptico, intensita' fissa a 20 mV
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

# sweep 2: intensita' sinaptica, delay fisso a 1.5 ms
strengths_to_test = [2, 5, 10, 20, 30]

print("\nEsperimento intensita' sinaptica:")
print(f"{'Intensita (mV)':>15} | {'Spike di A':>10} | {'Spike di B':>10}")
print("-" * 42)

results = []
for strength in strengths_to_test:
    t, i, n_a, n_b = run_simulation(1.5, strength)
    results.append((strength, n_b))
    print(f"{strength:>15} | {n_a:>10} | {n_b:>10}")

figure(figsize=(7, 4))
strengths, n_b_values = zip(*results)
plot(strengths, n_b_values, 'o-', markersize=8)
xlabel('Intensita\' sinaptica (mV)')
ylabel('Numero di spike di B')
title('Soglia sinaptica: sotto un certo valore, B non risponde piu\'')
axhline(0, color='gray', linewidth=0.5)
grid(alpha=0.3)
savefig('esperimento_soglia.png', dpi=150)
print("\nGrafico 2 salvato: esperimento_soglia.png")
