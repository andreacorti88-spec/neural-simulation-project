"""
Convergenza 3:1 (N0,N1,N2 -> N3): confronto input sincroni vs sfasati,
per isolare l'effetto di sommazione temporale sulla risposta di N3.
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

def run_convergent(offsets_ms, synapse_strength_mV, label):
    """N0,N1,N2 sparano una volta ciascuno ai tempi offsets_ms, convergono su N3."""
    start_scope()

    neurons = NeuronGroup(4, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                           method='euler')
    neurons.v = -65*mV
    neurons.u = b * neurons.v
    neurons.I = 0*mV/ms

    # SpikeGeneratorGroup per fissare i tempi esatti di N0,N1,N2
    driver = SpikeGeneratorGroup(3, indices=[0, 1, 2],
                                   times=[t*ms for t in offsets_ms])

    drive_syn = Synapses(driver, neurons, on_pre='v_post += 35*mV')
    drive_syn.connect(i=[0, 1, 2], j=[0, 1, 2])  # driver -> N0,N1,N2

    converge_syn = Synapses(neurons, neurons,
                              on_pre=f'v_post += {synapse_strength_mV}*mV')
    converge_syn.connect(i=[0, 1, 2], j=[3, 3, 3])  # N0,N1,N2 -> N3
    converge_syn.delay = 1.5*ms

    spikes = SpikeMonitor(neurons)
    mon = StateMonitor(neurons, 'v', record=3)

    run(50*ms)

    n_target_spikes = sum(spikes.i == 3)
    print(f"{label}: N3 ha sparato {n_target_spikes} volte "
          f"(input a t={offsets_ms} ms)")

    return mon.t/ms, mon.v[0]/mV, spikes.t/ms, spikes.i


t_sync, v_sync, spk_t_sync, spk_i_sync = run_convergent(
    [10, 10, 10], 10, "SINCRONO")

t_async, v_async, spk_t_async, spk_i_async = run_convergent(
    [10, 15, 20], 10, "SFASATO")

figure(figsize=(10, 7))

subplot(2, 2, 1)
plot(t_sync, v_sync, color='C0')
axhline(30, color='red', linestyle='--', linewidth=0.8, label='soglia')
title('N3: input SINCRONI (t=10,10,10 ms)')
ylabel('Potenziale (mV)')
legend(fontsize=8)

subplot(2, 2, 2)
plot(t_async, v_async, color='C1')
axhline(30, color='red', linestyle='--', linewidth=0.8, label='soglia')
title('N3: input SFASATI (t=10,15,20 ms)')
legend(fontsize=8)

subplot(2, 2, 3)
for idx in range(4):
    mask = spk_i_sync == idx
    plot(spk_t_sync[mask], spk_i_sync[mask], 'o', markersize=8)
yticks(range(4), ['N0', 'N1', 'N2', 'N3 (target)'])
xlabel('Tempo (ms)')
title('Raster: caso sincrono')

subplot(2, 2, 4)
for idx in range(4):
    mask = spk_i_async == idx
    plot(spk_t_async[mask], spk_i_async[mask], 'o', markersize=8)
yticks(range(4), ['N0', 'N1', 'N2', 'N3 (target)'])
xlabel('Tempo (ms)')
title('Raster: caso sfasato')

suptitle('Sommazione di input: la sincronia conta')
tight_layout()
savefig('rete_convergente.png', dpi=150)
print("\nGrafico salvato: rete_convergente.png")
