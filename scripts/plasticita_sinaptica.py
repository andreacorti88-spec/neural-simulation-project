"""
Facilitazione sinaptica a breve termine (STF): g si incrementa ad ogni
spike di A e decade con tau_facilitation. A -> B, stimolo a raffica su A.
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

# corrente su A piu' alta del solito -> raffica di spike, effetto facilitazione visibile
neurons.I = [25*mV/ms, 0*mV/ms]

# g: peso corrente della sinapsi, rilassa verso g_base con costante tau_facilitation
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
syn.g = 8*mV                      # sotto soglia da solo
syn.g_base = 8*mV
syn.tau_facilitation = 50*ms
syn.facilitation_increment = 6*mV

mon = StateMonitor(neurons, 'v', record=True)
g_mon = StateMonitor(syn, 'g', record=0)
spikes = SpikeMonitor(neurons)

run(200*ms)

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
print(f"Peso sinaptico iniziale: {syn.g[0]/mV:.1f} mV (vedi andamento nel grafico centrale)")
print("\nI primi spike di A da soli non bastano a far sparare B; dopo alcuni")
print("spike ravvicinati la sinapsi facilitata fa scattare B.")
