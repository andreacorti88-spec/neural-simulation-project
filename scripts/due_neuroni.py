"""
A -> B, sinapsi eccitatoria singola, Izhikevich RS.
Sanity check di base sulla toolchain prima di passare a reti piu' grandi.
"""

from brian2 import *

# parametri Izhikevich, regime RS
a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

# I in mV/ms, come nei tutorial Brian2 standard (niente C/R espliciti)
eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
du/dt = a*(b*v - u) : volt/second
I : volt/second
'''

neurons = NeuronGroup(2, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV
neurons.u = b * neurons.v

# solo A (idx 0) riceve corrente esterna; B dipende solo dalla sinapsi
neurons.I = [15*mV/ms, 0*mV/ms]

syn = Synapses(neurons, neurons, on_pre='v_post += 20*mV')
syn.connect(i=0, j=1)
syn.delay = 1.5*ms

mon = StateMonitor(neurons, 'v', record=True)
spikes = SpikeMonitor(neurons)

run(200*ms)

figure(figsize=(10, 5))

subplot(2, 1, 1)
plot(mon.t/ms, mon.v[0]/mV, label='Neurone A (stimolato)')
plot(mon.t/ms, mon.v[1]/mV, label='Neurone B (riceve da A)')
ylabel('Potenziale (mV)')
legend(loc='upper right')
title('Potenziale di membrana: comunicazione A -> B')

subplot(2, 1, 2)
plot(spikes.t/ms, spikes.i, 'o', markersize=8)
yticks([0, 1], ['A', 'B'])
xlabel('Tempo (ms)')
ylabel('Neurone')
title('Raster degli spike (i puntini = quando ogni neurone spara)')

tight_layout()
savefig('due_neuroni_output.png', dpi=150)
print("Simulazione completata. Grafico salvato in due_neuroni_output.png")
print(f"Numero di spike del Neurone A: {sum(spikes.i == 0)}")
print(f"Numero di spike del Neurone B: {sum(spikes.i == 1)}")
