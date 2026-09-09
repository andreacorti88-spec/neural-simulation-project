"""
STDP pair-based, traccia esponenziale pre/post (Song & Abbott 2001).
Pattern: 5 ripetizioni di A a t, B a t+5ms (A causa B), verifica che w salga
monotonicamente. Da invertire manualmente (B prima di A) per il caso LTD.
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

# parametri STDP
tau_stdp = 20*ms
A_ltp = 1.0*mV
A_ltd = 1.0*mV
w_max = 40*mV
w_min = 0*mV

start_scope()

neurons = NeuronGroup(2, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV
neurons.u = b * neurons.v
neurons.I = 0*mV/ms

# SpikeGeneratorGroup per fissare i tempi di spike di A/B, indipendenti
# dalla dinamica del neurone -- 5 ripetizioni "A a t, B a t+5ms"
pattern_causale = []
pattern_indices = []
for rep in range(5):
    t0 = 20 + rep*30
    pattern_causale += [t0, t0+5]
    pattern_indices += [0, 1]   # 0 = A, 1 = B

driver = SpikeGeneratorGroup(2, indices=pattern_indices,
                               times=[t*ms for t in pattern_causale])

# forza gli spike di A/B direttamente (impulso fortissimo) per tempi certi
force_syn = Synapses(driver, neurons, on_pre='v_post += 50*mV')
force_syn.connect(j='i')

# sinapsi plastica A->B
stdp_eqs = '''
w : volt
dapre/dt = -apre / tau_stdp : volt (event-driven)
dapost/dt = -apost / tau_stdp : volt (event-driven)
'''

syn = Synapses(neurons, neurons, model=stdp_eqs,
                on_pre='''
                apre += A_ltp
                w = clip(w + apost, w_min, w_max)
                ''',
                on_post='''
                apost -= A_ltd
                w = clip(w + apre, w_min, w_max)
                ''',
                method='euler')
syn.connect(i=0, j=1)
syn.w = 5*mV
initial_w = syn.w[0]

w_mon = StateMonitor(syn, 'w', record=0)
spikes = SpikeMonitor(neurons)

run(200*ms)

figure(figsize=(10, 6))

subplot(2, 1, 1)
plot(w_mon.t/ms, w_mon.w[0]/mV, color='purple', linewidth=2)
ylabel('Peso sinaptico w (mV)')
title('STDP: il peso cresce ad ogni ripetizione "A prima di B" (causale)')
grid(alpha=0.3)

subplot(2, 1, 2)
plot(spikes.t[spikes.i == 0]/ms, [0]*sum(spikes.i == 0), 'o', color='C0', label='A')
plot(spikes.t[spikes.i == 1]/ms, [1]*sum(spikes.i == 1), 'o', color='C1', label='B')
yticks([0, 1], ['A', 'B'])
xlabel('Tempo (ms)')
title('Pattern di stimolazione: A spara sempre 5ms prima di B')
legend(loc='upper right', fontsize=8)

tight_layout()
savefig('stdp.png', dpi=150)

print("Simulazione completata. Grafico salvato in stdp.png")
print(f"Peso sinaptico iniziale: {initial_w/mV:.2f} mV")
print(f"Peso sinaptico finale:   {w_mon.w[0][-1]/mV:.2f} mV")
print(f"Rafforzamento totale dopo 5 ripetizioni 'A causa B': "
      f"+{(w_mon.w[0][-1] - initial_w)/mV:.2f} mV")
print("\nCaso LTD (B prima di A) da testare invertendo l'ordine nel pattern.")
