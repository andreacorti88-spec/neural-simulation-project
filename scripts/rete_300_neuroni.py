"""
Rete E/I bilanciata, 300 neuroni, stile Brunel: 80/20 E/I, connettivita'
sparsa casuale p=0.1, w_inh > w_exc per tenere il regime asincrono-irregolare
senza saturare o spegnersi. Rumore di fondo Poisson su tutti i neuroni.
"""

from brian2 import *

start_scope()

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
neurons.v = -65*mV + 5*mV*randn(N)   # rompe simmetrie iniziali
neurons.u = b * neurons.v[:]
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

spikes = SpikeMonitor(neurons)
rate_mon = PopulationRateMonitor(neurons)

run(500*ms)

n_spikes_total = sum(spikes.count)
rate_medio = n_spikes_total / N / 0.5  # Hz su 500ms

spikes_E = sum(spikes.count[:N_E])
spikes_I = sum(spikes.count[N_E:])
rate_E = spikes_E / N_E / 0.5
rate_I = spikes_I / N_I / 0.5

print("=" * 55)
print(f"Rete: {N_E} eccitatori + {N_I} inibitori = {N} neuroni totali")
print(f"Connettivita': {p_conn*100:.0f}% (sparsa e casuale)")
print("=" * 55)
print(f"Spike totali in 500ms: {n_spikes_total}")
print(f"Frequenza media della rete: {rate_medio:.2f} Hz")
print(f"Frequenza media eccitatori: {rate_E:.2f} Hz")
print(f"Frequenza media inibitori:  {rate_I:.2f} Hz")

figure(figsize=(11, 7))

subplot(2, 1, 1)
exc_mask = spikes.i < N_E
plot(spikes.t[exc_mask]/ms, spikes.i[exc_mask], '.', color='C0',
     markersize=2, label=f'Eccitatori (n={N_E})')
plot(spikes.t[~exc_mask]/ms, spikes.i[~exc_mask], '.', color='C3',
     markersize=2, label=f'Inibitori (n={N_I})')
ylabel('Indice neurone')
title(f'Raster della rete ({N} neuroni, connettivita\' sparsa al {p_conn*100:.0f}%)')
legend(loc='upper right', markerscale=5, fontsize=8)

subplot(2, 1, 2)
plot(rate_mon.t/ms, rate_mon.smooth_rate(window='flat', width=5*ms)/Hz,
     color='black')
xlabel('Tempo (ms)')
ylabel('Frequenza di scarica\ndella popolazione (Hz)')
title('Attivita\' collettiva della rete nel tempo')
grid(alpha=0.3)

tight_layout()
savefig('rete_300_neuroni.png', dpi=150)
print("\nGrafico salvato in rete_300_neuroni.png")
