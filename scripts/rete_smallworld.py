"""
Connettivita' small-world (Watts-Strogatz 1998) su anello, 1000 neuroni.
p(connessione) decade esponenzialmente con la distanza circolare + probabilita'
costante di shortcut a lunga distanza. Stimolo su un gruppo spazialmente
contiguo -> verifica di propagazione ondulatoria (distanza media dal punto
di stimolo che cresce monotonicamente nel tempo), assente con connettivita' casuale.
"""

from brian2 import *
import numpy as np

start_scope()

N = 1000
N_E = 800
N_I = 200

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
du/dt = a*(b*v - u) : volt/second
I : volt/second
x : 1
'''

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV + 15*mV*rand(N) - 5*mV
neurons.u = b * neurons.v[:] + 2*mV/ms*randn(N)
neurons.I = 0*mV/ms
neurons.x = arange(N)/N   # posizione sull'anello

E = neurons[:N_E]
I_pop = neurons[N_E:]

background = PoissonGroup(N, rates=1500*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
bg_syn.connect(j='i')

w_exc = 1.0*mV
w_inh = 4.0*mV

# connettivita' dipendente dalla distanza
xi = 0.03             # lunghezza caratteristica di localita'
p_local = 0.5         # prob. massima per vicini immediati
p_shortcut = 0.005    # prob. costante shortcut a lunga distanza

# distanza circolare (0..1, con avvolgimento)
dist_expr = '(0.5 - abs(abs(x_pre-x_post) - 0.5))'
prob_expr = f'{p_local}*exp(-({dist_expr})/{xi}) + {p_shortcut}'

syn_ee = Synapses(E, E, on_pre='v_post += w_exc')
syn_ee.connect(condition='i!=j', p=prob_expr)
syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
syn_ei.connect(p=prob_expr)
syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
syn_ie.connect(p=prob_expr)
syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
syn_ii.connect(condition='i!=j', p=prob_expr)

print(f"Sinapsi create -- EE: {len(syn_ee)}, EI: {len(syn_ei)}, "
      f"IE: {len(syn_ie)}, II: {len(syn_ii)}")

run(200*ms)
print("Burn-in completato.")

spikes = SpikeMonitor(neurons)
rate_mon = PopulationRateMonitor(neurons)

run(100*ms)
print("Baseline registrata.")

# stimolo: 20 neuroni contigui nello spazio (vicino a x=0)
stim_group = E[0:20]
stim_time = 300*ms
stim_group.I = 40*mV/ms
run(15*ms)
stim_group.I = 0*mV/ms
print("Stimolo applicato a una regione LOCALE della rete (non sparsa a caso).")

run(150*ms)
print("Simulazione completata.")

# la distanza media degli spike dal punto di stimolo cresce nel tempo?
t_arr = np.array(spikes.t/ms)
i_arr = np.array(spikes.i)

bins = np.arange(300, 350, 2)
bin_centers, mean_dists, spike_counts = [], [], []
for k in range(len(bins)-1):
    mask = (t_arr >= bins[k]) & (t_arr < bins[k+1])
    if mask.sum() > 0:
        x_vals = i_arr[mask] / N
        circ_dist = 0.5 - np.abs(np.abs(x_vals - 0) - 0.5)
        bin_centers.append((bins[k]+bins[k+1])/2)
        mean_dists.append(circ_dist.mean())
        spike_counts.append(mask.sum())

print("\n" + "=" * 55)
print("ANALISI DELL'ONDA DI PROPAGAZIONE")
print("=" * 55)
print(f"Distanza media dal punto di stimolo nei primi 2ms dopo lo stimolo: "
      f"{mean_dists[0]:.3f}")
print(f"Distanza media al picco (massima estensione dell'onda): "
      f"{max(mean_dists):.3f}")
if mean_dists[-1] > mean_dists[0]:
    print(">>> La distanza cresce nel tempo: confermata una vera "
          "propagazione spaziale (onda). <<<")

figure(figsize=(11, 8))

subplot(3, 1, 1)
plot(spikes.t/ms, spikes.i, '.', color='C0', markersize=1.5, alpha=0.6)
axvspan(stim_time/ms, (stim_time+15*ms)/ms, color='red', alpha=0.15)
axhspan(0, 20, color='gold', alpha=0.3)
ylabel('Indice neurone\n(= posizione sull\'anello)')
title('Raster: propagazione spaziale dopo stimolo su una regione locale '
      '(fascia gialla = neuroni stimolati)')

subplot(3, 1, 2)
plot(bin_centers, mean_dists, 'o-', color='C3')
axhline(0.25, color='gray', linestyle='--', linewidth=1,
        label='distanza media attesa se random (nessuna onda)')
xlabel('Tempo (ms)')
ylabel('Distanza media\ndal punto di stimolo')
title('L\'onda si allontana dal punto di stimolo nel tempo')
legend(loc='lower right', fontsize=8)
grid(alpha=0.3)

subplot(3, 1, 3)
r = array(rate_mon.smooth_rate(window='flat', width=3*ms)/Hz)
t = array(rate_mon.t/ms)
plot(t, r, color='black', linewidth=1)
axvspan(stim_time/ms, (stim_time+15*ms)/ms, color='red', alpha=0.15, label='stimolo')
xlabel('Tempo (ms)')
ylabel('Frequenza popolazione (Hz)')
title('Attivita\' collettiva della rete')
legend(loc='upper right', fontsize=8)
grid(alpha=0.3)

tight_layout()
savefig('rete_smallworld.png', dpi=150)
print("\nGrafico salvato in rete_smallworld.png")
