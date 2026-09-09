"""
Rete E/I a 300 neuroni con protocollo burn-in / baseline / stimolo / recovery.
Burn-in scarta il transitorio sincronizzato iniziale; poi impulso su 20
eccitatori e verifica se la frequenza di popolazione torna al baseline.
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
neurons.v = -65*mV + 15*mV*rand(N) - 5*mV   # dispersione iniziale piu' ampia
neurons.u = b * neurons.v[:] + 2*mV/ms*randn(N)
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

# burn-in: 200ms senza monitor, scarta il transitorio iniziale
run(200*ms)
print("Burn-in completato (200ms, non registrato).")

# baseline: monitor collegati, 200ms di attivita' spontanea pre-stimolo
rate_mon = PopulationRateMonitor(neurons)
spikes = SpikeMonitor(neurons)

run(200*ms)
print("Baseline registrata (200ms).")

# stimolo: impulso di 20ms su 20 eccitatori (su 300)
N_STIM = 20
stim_group = E[:N_STIM]
stim_group.I = 40*mV/ms
stim_time = 400*ms   # 200 burn-in + 200 baseline

run(20*ms)
stim_group.I = 0*mV/ms
print(f"Stimolo applicato a {N_STIM} neuroni eccitatori (20ms di impulso).")

# recovery
run(280*ms)
print("Fase di recovery completata.")

r = array(rate_mon.smooth_rate(window='flat', width=5*ms)/Hz)
t = array(rate_mon.t/ms)

baseline_rate = r[(t > 200) & (t < 400)].mean()
peak_rate = r[(t >= 400) & (t < 450)].max()
recovery_rate = r[t > 650].mean()

print("\n" + "=" * 55)
print("RISULTATI")
print("=" * 55)
print(f"Frequenza baseline (prima dello stimolo): {baseline_rate:.2f} Hz")
print(f"Picco di frequenza dopo lo stimolo:        {peak_rate:.2f} Hz")
print(f"Frequenza finale (dopo il recovery):       {recovery_rate:.2f} Hz")
if abs(recovery_rate - baseline_rate) < 3:
    print(">>> La rete e' tornata al livello di attivita' di partenza. <<<")

figure(figsize=(11, 8))

subplot(2, 1, 1)
exc_mask = spikes.i < N_E
stim_mask = spikes.i < N_STIM   # sottoinsieme stimolato di E

plot(spikes.t[exc_mask & ~stim_mask]/ms, spikes.i[exc_mask & ~stim_mask], '.',
     color='C0', markersize=2, label=f'Eccitatori (n={N_E-N_STIM})')
plot(spikes.t[~exc_mask]/ms, spikes.i[~exc_mask], '.',
     color='C3', markersize=2, label=f'Inibitori (n={N_I})')
plot(spikes.t[stim_mask]/ms, spikes.i[stim_mask], '.',
     color='gold', markersize=4, label=f'Stimolati (n={N_STIM})')

axvspan(stim_time/ms, (stim_time+20*ms)/ms, color='red', alpha=0.15)
ylabel('Indice neurone')
title('Raster: attivita\' spontanea, poi stimolo (fascia rossa), poi recovery')
legend(loc='upper right', markerscale=4, fontsize=8)

subplot(2, 1, 2)
plot(t, r, color='black', linewidth=1)
axvspan(stim_time/ms, (stim_time+20*ms)/ms, color='red', alpha=0.15, label='stimolo')
axhline(baseline_rate, color='gray', linestyle='--', linewidth=1, label='baseline')
xlabel('Tempo (ms)')
ylabel('Frequenza popolazione (Hz)')
title('Risposta collettiva della rete allo stimolo')
legend(loc='upper right', fontsize=8)
grid(alpha=0.3)

tight_layout()
savefig('rete_300_stimolo.png', dpi=150)
print("\nGrafico salvato in rete_300_stimolo.png")
