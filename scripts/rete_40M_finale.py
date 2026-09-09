"""
Rete E/I bilanciata a 40M neuroni (57.1% scala topo), massima dimensione
verificata praticabile su questa macchina -- stesso protocollo
baseline/stimolo/recovery usato a 300 neuroni.

Memoria: niente spike-per-spike su tutti i 40M (PopulationRateMonitor per
l'aggregato + raster su un sottoinsieme di 2000). Stimolo sullo 0.5% degli
eccitatori (~160k su 32M), sufficiente per un effetto misurabile a livello
di popolazione. Costruzione rete + simulazione: ~20-30 min stimati.
"""

from brian2 import *
import time as pytime
import numpy as np

print("Costruzione della rete (40 milioni di neuroni)...")

start_scope()

N = 40_000_000
N_E = int(N * 0.8)
N_I = N - N_E

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
du/dt = a*(b*v - u) : volt/second
I : volt/second
'''

t_build_start = pytime.time()

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV + 15*mV*rand(N) - 5*mV
neurons.u = b * neurons.v[:] + 2*mV/ms*randn(N)
neurons.I = 0*mV/ms

E = neurons[:N_E]
I_pop = neurons[N_E:]

background = PoissonGroup(N, rates=1500*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
bg_syn.connect(j='i')

K_local = 15
w_exc = 1.0*mV
w_inh = 4.0*mV
K_ei = min(K_local, N_I)
K_ii = min(K_local, N_I)

syn_ee = Synapses(E, E, on_pre='v_post += w_exc')
syn_ee.connect(j=f'k for k in sample(N_E, size={K_local})', skip_if_invalid=True)
syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
syn_ei.connect(j=f'k for k in sample(N_I, size={K_ei})', skip_if_invalid=True)
syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
syn_ie.connect(j=f'k for k in sample(N_E, size={K_local})', skip_if_invalid=True)
syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
syn_ii.connect(j=f'k for k in sample(N_I, size={K_ii})', skip_if_invalid=True)

print(f"Rete costruita in {pytime.time()-t_build_start:.1f}s. Inizio burn-in...")

run(200*ms)
print("Burn-in completato.")

rate_mon = PopulationRateMonitor(neurons)
raster_subset = SpikeMonitor(neurons[:2000], record=True)  # sottoinsieme, non tutti i 40M

run(200*ms)
print("Baseline registrata.")

# stimolo: 0.5% degli eccitatori (~160k su 32M)
N_STIM = int(N_E * 0.005)
stim_group = E[:N_STIM]
stim_time = 400*ms
stim_group.I = 40*mV/ms
run(20*ms)
stim_group.I = 0*mV/ms
print(f"Stimolo applicato a {N_STIM:,} neuroni ({N_STIM/N*100:.3f}% della rete totale).")

run(280*ms)
print("Recovery completato.")

r = np.array(rate_mon.smooth_rate(window='flat', width=5*ms)/Hz)
t = np.array(rate_mon.t/ms)

baseline_rate = r[(t > 200) & (t < 400)].mean()
peak_rate = r[(t >= 400) & (t < 450)].max()
recovery_rate = r[t > 650].mean()

MOUSE_BRAIN_NEURONS = 70_000_000
pct_mouse = 100 * N / MOUSE_BRAIN_NEURONS

print("\n" + "=" * 60)
print("RISULTATI FINALI")
print("=" * 60)
print(f"Scala della rete: {N:,} neuroni ({pct_mouse:.1f}% di un cervello di topo)")
print(f"Frequenza baseline: {baseline_rate:.2f} Hz")
print(f"Picco dopo lo stimolo: {peak_rate:.2f} Hz")
print(f"Frequenza dopo il recovery: {recovery_rate:.2f} Hz")
if abs(recovery_rate - baseline_rate) < 3:
    print(">>> La rete e' tornata al livello di attivita' di partenza "
          "anche a questa scala. <<<")

figure(figsize=(11, 7))

subplot(2, 1, 1)
plot(raster_subset.t/ms, raster_subset.i, '.', color='C0', markersize=2)
axvspan(stim_time/ms, (stim_time+20*ms)/ms, color='red', alpha=0.15)
ylabel('Indice neurone\n(sottoinsieme di 2000 su 40M)')
title(f'Raster di un sottoinsieme -- rete a {N:,} neuroni '
      f'({pct_mouse:.1f}% di un cervello di topo)')

subplot(2, 1, 2)
plot(t, r, color='black', linewidth=1)
axvspan(stim_time/ms, (stim_time+20*ms)/ms, color='red', alpha=0.15, label='stimolo')
axhline(baseline_rate, color='gray', linestyle='--', linewidth=1, label='baseline')
xlabel('Tempo (ms)')
ylabel('Frequenza INTERA popolazione (Hz)')
title(f'Risposta collettiva di tutti i {N:,} neuroni allo stimolo')
legend(loc='upper right', fontsize=8)
grid(alpha=0.3)

tight_layout()
savefig('rete_40M_finale.png', dpi=150)
print("\nGrafico salvato in rete_40M_finale.png")
