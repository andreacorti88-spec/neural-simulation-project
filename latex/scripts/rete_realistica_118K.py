"""
Rete E/I con densita' sinaptica realistica (~2600 sinapsi/neurone, livello
Allen Institute) e sinapsi a conduttanza (invece del salto di tensione
v_post += w usato negli script precedenti): ge/gi decadono esponenzialmente
e spingono v verso un potenziale di inversione (E_exc=0mV tipo AMPA,
E_inh=-80mV tipo GABA-A). Budget di memoria fisso -> aumentare la densita'
21x rispetto a rete_40M_finale costringe a ridurre N a ~118K.

Pesi w_e/w_i ricalibrati empiricamente per K=1300: la regola "pesi
inversamente proporzionali a K" non basta, la rete e' molto piu' sensibile
del previsto quando K cresce di un ordine di grandezza (target ~18Hz stabile).

Costo: ~19 min per blocco di 200ms, run completo stimato 45-90 min.
"""

from brian2 import *
import numpy as np
import time as pytime

print("Costruzione rete: 118.000 neuroni, ~2600 sinapsi/neurone...")
print("La generazione della connettivita' e la compilazione richiedono qualche minuto.")

start_scope()

N = 118_000
N_E = int(N * 0.8)
N_I = N - N_E
K = 1300  # connessioni/neurone per tipo -> ~2600 sinapsi/neurone totali

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms
E_exc = 0*mV      # AMPA
E_inh = -80*mV    # GABA-A
tau_e = 5*ms
tau_i = 10*ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + ge*(E_exc-v) + gi*(E_inh-v) + I_ext : volt
du/dt = a*(b*v - u) : volt/second
dge/dt = -ge/tau_e : 1/second
dgi/dt = -gi/tau_i : 1/second
I_ext : volt/second
'''

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV + 10*mV*rand(N)
neurons.u = b * neurons.v[:]
neurons.I_ext = 0*mV/ms

E = neurons[:N_E]
I_pop = neurons[N_E:]

background = PoissonGroup(N, rates=1200*Hz)
bg_syn = Synapses(background, neurons, on_pre='ge_post += 15/second')
bg_syn.connect(j='i')

# pesi calibrati empiricamente per K=1300
w_e = 0.22/second
w_i = 0.9/second


def gen_fixed_degree(n_source, n_target, K, exclude_self=False):
    """Connettivita' a grado fisso, vettorializzata (piu' efficiente della
    sintassi p= di Brian2 quando K e' dell'ordine delle migliaia)."""
    src = np.repeat(np.arange(n_source), K)
    tgt = np.random.randint(0, n_target, size=n_source*K)
    if exclude_self:
        mask = src != tgt
        src, tgt = src[mask], tgt[mask]
    return src, tgt


t0 = pytime.time()
src, tgt = gen_fixed_degree(N_E, N_E, K, exclude_self=True)
syn_ee = Synapses(E, E, on_pre='ge_post += w_e')
syn_ee.connect(i=src, j=tgt)

src, tgt = gen_fixed_degree(N_E, N_I, K)
syn_ei = Synapses(E, I_pop, on_pre='ge_post += w_e')
syn_ei.connect(i=src, j=tgt)

src, tgt = gen_fixed_degree(N_I, N_E, K)
syn_ie = Synapses(I_pop, E, on_pre='gi_post += w_i')
syn_ie.connect(i=src, j=tgt)

src, tgt = gen_fixed_degree(N_I, N_I, K, exclude_self=True)
syn_ii = Synapses(I_pop, I_pop, on_pre='gi_post += w_i')
syn_ii.connect(i=src, j=tgt)

n_syn_total = len(syn_ee) + len(syn_ei) + len(syn_ie) + len(syn_ii) + len(bg_syn)
print(f"Rete costruita in {pytime.time()-t0:.1f}s.")
print(f"Sinapsi totali: {n_syn_total:,} ({n_syn_total/N:.0f} per neurone in media)")

net = Network(collect())

print("\nInizio burn-in (200ms simulati)...")
t0 = pytime.time()
net.run(200*ms)
print(f"Burn-in completato in {pytime.time()-t0:.1f}s.")

rate_mon = PopulationRateMonitor(neurons)
raster_subset = SpikeMonitor(neurons[:1500], record=True)
net.add(rate_mon, raster_subset)

print("Inizio baseline...")
t0 = pytime.time()
net.run(150*ms)
print(f"Baseline completata in {pytime.time()-t0:.1f}s.")

N_STIM = int(N_E * 0.005)
stim_idx = np.arange(N_STIM)
stim_time = 350*ms
neurons.I_ext[stim_idx] = 40*mV/ms
net.run(20*ms)
neurons.I_ext[stim_idx] = 0*mV/ms
print(f"Stimolo applicato a {N_STIM:,} neuroni ({N_STIM/N*100:.3f}% della rete).")

print("Inizio recovery...")
t0 = pytime.time()
net.run(150*ms)
print(f"Recovery completata in {pytime.time()-t0:.1f}s.")

r = np.array(rate_mon.smooth_rate(window='flat', width=5*ms)/Hz)
t = np.array(rate_mon.t/ms)

baseline_rate = r[(t > 200) & (t < 350)].mean()
peak_rate = r[(t >= 350) & (t < 400)].max()
recovery_rate = r[t > 490].mean()

MOUSE_BRAIN_NEURONS = 70_000_000
pct_mouse = 100 * N / MOUSE_BRAIN_NEURONS

print("\n" + "=" * 60)
print("RISULTATI FINALI")
print("=" * 60)
print(f"Scala: {N:,} neuroni, {n_syn_total:,} sinapsi ({n_syn_total/N:.0f}/neurone)")
print(f"Confronto densita': Allen Institute ~2600/neurone, "
      f"cervello di topo reale: migliaia/neurone")
print(f"Percentuale di un cervello di topo (per numero di neuroni): {pct_mouse:.3f}%")
print(f"Frequenza baseline: {baseline_rate:.2f} Hz")
print(f"Picco dopo lo stimolo: {peak_rate:.2f} Hz")
print(f"Frequenza dopo il recovery: {recovery_rate:.2f} Hz")
if abs(recovery_rate - baseline_rate) < 3:
    print(">>> La rete e' tornata al livello di attivita' di partenza "
          "anche con sinapsi a conduttanza e densita' realistica. <<<")

figure(figsize=(11, 7))

subplot(2, 1, 1)
plot(raster_subset.t/ms, raster_subset.i, '.', color='C0', markersize=1.5)
axvspan(stim_time/ms, (stim_time+20*ms)/ms, color='red', alpha=0.15)
ylabel('Indice neurone\n(sottoinsieme di 1500)')
title(f'Rete a {N:,} neuroni, {n_syn_total/N:.0f} sinapsi/neurone '
      f'(densita\' Allen Institute), sinapsi a conduttanza')

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
savefig('rete_realistica_118K.png', dpi=150)
print("\nGrafico salvato in rete_realistica_118K.png")
