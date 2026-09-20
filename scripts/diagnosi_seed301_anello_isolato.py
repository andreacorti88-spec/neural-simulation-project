"""
Il seed 301 esplode identicamente sia in dialogo_spiking_bidirezionale_
omeostasi.py sia nella v2, nonostante l'omeostasi agisca SOLO sulle
sinapsi inter-anello. Ipotesi: il problema non e' la bidirezionalita' ne'
l'omeostasi, ma la connettivita' RICORRENTE INTERNA di uno dei due anelli
per questo particolare seed (costruisci_anello usa seed=1000+SEED per A,
seed=2000+SEED per B -- per SEED=301, seed=1301 e seed=2301) -- una
particolare estrazione casuale della connettivita' a soglia gaussiana
che, per puro caso, produce un cluster locale troppo denso, instabile di
suo indipendentemente da qualunque interazione tra anelli.

Qui si isola UN SOLO anello (nessun anello B, nessuna sinapsi inter-
anello, nessuno stimolo esterno oltre al rumore di fondo) con la stessa
identica connettivita' usata per l'anello A del seed 301, e si osserva
se esplode da solo in 60 secondi di solo rumore di fondo.
"""

from brian2 import *

start_scope()

N_E = 100
N_I = 25
N_ring = N_E + N_I

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms
tau_nmda = 100*ms
tau_adapt_pop = 200*ms
gain_adapt_pop = 0.15*mV/ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I + g_nmda - adapt : volt
du/dt = a*(b*v - u) : volt/second
dg_nmda/dt = -g_nmda/tau_nmda : volt/second
dadapt/dt = -adapt/tau_adapt_pop : volt/second
I : volt/second
x : 1
'''

neurons = NeuronGroup(N_ring, eqs, threshold='v > 30*mV',
                       reset='v = c; u += d; adapt += gain_adapt_pop',
                       method='euler')
neurons.v = -65*mV + 5*mV*randn(N_ring)
neurons.u = b * neurons.v[:]
neurons.I = 0*mV/ms
neurons.g_nmda = 0*mV/ms
neurons.adapt = 0*mV/ms

E = neurons[:N_E]
I_pop = neurons[N_E:]
xi = np.arange(N_E) / N_E
E.x = xi

background = PoissonGroup(N_ring, rates=400*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 2.0*mV')
bg_syn.connect(j='i')

sigma_exc = 0.05
CUTOFF_LOCALE = 3 * sigma_exc
w_ee_nmda = 0.325*mV/ms
p_conn = 0.15
w_exc = 1.0*mV
w_inh = 10.0*mV


def dist_circ(dx):
    dx = abs(dx)
    return np.minimum(dx, 1 - dx)


SEED_TEST = 301
seed_anello = 1000 + SEED_TEST  # stesso seed usato per l'anello A quando SEED=301
print(f"Costruzione anello con seed={seed_anello} (identico all'anello A del seed MC 301)")

dmat = dist_circ(xi[:, None] - xi[None, :])
prob_ee = np.exp(-dmat**2 / (2*sigma_exc**2))
prob_ee[dmat > CUTOFF_LOCALE] = 0.0
np.fill_diagonal(prob_ee, 0)
rng = np.random.RandomState(seed_anello)
mask_ee = rng.rand(N_E, N_E) < prob_ee
i_idx, j_idx = np.nonzero(mask_ee)

# controllo diagnostico: distribuzione del grado in-entrata (quanti input
# ricorrenti riceve ciascun neurone E) -- un cluster anomalo dovrebbe
# mostrarsi come una coda anomala rispetto al resto della distribuzione
grado_entrata = np.bincount(j_idx, minlength=N_E)
print(f"Grado di entrata ricorrente per neurone: media={grado_entrata.mean():.1f}, "
      f"std={grado_entrata.std():.1f}, min={grado_entrata.min()}, max={grado_entrata.max()}")
soglia_anomala = grado_entrata.mean() + 3*grado_entrata.std()
anomali = np.where(grado_entrata > soglia_anomala)[0]
print(f"Neuroni con grado > media+3*std ({soglia_anomala:.1f}): {list(anomali)}")

syn_ee = Synapses(E, E, on_pre='g_nmda_post += w_ee_nmda')
syn_ee.connect(i=i_idx, j=j_idx)
syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
syn_ei.connect(p=p_conn)
syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
syn_ie.connect(p=p_conn)
syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
syn_ii.connect(condition='i!=j', p=p_conn)

spikes = SpikeMonitor(neurons)

print("\nSimulazione: 60 secondi di solo rumore di fondo, nessuno stimolo esterno, nessun altro anello")
run(60*second, report='text', report_period=10*second)

n_tot = int((spikes.i < N_E).sum())
print(f"\nSpike totali E in 60s di solo rumore: {n_tot}")
print(f"(atteso: qualche centinaio/migliaio se stabile come nella sezione 5.6; "
      f"milioni se questo singolo anello esplode gia' da solo)")
if n_tot > 200000:
    print(">>> CONFERMATO: l'anello con questa connettivita' (seed 301) esplode "
          "GIA' DA SOLO, con solo rumore di fondo, senza alcuna interazione "
          "bidirezionale. Il problema e' nella connettivita' ricorrente interna "
          "per questo seed, non nella bidirezionalita' ne' nell'omeostasi. <<<")
else:
    print(">>> Il singolo anello isolato resta stabile: il problema del seed 301 "
          "deve quindi dipendere specificamente dall'interazione bidirezionale, "
          "non dalla sola connettivita' interna. <<<")

fig = figure(figsize=(11, 3))
mask_e = spikes.i < N_E
plot(spikes.t/ms, spikes.i, '.', markersize=1.5)
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title(f"Anello isolato, connettivita' del seed 301, solo rumore di fondo (60s)")
tight_layout()
savefig('diagnosi_seed301_anello_isolato.png', dpi=150)
print("\nGrafico salvato in diagnosi_seed301_anello_isolato.png")
