"""
Test diagnostico rapido: l'esplosione vista in ogni tentativo di
bidirezionalita' dalla sezione 5.9 in poi e' davvero causata dalla
bidirezionalita' (due anelli mutuamente accoppiati), o semplicemente dal
fatto che w_init/w_max sono stati alzati (3.5mV/8.0mV, dal tentativo 3 in
poi) rispetto ai valori originali stabili di dialogo_spiking.py (sezione
5.7: w_init=1.5mV, w_max=5.0mV) e mai piu' riportati indietro nei
tentativi successivi?

Qui si ripete l'identica configurazione UNIDIREZIONALE di dialogo_
spiking.py (solo A->B, nessun canale di ritorno, esattamente come nella
sezione 5.7 gia' verificata stabile), ma con i pesi "forti" (w_init=3.5mV,
w_max=8.0mV) e la stessa cadenza piu' lenta usata nei tentativi di
bidirezionalita' (T_STIM=60ms, T_PAUSA=940ms, 60 ripetizioni, 60 secondi
totali). Se ANCHE questo canale puramente feedforward esplode, la causa
vera non e' la retroazione tra i due anelli ma semplicemente l'aver
alzato il tetto dei pesi senza piu' riabbassarlo.
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

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I + g_nmda - adapt : volt
du/dt = a*(b*v - u) : volt/second
dg_nmda/dt = -g_nmda/tau_nmda : volt/second
dadapt/dt = -adapt/tau_adapt_pop : volt/second
I : volt/second
x : 1
'''
tau_adapt_pop = 200*ms
gain_adapt_pop = 0.15*mV/ms

neurons = NeuronGroup(2*N_ring, eqs, threshold='v > 30*mV',
                       reset='v = c; u += d; adapt += gain_adapt_pop',
                       method='euler')
neurons.v = -65*mV + 5*mV*randn(2*N_ring)
neurons.u = b * neurons.v[:]
neurons.I = 0*mV/ms
neurons.g_nmda = 0*mV/ms
neurons.adapt = 0*mV/ms

A_E = neurons[:N_E]
A_I = neurons[N_E:N_ring]
B_E = neurons[N_ring:N_ring+N_E]
B_I = neurons[N_ring+N_E:]

xi = np.arange(N_E) / N_E
A_E.x = xi
B_E.x = xi

background = PoissonGroup(2*N_ring, rates=400*Hz)
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


def costruisci_anello(E, I_pop, seed):
    dmat = dist_circ(xi[:, None] - xi[None, :])
    prob_ee = np.exp(-dmat**2 / (2*sigma_exc**2))
    prob_ee[dmat > CUTOFF_LOCALE] = 0.0
    np.fill_diagonal(prob_ee, 0)
    rng = np.random.RandomState(seed)
    mask_ee = rng.rand(N_E, N_E) < prob_ee
    i_idx, j_idx = np.nonzero(mask_ee)
    syn_ee = Synapses(E, E, on_pre='g_nmda_post += w_ee_nmda')
    syn_ee.connect(i=i_idx, j=j_idx)
    syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
    syn_ei.connect(p=p_conn)
    syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
    syn_ie.connect(p=p_conn)
    syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
    syn_ii.connect(condition='i!=j', p=p_conn)
    return syn_ee, syn_ei, syn_ie, syn_ii


syn_ee_A, syn_ei_A, syn_ie_A, syn_ii_A = costruisci_anello(A_E, A_I, seed=0)
syn_ee_B, syn_ei_B, syn_ie_B, syn_ii_B = costruisci_anello(B_E, B_I, seed=1)

tau_stdp = 20*ms
A_ltp = 0.15*mV
A_ltd = 0.15*mV
w_max = 8.0*mV   # "forte", come nei tentativi bidirezionali con interruttore
w_min = 0*mV
w_init = 3.5*mV  # "forte", idem

stdp_model = '''
w : volt
dapre/dt = -apre/tau_stdp : volt (event-driven)
dapost/dt = -apost/tau_stdp : volt (event-driven)
'''
stdp_on_pre = '''
v_post += w
apre += A_ltp
w = clip(w+apost, w_min, w_max)
'''
stdp_on_post = '''
apost -= A_ltd
w = clip(w+apre, w_min, w_max)
'''

syn_ab = Synapses(A_E, B_E, model=stdp_model, on_pre=stdp_on_pre, on_post=stdp_on_post)
dmat_inter = dist_circ(xi[:, None] - xi[None, :])
sigma_inter = 0.08
p_inter_picco = 0.5
prob_inter = p_inter_picco * np.exp(-dmat_inter**2 / (2*sigma_inter**2))
rng_inter = np.random.RandomState(2)
mask_inter = rng_inter.rand(N_E, N_E) < prob_inter
i_inter, j_inter = np.nonzero(mask_inter)
syn_ab.connect(i=i_inter, j=j_inter)
syn_ab.w = w_init
print(f"Sinapsi inter-anello A->B (SOLO questo verso, nessun B->A): {len(syn_ab)} "
      f"(peso iniziale {w_init}, tetto {w_max})")

spikes = SpikeMonitor(neurons)

MESSAGGIO = 0.3
larghezza_stimolo = 0.05
peso_stim = np.exp(-dist_circ(xi - MESSAGGIO)**2 / (2*larghezza_stimolo**2))
gruppo_stim = np.where(peso_stim > 0.5)[0]

N_RIPETIZIONI = 60
T_STIM = 60*ms
T_PAUSA = 940*ms   # stessa cadenza lenta dei tentativi di bidirezionalita'

print(f"\nTraining: {N_RIPETIZIONI} ripetizioni, SOLO canale A->B, "
      f"stessa cadenza dei test bidirezionali ({(T_STIM+T_PAUSA)*N_RIPETIZIONI/ms:.0f}ms totali)")

for rep in range(N_RIPETIZIONI):
    A_E.I[gruppo_stim] = 15*mV/ms
    run(T_STIM)
    A_E.I[:] = 0*mV/ms
    run(T_PAUSA)
    if rep % 10 == 0 or rep == N_RIPETIZIONI-1:
        n_a = int((spikes.i < N_E).sum())
        n_b = int(((spikes.i >= N_ring) & (spikes.i < N_ring+N_E)).sum())
        print(f"  ripetizione {rep:3d}/{N_RIPETIZIONI}: "
              f"A->B peso medio={np.mean(syn_ab.w)/mV:.3f}mV max={np.max(syn_ab.w)/mV:.3f}mV  |  "
              f"spike finora A={n_a} B={n_b}")

print("\nTraining completato.")
n_a_tot = int((spikes.i < N_E).sum())
n_b_tot = int(((spikes.i >= N_ring) & (spikes.i < N_ring+N_E)).sum())
print(f"\nSpike totali: A={n_a_tot}, B={n_b_tot} "
      f"(atteso migliaia se stabile come in 5.7; milioni se esplosione)")
if n_a_tot > 200000 or n_b_tot > 200000:
    print(">>> ESPLOSIONE anche in un canale PURAMENTE UNIDIREZIONALE: "
          "la causa e' il tetto dei pesi alzato (w_max=8mV), NON la "
          "bidirezionalita' in se'. <<<")
else:
    print(">>> Stabile: il canale unidirezionale con pesi forti regge. "
          "L'esplosione nei test bidirezionali e' quindi davvero legata "
          "alla retroazione tra i due anelli, non solo al tetto dei pesi. <<<")
