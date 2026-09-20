"""
Verifica Monte Carlo di dialogo_spiking_multi.py (sezione 5.10): la singola
run mostrava nessuna interferenza tra due messaggi in dominio spiking, ma
si basava su UN SOLO seed casuale -- stesso limite gia' affrontato nel
modello a rate con dialogo_multi_scoperta_locale_montecarlo.py (sezione
5.5.10). Qui si ripete lo stesso esperimento su piu' seed indipendenti
(rumore di fondo, connettivita' intra- e inter-anello tutti derivati dal
seed passato da riga di comando) per verificare che il risultato positivo
non sia un colpo di fortuna di un singolo run.

Uso: python3 dialogo_spiking_multi_montecarlo.py <SEED>
Stampa una riga "MC_RIGA" parsabile con l'esito per ciascun messaggio.
"""

import sys
from brian2 import *

SEED = int(sys.argv[1])

start_scope()
seed(SEED)

N_E = 100
N_I = 25
N_ring = N_E + N_I

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms
tau_nmda = 100*ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I + g_nmda : volt
du/dt = a*(b*v - u) : volt/second
dg_nmda/dt = -g_nmda/tau_nmda : volt/second
I : volt/second
x : 1
'''

neurons = NeuronGroup(2*N_ring, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV + 5*mV*randn(2*N_ring)
neurons.u = b * neurons.v[:]
neurons.I = 0*mV/ms
neurons.g_nmda = 0*mV/ms

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


def costruisci_anello(E, I_pop, seed_locale):
    dmat = dist_circ(xi[:, None] - xi[None, :])
    prob_ee = np.exp(-dmat**2 / (2*sigma_exc**2))
    prob_ee[dmat > CUTOFF_LOCALE] = 0.0
    np.fill_diagonal(prob_ee, 0)
    rng = np.random.RandomState(seed_locale)
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


syn_A = costruisci_anello(A_E, A_I, seed_locale=1000+SEED)
syn_B = costruisci_anello(B_E, B_I, seed_locale=2000+SEED)

tau_stdp = 20*ms
A_ltp = 0.15*mV
A_ltd = 0.15*mV
w_max = 5.0*mV
w_min = 0*mV
w_init = 1.5*mV
sigma_inter = 0.08
p_inter_picco = 0.5

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
prob_inter = p_inter_picco * np.exp(-dmat_inter**2 / (2*sigma_inter**2))
rng_inter = np.random.RandomState(3000+SEED)
mask_inter = rng_inter.rand(N_E, N_E) < prob_inter
i_inter, j_inter = np.nonzero(mask_inter)
syn_ab.connect(i=i_inter, j=j_inter)
syn_ab.w = w_init

spikes = SpikeMonitor(neurons)

MESSAGGI = [0.15, 0.65]
larghezza_stimolo = 0.05
gruppi_stim = []
for m in MESSAGGI:
    peso_stim = np.exp(-dist_circ(xi - m)**2 / (2*larghezza_stimolo**2))
    gruppi_stim.append(np.where(peso_stim > 0.5)[0])

N_RIPETIZIONI_TOT = 60
T_STIM = 60*ms
T_PAUSA = 340*ms

storico = []
t_corrente = 0*ms
for rep in range(N_RIPETIZIONI_TOT):
    k = rep % len(MESSAGGI)
    storico.append((k, t_corrente/ms))
    A_E.I[gruppi_stim[k]] = 15*mV/ms
    run(T_STIM)
    A_E.I[:] = 0*mV/ms
    run(T_PAUSA)
    t_corrente += T_STIM + T_PAUSA

ciclo_ms = (T_STIM+T_PAUSA)/ms
mask_b = (spikes.i >= N_ring) & (spikes.i < N_ring+N_E)


def valuta_messaggio(k, n_eval):
    tempi = [t0 for (kk, t0) in storico if kk == k]
    primi, ultimi = tempi[:n_eval], tempi[-n_eval:]

    def conta(lista):
        n_tot, posizioni = 0, []
        for t0 in lista:
            fin = t0 + ciclo_ms
            m = mask_b & (spikes.t/ms >= t0) & (spikes.t/ms < fin)
            idx = spikes.i[m] - N_ring
            n_tot += m.sum()
            posizioni.extend(list(xi[idx]))
        centro = None
        if posizioni:
            centro = float(np.angle(np.mean(np.exp(2j*np.pi*np.array(posizioni)))) / (2*np.pi) % 1.0)
        return n_tot, centro
    return conta(primi), conta(ultimi)


SOGLIA_SUCCESSO = 0.10  # errore di localizzazione entro 2x la larghezza dello stimolo
risultati = []
for k, m in enumerate(MESSAGGI):
    (n_i, c_i), (n_f, c_f) = valuta_messaggio(k, 5)
    err_f = dist_circ(c_f - m) if c_f is not None else 1.0
    cresciuto = n_f > n_i
    successo = (err_f < SOGLIA_SUCCESSO) and cresciuto and n_f > 0
    risultati.append((m, n_i, n_f, err_f, successo))

esiti = ''.join('1' if r[4] else '0' for r in risultati)
print(f"MC_RIGA seed={SEED} esiti={esiti} " +
      " ".join(f"m{k}(x={m}):n_i={n_i},n_f={n_f},err={err_f:.3f},ok={int(ok)}"
               for k, (m, n_i, n_f, err_f, ok) in enumerate(risultati)))
