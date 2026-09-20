"""
Prima di accettare la sezione 5.9 come conclusione definitiva, un controllo
mai fatto finora: il canale B->A produce ALMENO un effetto sotto-soglia su
A (depolarizzazione misurabile del potenziale di membrana, anche senza mai
arrivare a un vero spike)? Se anche questo e' assente, la conclusione
negativa e' ancora piu' solida (il canale non fa proprio nulla). Se invece
c'e' un effetto sotto-soglia, vuol dire che il segnale esiste ma non basta
a superare la soglia -- un problema di intensita' assoluta, non di
principio, che potrebbe essere risolto rinforzando la convergenza
sinaptica invece del solo peso per sinapsi.

Setup identico a dialogo_spiking_bidirezionale_gate.py (turni A/B alternati,
plasticita' B->A gated), ma con un StateMonitor sul potenziale di membrana
di un gruppo di neuroni A_E vicino alla posizione del messaggio, campionato
durante i turni 'B' (nessuno stimolo diretto su A in quei turni).
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


syn_A = costruisci_anello(A_E, A_I, seed=0)
syn_B = costruisci_anello(B_E, B_I, seed=1)

tau_stdp = 20*ms
A_ltp = 0.15*mV
A_ltd = 0.15*mV
w_max = 8.0*mV
w_min = 0*mV
w_init = 3.5*mV
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
rng_inter = np.random.RandomState(2)
mask_inter = rng_inter.rand(N_E, N_E) < prob_inter
i_inter, j_inter = np.nonzero(mask_inter)
syn_ab.connect(i=i_inter, j=j_inter)
syn_ab.w = w_init

stdp_model_ba = '''
w : volt
dapre/dt = -apre/tau_stdp : volt (event-driven)
dapost/dt = -apost/tau_stdp : volt (event-driven)
ltp_gate : volt (shared)
ltd_gate : volt (shared)
'''
stdp_on_pre_ba = '''
v_post += w
apre += ltp_gate
w = clip(w+apost, w_min, w_max)
'''
stdp_on_post_ba = '''
apost -= ltd_gate
w = clip(w+apre, w_min, w_max)
'''

syn_ba = Synapses(B_E, A_E, model=stdp_model_ba, on_pre=stdp_on_pre_ba, on_post=stdp_on_post_ba)
mask_inter_ba = rng_inter.rand(N_E, N_E) < prob_inter
i_inter_ba, j_inter_ba = np.nonzero(mask_inter_ba)
syn_ba.connect(i=i_inter_ba, j=j_inter_ba)
syn_ba.w = w_init
syn_ba.ltp_gate = A_ltp
syn_ba.ltd_gate = A_ltd

MESSAGGIO = 0.3
larghezza_stimolo = 0.05
peso_stim = np.exp(-dist_circ(xi - MESSAGGIO)**2 / (2*larghezza_stimolo**2))
gruppo_stim = np.where(peso_stim > 0.5)[0]

# ATTENZIONE: gruppo_stim (i neuroni con peso_stim>0.5) riceve stimolo
# ESTERNO diretto durante i turni 'A' -- includerlo nel gruppo "vicino"
# confonderebbe il loro adattamento residuo (post-scarica, dal proprio
# stimolo diretto) con un eventuale segnale del canale B->A. Il gruppo
# "vicino" deve quindi essere vicino al messaggio (dentro il raggio di
# convergenza sigma_inter=0.08 delle sinapsi B->A) ma ESCLUDERE
# esplicitamente gruppo_stim.
dentro_convergenza = dist_circ(xi - MESSAGGIO) < sigma_inter
gruppo_vicino = np.where(dentro_convergenza & ~np.isin(np.arange(N_E), gruppo_stim))[0]
gruppo_lontano = np.where(dist_circ(xi - (MESSAGGIO + 0.5) % 1.0) < larghezza_stimolo)[0]
print(f"Neuroni A_E stimolati direttamente (esclusi dal confronto): {gruppo_stim}")
print(f"Neuroni A_E vicini al messaggio (mai stimolati direttamente): {gruppo_vicino}")
print(f"Neuroni A_E di controllo (lontani): {gruppo_lontano}")

n_vicino = len(gruppo_vicino)
v_mon = StateMonitor(A_E, 'v', record=np.concatenate([gruppo_vicino, gruppo_lontano]), dt=2*ms)
spikes = SpikeMonitor(neurons)

N_CICLI = 40
T_STIM = 60*ms
T_PAUSA = 940*ms

turni = []
t_corrente = 0*ms
for ciclo in range(N_CICLI):
    turno = 'A' if ciclo % 2 == 0 else 'B'
    turni.append((turno, t_corrente/ms))
    if turno == 'A':
        syn_ba.ltp_gate = 0*mV
        syn_ba.ltd_gate = 0*mV
        A_E.I[gruppo_stim] = 15*mV/ms
        run(T_STIM)
        A_E.I[:] = 0*mV/ms
    else:
        syn_ba.ltp_gate = A_ltp
        syn_ba.ltd_gate = A_ltd
        B_E.I[gruppo_stim] = 15*mV/ms
        run(T_STIM)
        B_E.I[:] = 0*mV/ms
    run(T_PAUSA)
    t_corrente += T_STIM + T_PAUSA
    if ciclo % 10 == 0 or ciclo == N_CICLI-1:
        print(f"  ciclo {ciclo:3d}/{N_CICLI} (turno {turno}): "
              f"B->A peso medio={np.mean(syn_ba.w)/mV:.3f}mV max={np.max(syn_ba.w)/mV:.3f}mV")

print("\nSimulazione completata. Analisi sotto-soglia in corso...")

ciclo_ms = (T_STIM+T_PAUSA)/ms
t_arr = np.asarray(v_mon.t/ms)
v_arr = np.asarray(v_mon.v/mV)  # shape (n_neuroni_monitorati, n_tempi)
n_vicino = len(gruppo_vicino)

turni_b = [t0 for (l, t0) in turni if l == 'B']
primi_b = turni_b[:5]
ultimi_b = turni_b[-5:]


def v_medio_finestra(lista_turni, indici_neuroni, offset_da_t0, durata):
    valori = []
    for t0 in lista_turni:
        mask_t = (t_arr >= t0+offset_da_t0) & (t_arr < t0+offset_da_t0+durata)
        if mask_t.sum() == 0:
            continue
        valori.append(v_arr[indici_neuroni][:, mask_t].mean())
    return np.mean(valori) if valori else None


idx_vicino = np.arange(n_vicino)
idx_lontano = np.arange(n_vicino, v_arr.shape[0])

for nome, lista in [('primi 5 turni B', primi_b), ('ultimi 5 turni B', ultimi_b)]:
    v_vicino_durante = v_medio_finestra(lista, idx_vicino, 20, 100)
    v_lontano_durante = v_medio_finestra(lista, idx_lontano, 20, 100)
    print(f"\n{nome} (finestra 20-120ms dopo inizio turno, esclude impulso diretto):")
    print(f"  v medio A_E VICINO al messaggio (potenziale convergenza B->A): {v_vicino_durante:.3f} mV")
    print(f"  v medio A_E LONTANO (controllo, nessuna convergenza sistematica): {v_lontano_durante:.3f} mV")
    print(f"  differenza (vicino - lontano): {v_vicino_durante - v_lontano_durante:.4f} mV")

fig = figure(figsize=(11, 5))
for t0, colore, label in [(primi_b[0], 'C0', 'primo turno B'), (ultimi_b[-1], 'C3', 'ultimo turno B')]:
    mask_t = (t_arr >= t0) & (t_arr < t0+300)
    v_vicino_traccia = v_arr[idx_vicino][:, mask_t].mean(axis=0)
    plot(t_arr[mask_t]-t0, v_vicino_traccia, color=colore, label=f'{label}, v medio vicino')
xlabel('Tempo dall\'inizio del turno B (ms)')
ylabel('Potenziale di membrana medio (mV), neuroni A vicino al messaggio')
title('Traccia sotto-soglia di A durante i turni B: primo vs ultimo turno')
axvline(60, color='gray', linestyle='--', linewidth=0.8, label='fine impulso diretto su B')
legend(fontsize=8)
tight_layout()
savefig('diagnosi_sottosoglia_ba.png', dpi=150)
print("\nGrafico salvato in diagnosi_sottosoglia_ba.png")
