"""
Estensione di ring_attractor_spiking.py (il primo ring attractor spiking
stabile del progetto, sezione 5.6) verso la comunicazione appresa: due
anelli, A e B, collegati da sinapsi STDP VERE (basate su tempi di spike
reali, Song & Abbott 2001 -- stesso schema di stdp.py, non piu'
l'approssimazione Hebbiana a grana di turno usata in tutta la sezione 5.5).

Si comincia da un solo verso (A->B), come il progetto ha sempre fatto nel
dominio a rate (dialogo_stdp.py prima di dialogo_stdp_adattamento.py):
un messaggio ripetuto (stimolo periodico su A, sempre alla stessa
posizione) durante una simulazione continua, sinapsi inter-anello sparse e
deboli all'inizio, plasticita' STDP che dovrebbe rinforzare selettivamente
le sinapsi dai neuroni di A vicino alla posizione stimolata verso i
neuroni di B che, per pura variabilita' del rumore di fondo, capitano ad
accendersi in corrispondenza -- lo stesso principio di "scoperta tramite
rumore poi consolidamento" della sezione 5.5, ma qui il rumore e' quello
vero del sistema (il rumore di fondo Poisson), non un meccanismo aggiunto
apposta.
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
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I + g_nmda : volt
du/dt = a*(b*v - u) : volt/second
dg_nmda/dt = -g_nmda/tau_nmda : volt/second
I : volt/second
x : 1
'''

# due anelli indipendenti, stessa architettura di ring_attractor_spiking.py
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
    # IMPORTANTE: Brian2 scarta silenziosamente gli oggetti Synapses creati
    # dentro una funzione se non restano referenziati da qualche parte --
    # devono essere tutti restituiti e tenuti in una variabile dal chiamante,
    # altrimenti la connettivita' ricorrente semplicemente non esiste nella
    # simulazione (bug scoperto empiricamente nel primo tentativo: 1 solo
    # spike in tutta la corsa, perche' gli anelli non avevano dinamica
    # ricorrente reale)
    return syn_ee, syn_ei, syn_ie, syn_ii


syn_A = costruisci_anello(A_E, A_I, seed=0)
syn_B = costruisci_anello(B_E, B_I, seed=1)

# --- sinapsi inter-anello A->B, STDP vera (Song & Abbott 2001, stesso
# schema di stdp.py) --- connettivita' sparsa, peso iniziale debole: ne'
# innato ne' azzerato, semplicemente troppo poco per accendere un bump in
# B da solo, cosi' la scoperta iniziale dipende dal rumore di fondo di B
tau_stdp = 20*ms
A_ltp = 0.15*mV
A_ltd = 0.15*mV
w_max = 5.0*mV
w_min = 0*mV
w_init = 1.5*mV  # tentativo 1: 0.3mV, B non ha mai sparato (0 spike in
                  # tutta la corsa) -- senza un primo successo, non c'e'
                  # mai uno spike "post" con cui la STDP possa accoppiarsi
sigma_inter = 0.08  # tentativo 2: connettivita' UNIFORME (p_inter=0.3 su
                     # tutta la coppia A-B) faceva rispondere B su tutto
                     # l'anello, diffuso, non localizzato -- niente struttura
                     # spaziale su cui la dinamica propria di B potesse
                     # concentrarsi. Qui la probabilita' di connessione cade
                     # con la distanza circolare (offset zero, come un debole
                     # bias innato "identita'"), un po' piu' larga del
                     # sigma=0.05 intra-anello per lasciare margine di deriva
p_inter_picco = 0.5  # probabilita' di connessione al centro del kernel

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
dmat_inter = dist_circ(xi[:, None] - xi[None, :])  # [i=A, j=B], offset zero
prob_inter = p_inter_picco * np.exp(-dmat_inter**2 / (2*sigma_inter**2))
rng_inter = np.random.RandomState(2)
mask_inter = rng_inter.rand(N_E, N_E) < prob_inter
i_inter, j_inter = np.nonzero(mask_inter)
syn_ab.connect(i=i_inter, j=j_inter)
syn_ab.w = w_init
print(f"Sinapsi inter-anello A->B create: {len(syn_ab)} (peso iniziale {w_init}, "
      f"connettivita' a caduta gaussiana sigma={sigma_inter})")

spikes = SpikeMonitor(neurons)
w_mon = StateMonitor(syn_ab, 'w', record=np.arange(min(200, len(syn_ab))), dt=50*ms)

MESSAGGIO = 0.3
larghezza_stimolo = 0.05
peso_stim = np.exp(-dist_circ(xi - MESSAGGIO)**2 / (2*larghezza_stimolo**2))
gruppo_stim = np.where(peso_stim > 0.5)[0]
print(f"Neuroni di A stimolati (vicino a x={MESSAGGIO}): {len(gruppo_stim)}")

N_RIPETIZIONI = 60
T_STIM = 60*ms
T_PAUSA = 340*ms

print(f"\nTraining: {N_RIPETIZIONI} ripetizioni dello stesso messaggio, "
      f"simulazione continua ({(T_STIM+T_PAUSA)*N_RIPETIZIONI/ms:.0f}ms totali)")

for rep in range(N_RIPETIZIONI):
    A_E.I[gruppo_stim] = 15*mV/ms
    run(T_STIM)
    A_E.I[:] = 0*mV/ms
    run(T_PAUSA)
    if rep % 10 == 0 or rep == N_RIPETIZIONI-1:
        print(f"  ripetizione {rep:3d}/{N_RIPETIZIONI}: peso medio sinapsi A->B = "
              f"{np.mean(syn_ab.w)/mV:.3f} mV, peso massimo = {np.max(syn_ab.w)/mV:.3f} mV")

print("\nTraining completato.")

fig, axes = subplots(3, 1, figsize=(11, 9))

subplot(3, 1, 1)
mask_a = spikes.i < N_E
mask_b = (spikes.i >= N_ring) & (spikes.i < N_ring+N_E)
plot(spikes.t[mask_a]/ms, spikes.i[mask_a], '.', color='C0', markersize=1.5, label='A (E)')
plot(spikes.t[mask_b]/ms, spikes.i[mask_b]-N_ring, '.', color='C3', markersize=1.5, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title(f'Raster completo: {N_RIPETIZIONI} ripetizioni del messaggio su A (x={MESSAGGIO})')
legend(loc='upper right', fontsize=8, markerscale=6)

subplot(3, 1, 2)
t_fine = min(5, N_RIPETIZIONI) * (T_STIM+T_PAUSA)/ms
mask_a_fine = mask_a & (spikes.t/ms < t_fine)
mask_b_fine = mask_b & (spikes.t/ms < t_fine)
plot(spikes.t[mask_a_fine]/ms, spikes.i[mask_a_fine], '.', color='C0', markersize=3, label='A (E)')
plot(spikes.t[mask_b_fine]/ms, spikes.i[mask_b_fine]-N_ring, '.', color='C3', markersize=3, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title('Dettaglio: prime ripetizioni (inizio training)')
legend(loc='upper right', fontsize=8, markerscale=4)

subplot(3, 1, 3)
t_tardi = max(0, N_RIPETIZIONI-5) * (T_STIM+T_PAUSA)/ms
mask_a_tardi = mask_a & (spikes.t/ms > t_tardi)
mask_b_tardi = mask_b & (spikes.t/ms > t_tardi)
plot(spikes.t[mask_a_tardi]/ms, spikes.i[mask_a_tardi], '.', color='C0', markersize=3, label='A (E)')
plot(spikes.t[mask_b_tardi]/ms, spikes.i[mask_b_tardi]-N_ring, '.', color='C3', markersize=3, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title('Dettaglio: ultime ripetizioni (fine training)')
legend(loc='upper right', fontsize=8, markerscale=4)

tight_layout()
savefig('dialogo_spiking.png', dpi=150)
print("\nGrafico salvato in dialogo_spiking.png")

# verifica quantitativa: B risponde piu' spesso, e nella zona giusta, verso
# la fine del training rispetto all'inizio?
def valuta_risposta_B(t0, t1):
    m = mask_b & (spikes.t/ms > t0) & (spikes.t/ms < t1)
    n = m.sum()
    if n == 0:
        return 0, None
    posizioni = xi[spikes.i[m] - N_ring]
    centro = np.angle(np.mean(np.exp(2j*np.pi*posizioni))) / (2*np.pi) % 1.0
    return n, centro

print("\n" + "="*70)
print("VALUTAZIONE: la risposta di B migliora durante il training?")
print("="*70)
n_inizio, c_inizio = valuta_risposta_B(0, 5*(T_STIM+T_PAUSA)/ms)
n_fine, c_fine = valuta_risposta_B((N_RIPETIZIONI-5)*(T_STIM+T_PAUSA)/ms, N_RIPETIZIONI*(T_STIM+T_PAUSA)/ms)
print(f"  Prime 5 ripetizioni: {n_inizio} spike di B, centro = {c_fine if c_inizio is None else f'{c_inizio:.3f}'}")
print(f"  Ultime 5 ripetizioni: {n_fine} spike di B, centro = {'N/A' if c_fine is None else f'{c_fine:.3f}'}")
print(f"  Peso medio finale sinapsi A->B: {np.mean(syn_ab.w)/mV:.3f} mV (iniziale: {w_init/mV:.3f} mV)")
