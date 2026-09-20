"""
Collega il progetto Geant4 cnaoRingImpact (sezione 5.21) alla rete
spiking: i 10 neuroni con la dose fisica piu' alta nello scenario
clinico piu' severo simulato (carbonio-12 a 120 MeV/u, quello con la
dose piu' concentrata) vengono "spenti" -- tutte le loro sinapsi in
USCITA rimosse (intra-anello E-E/E-I e, per l'anello A, anche il
canale inter-anello A->B), il modello funzionale piu' diretto di un
neurone che non puo' piu' influenzare nessuno, indipendentemente da
cosa riceva in ingresso. Applicato simmetricamente a entrambi gli
anelli A e B (rappresentano la stessa architettura fisica).

Base: dialogo_spiking_multi.py (sezione 5.10, due messaggi
unidirezionali A->B, gia' verificato robusto su 8 seed in 5.11).
Domanda: la comunicazione multi-messaggio regge con questi 10 neuroni
(su 100, 10% dell'anello) resi inerti? Nota di rilievo: il neurone 65
(x=0.65) e' tra i dieci -- coincide esattamente con una delle
posizioni usate per i messaggi.
"""

KNOCKOUT = sorted([8, 16, 45, 50, 58, 59, 60, 65, 66, 72])  # top-10 dose,
                                                              # carbonio-12 120 MeV/u
                                                              # (geant4_cnao_ring/build/results)

"""
Estensione di dialogo_spiking.py (sezione 5.7, comunicazione appresa A->B
con STDP vera su un solo messaggio) a DUE messaggi, sulle stesse posizioni
usate nel primo test di interferenza del modello a rate (dialogo_multi_
messaggio.py, sezione 5.5.4): x=0.15 e x=0.65.

Nota di partenza importante, diversa dal modello a rate: li' i pesi
inter-popolazione partivano CASUALI e la scoperta di una mappa richiedeva
un processo di ricerca (prima "rumore puro", poi "rumore + repulsione",
poi "rumore + scoperta locale") -- da cui nasceva l'interferenza quando
due messaggi cercavano di scoprire la loro mappa nello stesso spazio
casuale contemporaneamente. Qui invece la connettivita' inter-anello parte
gia' strutturata (caduta gaussiana sulla distanza circolare, offset zero,
sezione 5.7): non serve alcuna scoperta, la mappa identita' e' gia'
presente prima che la STDP inizi ad allenare qualunque cosa. Non e' quindi
scontato che la stessa interferenza si ripresenti -- il meccanismo che la
causava nel modello a rate potrebbe semplicemente non esistere qui.
Bidirezionalita' non ritestata in questo script: la sezione 5.9 ha gia'
concluso che il canale B->A non trasmette causalmente in questo sistema,
quindi ci si concentra sul canale A->B, l'unico verificato funzionante.
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


def costruisci_anello(E, I_pop, seed, knockout=()):
    dmat = dist_circ(xi[:, None] - xi[None, :])
    prob_ee = np.exp(-dmat**2 / (2*sigma_exc**2))
    prob_ee[dmat > CUTOFF_LOCALE] = 0.0
    np.fill_diagonal(prob_ee, 0)
    rng = np.random.RandomState(seed)
    mask_ee = rng.rand(N_E, N_E) < prob_ee
    i_idx, j_idx = np.nonzero(mask_ee)

    ko = np.asarray(list(knockout), dtype=int)
    # "spento" = nessuna sinapsi in USCITA sopravvive, da nessun
    # neurone in knockout -- puo' ancora ricevere input, ma non puo'
    # piu' influenzare nessuno (l'unica cosa che conta per una rete
    # a spike: se non spara mai in modo efficace a valle, e' come
    # se non ci fosse)
    tieni_ee = ~np.isin(i_idx, ko)
    i_idx, j_idx = i_idx[tieni_ee], j_idx[tieni_ee]

    syn_ee = Synapses(E, E, on_pre='g_nmda_post += w_ee_nmda')
    syn_ee.connect(i=i_idx, j=j_idx)
    syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
    # Costruita esplicitamente (invece di .connect(p=...)) per poter
    # filtrare via i pre-sinaptici spenti, come sopra per E->E.
    rng_ei = np.random.RandomState(seed + 100)
    mask_ei = rng_ei.rand(len(E), len(I_pop)) < p_conn
    i_ei, j_ei = np.nonzero(mask_ei)
    tieni_ei = ~np.isin(i_ei, ko)
    syn_ei.connect(i=i_ei[tieni_ei], j=j_ei[tieni_ei])
    syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
    syn_ie.connect(p=p_conn)
    syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
    syn_ii.connect(condition='i!=j', p=p_conn)
    return syn_ee, syn_ei, syn_ie, syn_ii


syn_A = costruisci_anello(A_E, A_I, seed=0, knockout=KNOCKOUT)
syn_B = costruisci_anello(B_E, B_I, seed=1, knockout=KNOCKOUT)
print(f"Knockout applicato: {len(KNOCKOUT)} neuroni per anello ({KNOCKOUT}) -- "
      f"tutte le sinapsi in uscita rimosse (E->E e E->I)")

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
rng_inter = np.random.RandomState(2)
mask_inter = rng_inter.rand(N_E, N_E) < prob_inter
i_inter, j_inter = np.nonzero(mask_inter)
ko_arr = np.asarray(KNOCKOUT, dtype=int)
tieni_inter = ~np.isin(i_inter, ko_arr)  # un neurone A spento non puo'
                                            # trasmettere neanche a B
i_inter, j_inter = i_inter[tieni_inter], j_inter[tieni_inter]
syn_ab.connect(i=i_inter, j=j_inter)
syn_ab.w = w_init
print(f"Sinapsi inter-anello A->B create: {len(syn_ab)} (peso iniziale {w_init}, "
      f"dopo aver rimosso le uscite dai {len(KNOCKOUT)} neuroni spenti)")

spikes = SpikeMonitor(neurons)

MESSAGGI = [0.15, 0.65]
larghezza_stimolo = 0.05
gruppi_stim = []
for m in MESSAGGI:
    peso_stim = np.exp(-dist_circ(xi - m)**2 / (2*larghezza_stimolo**2))
    gruppi_stim.append(np.where(peso_stim > 0.5)[0])
    print(f"Messaggio x={m}: {len(gruppi_stim[-1])} neuroni di A stimolati")

N_RIPETIZIONI_TOT = 60  # round-robin: 30 ripetizioni per messaggio
T_STIM = 60*ms
T_PAUSA = 340*ms

print(f"\nTraining: {N_RIPETIZIONI_TOT} ripetizioni round-robin tra "
      f"{len(MESSAGGI)} messaggi ({(T_STIM+T_PAUSA)*N_RIPETIZIONI_TOT/ms:.0f}ms totali)")

storico = []  # (indice_messaggio, tempo_inizio_ms)
t_corrente = 0*ms
for rep in range(N_RIPETIZIONI_TOT):
    k = rep % len(MESSAGGI)
    storico.append((k, t_corrente/ms))
    A_E.I[gruppi_stim[k]] = 15*mV/ms
    run(T_STIM)
    A_E.I[:] = 0*mV/ms
    run(T_PAUSA)
    t_corrente += T_STIM + T_PAUSA
    if rep % 10 == 0 or rep == N_RIPETIZIONI_TOT-1:
        print(f"  ripetizione {rep:3d}/{N_RIPETIZIONI_TOT} (messaggio x={MESSAGGI[k]}): "
              f"peso medio={np.mean(syn_ab.w)/mV:.3f}mV max={np.max(syn_ab.w)/mV:.3f}mV")

print("\nTraining completato.")

fig, axes = subplots(2, 1, figsize=(11, 7))
subplot(2, 1, 1)
mask_a = spikes.i < N_E
mask_b = (spikes.i >= N_ring) & (spikes.i < N_ring+N_E)
plot(spikes.t[mask_a]/ms, spikes.i[mask_a], '.', color='C0', markersize=1.5, label='A (E)')
plot(spikes.t[mask_b]/ms, spikes.i[mask_b]-N_ring, '.', color='C3', markersize=1.5, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title(f'Raster completo: {N_RIPETIZIONI_TOT} ripetizioni, messaggi x={MESSAGGI}, '
      f'{len(KNOCKOUT)} neuroni spenti (dose CNAO piu\' alta)')
legend(loc='upper right', fontsize=8, markerscale=6)

subplot(2, 1, 2)
t_tardi = (N_RIPETIZIONI_TOT-10) * (T_STIM+T_PAUSA)/ms
mask_a_tardi = mask_a & (spikes.t/ms > t_tardi)
mask_b_tardi = mask_b & (spikes.t/ms > t_tardi)
plot(spikes.t[mask_a_tardi]/ms, spikes.i[mask_a_tardi], '.', color='C0', markersize=3, label='A (E)')
plot(spikes.t[mask_b_tardi]/ms, spikes.i[mask_b_tardi]-N_ring, '.', color='C3', markersize=3, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title('Dettaglio: ultime ripetizioni (fine training)')
legend(loc='upper right', fontsize=8, markerscale=4)

tight_layout()
savefig('dialogo_spiking_cnao_knockout.png', dpi=150)
print("\nGrafico salvato in dialogo_spiking_cnao_knockout.png")

ciclo_ms = (T_STIM+T_PAUSA)/ms


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


print("\n" + "="*70)
print("VALUTAZIONE PER MESSAGGIO: la risposta di B si localizza "
      "correttamente e si rafforza per ENTRAMBI i messaggi?")
print("="*70)
for k, m in enumerate(MESSAGGI):
    (n_i, c_i), (n_f, c_f) = valuta_messaggio(k, 5)
    err_i = dist_circ(c_i - m) if c_i is not None else None
    err_f = dist_circ(c_f - m) if c_f is not None else None
    print(f"  Messaggio x={m}:")
    print(f"    prime 5 ripetizioni: {n_i} spike, centro={('N/A' if c_i is None else f'{c_i:.3f}')}, "
          f"errore={('N/A' if err_i is None else f'{err_i:.3f}')}")
    print(f"    ultime 5 ripetizioni: {n_f} spike, centro={('N/A' if c_f is None else f'{c_f:.3f}')}, "
          f"errore={('N/A' if err_f is None else f'{err_f:.3f}')}")

print(f"\nPeso medio finale A->B: {np.mean(syn_ab.w)/mV:.3f}mV (iniziale {w_init/mV:.3f}mV)")
print(f"Peso massimo A->B: {np.max(syn_ab.w)/mV:.3f}mV")
