"""
Cambio di architettura, non solo di parametri, dopo che sei tentativi
indipendenti (5.8, 5.9 con turni/gating su due seed, sotto-soglia,
ablazione) hanno mostrato che il canale B->A non trasmette causalmente
nell'architettura "passiva" (A semplicemente non riceve stimolo diretto
durante i turni B, ma la sua dinamica ricorrente resta libera). Il test
sotto-soglia ha rivelato perche' quell'approccio resta confuso: la
dinamica propria di A (bump piu' largo della zona di stimolo, residui tra
un turno e l'altro) domina qualunque segnale sottile che B potrebbe star
mandando.

Nuova idea, mai provata finora: SILENZIAMENTO ATTIVO. Durante la finestra
di stimolo di ogni turno 'B', ad A viene applicata una forte corrente
iperpolarizzante (non solo l'assenza di stimolo diretto) -- una specie di
"modalita' ascolto forzata" che azzera la sua attivita' residua invece di
lasciarla semplicemente decadere. Alla fine di quella finestra il
silenziamento viene rilasciato, e A e' libera di evolvere secondo
qualunque input riceva dal canale B->A (ormai gia' rinforzato dai turni
precedenti) su una base davvero pulita, senza la propria dinamica interna
residua a confondere il quadro.
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
print(f"Sinapsi inter-anello A->B create: {len(syn_ab)} (peso iniziale {w_init})")

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
print(f"Sinapsi inter-anello B->A create: {len(syn_ba)} (peso iniziale {w_init})")

spikes = SpikeMonitor(neurons)

MESSAGGIO = 0.3
larghezza_stimolo = 0.05
peso_stim = np.exp(-dist_circ(xi - MESSAGGIO)**2 / (2*larghezza_stimolo**2))
gruppo_stim = np.where(peso_stim > 0.5)[0]
print(f"Neuroni stimolati (vicino a x={MESSAGGIO}): {len(gruppo_stim)}")

I_SILENZIO = -8*mV/ms  # abbastanza forte da spegnere qualunque attivita'
                        # residua di A senza distruggere i suoi valori
                        # interni (u, adapt) in modo innaturale -- solo v
                        # viene spinta sotto soglia

N_CICLI = 60
T_STIM = 60*ms
T_PAUSA = 940*ms

print(f"\nTraining: {N_CICLI} turni alternati A/B/A/B/..., A silenziata "
      f"attivamente durante lo stimolo dei turni 'B' "
      f"({(T_STIM+T_PAUSA)*N_CICLI/ms:.0f}ms totali)")

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
        A_E.I[:] = I_SILENZIO   # silenziamento attivo di TUTTA A, non solo
                                  # assenza di stimolo diretto
        A_I.I[:] = I_SILENZIO
        B_E.I[gruppo_stim] = 15*mV/ms
        run(T_STIM)
        B_E.I[:] = 0*mV/ms
        A_E.I[:] = 0*mV/ms      # rilascio: A e' libera di evolvere secondo
        A_I.I[:] = 0*mV/ms      # il solo input ricevuto (rumore + B->A)
    run(T_PAUSA)
    t_corrente += T_STIM + T_PAUSA
    if ciclo % 10 == 0 or ciclo == N_CICLI-1:
        print(f"  ciclo {ciclo:3d}/{N_CICLI} (turno {turno}): "
              f"A->B peso medio={np.mean(syn_ab.w)/mV:.3f}mV max={np.max(syn_ab.w)/mV:.3f}mV  |  "
              f"B->A peso medio={np.mean(syn_ba.w)/mV:.3f}mV max={np.max(syn_ba.w)/mV:.3f}mV")

print("\nTraining completato.")

fig, axes = subplots(3, 1, figsize=(11, 9))

subplot(3, 1, 1)
mask_a = spikes.i < N_E
mask_b = (spikes.i >= N_ring) & (spikes.i < N_ring+N_E)
plot(spikes.t[mask_a]/ms, spikes.i[mask_a], '.', color='C0', markersize=1.5, label='A (E)')
plot(spikes.t[mask_b]/ms, spikes.i[mask_b]-N_ring, '.', color='C3', markersize=1.5, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title(f'Raster completo: {N_CICLI} turni, A silenziata attivamente nei turni B (x={MESSAGGIO})')
legend(loc='upper right', fontsize=8, markerscale=6)

subplot(3, 1, 2)
t_fine = 8 * (T_STIM+T_PAUSA)/ms
mask_a_fine = mask_a & (spikes.t/ms < t_fine)
mask_b_fine = mask_b & (spikes.t/ms < t_fine)
plot(spikes.t[mask_a_fine]/ms, spikes.i[mask_a_fine], '.', color='C0', markersize=3, label='A (E)')
plot(spikes.t[mask_b_fine]/ms, spikes.i[mask_b_fine]-N_ring, '.', color='C3', markersize=3, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title('Dettaglio: primi turni (inizio training)')
legend(loc='upper right', fontsize=8, markerscale=4)

subplot(3, 1, 3)
t_tardi = (N_CICLI-8) * (T_STIM+T_PAUSA)/ms
mask_a_tardi = mask_a & (spikes.t/ms > t_tardi)
mask_b_tardi = mask_b & (spikes.t/ms > t_tardi)
plot(spikes.t[mask_a_tardi]/ms, spikes.i[mask_a_tardi], '.', color='C0', markersize=3, label='A (E)')
plot(spikes.t[mask_b_tardi]/ms, spikes.i[mask_b_tardi]-N_ring, '.', color='C3', markersize=3, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title('Dettaglio: ultimi turni (fine training)')
legend(loc='upper right', fontsize=8, markerscale=4)

tight_layout()
savefig('dialogo_spiking_bidirezionale_silenzia.png', dpi=150)
print("\nGrafico salvato in dialogo_spiking_bidirezionale_silenzia.png")

ciclo_ms = (T_STIM+T_PAUSA)/ms


def valuta_canale(lettera_turno, mask_dest, n_turni_eval, escludi_diretto):
    subset = [t0 for (l, t0) in turni if l == lettera_turno]
    primi = subset[:n_turni_eval]
    ultimi = subset[-n_turni_eval:]

    def conta(lista_inizio):
        n_tot = 0
        posizioni = []
        for t0 in lista_inizio:
            fin = t0 + ciclo_ms
            m = mask_dest & (spikes.t/ms >= t0) & (spikes.t/ms < fin)
            if escludi_diretto:
                m = m & (spikes.t/ms >= t0 + 20)
            idx = spikes.i[m]
            if lettera_turno == 'A':
                idx = idx - N_ring
            n_tot += m.sum()
            posizioni.extend(list(xi[idx % N_E]))
        centro = None
        if posizioni:
            centro = float(np.angle(np.mean(np.exp(2j*np.pi*np.array(posizioni)))) / (2*np.pi) % 1.0)
        return n_tot, centro

    return conta(primi), conta(ultimi)


print("\n" + "="*70)
print("VALUTAZIONE: turni 'A' -- B risponde? (canale A->B)")
print("="*70)
(n_i, c_i), (n_f, c_f) = valuta_canale('A', mask_b, 5, escludi_diretto=False)
print(f"  Primi 5 turni A: {n_i} spike di B, centro={('N/A' if c_i is None else f'{c_i:.3f}')}")
print(f"  Ultimi 5 turni A: {n_f} spike di B, centro={('N/A' if c_f is None else f'{c_f:.3f}')}")
print(f"  Peso medio A->B: {np.mean(syn_ab.w)/mV:.3f}mV (iniziale {w_init/mV:.3f}mV)")

print("\n" + "="*70)
print("VALUTAZIONE: turni 'B' -- A risponde DOPO il rilascio del "
      "silenziamento? (canale B->A, plastico solo in questi turni)")
print("="*70)
(n_i, c_i), (n_f, c_f) = valuta_canale('B', mask_a, 5, escludi_diretto=True)
print(f"  Primi 5 turni B: {n_i} spike di A (>20ms dall'inizio turno, dopo il rilascio), centro={('N/A' if c_i is None else f'{c_i:.3f}')}")
print(f"  Ultimi 5 turni B: {n_f} spike di A (>20ms dall'inizio turno, dopo il rilascio), centro={('N/A' if c_f is None else f'{c_f:.3f}')}")
print(f"  Peso medio B->A: {np.mean(syn_ba.w)/mV:.3f}mV (iniziale {w_init/mV:.3f}mV) max={np.max(syn_ba.w)/mV:.3f}mV")
if n_f > n_i and n_f > 2:
    print("  >>> A mostra attivita' crescente dopo il rilascio del "
          "silenziamento, partendo da una lavagna pulita: possibile "
          "segnale causale reale del canale B->A. <<<")
elif n_i == 0 and n_f == 0:
    print("  >>> Zero spike di A anche con silenziamento attivo e "
          "lavagna pulita: nessuna evidenza di trasmissione causale, "
          "nemmeno in queste condizioni. <<<")
else:
    print("  >>> Segnale debole o non chiaramente crescente. <<<")
