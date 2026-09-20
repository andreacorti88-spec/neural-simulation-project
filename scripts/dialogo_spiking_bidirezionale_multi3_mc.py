"""
Quattordicesimo tentativo: stress-test della sezione 5.18 (due messaggi
bidirezionali, 8/8 seed stabili) scalando a TRE messaggi -- lo stesso
irrigidimento del carico gia' affrontato dal modello a rate quando e'
passato da due a tre messaggi (dialogo_multi_repulsione_3msg.py, sezione
5.5.8), dove l'interferenza si era fatta sentire per la prima volta in
modo serio. Stesse tre posizioni di quel test: x=0.10, 0.40, 0.70.
Schema a turni esteso a 6 combinazioni cicliche (messaggio, direzione):
(m0,A),(m0,B),(m1,A),(m1,B),(m2,A),(m2,B), ripetuto 20 volte (120 turni,
120 secondi di simulazione).
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

tau_adapt_pop = 200*ms
gain_adapt_pop = 0.15*mV/ms

tau_rate_hat = 200*ms  # v2: finestra piu' corta (era 500ms) -- reagisce piu'
                        # in fretta a un'accensione che sta per scappare

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I + g_nmda - adapt : volt
du/dt = a*(b*v - u) : volt/second
dg_nmda/dt = -g_nmda/tau_nmda : volt/second
dadapt/dt = -adapt/tau_adapt_pop : volt/second
dratehat/dt = -ratehat/tau_rate_hat : Hz
I : volt/second
x : 1
'''

neurons = NeuronGroup(2*N_ring, eqs, threshold='v > 30*mV',
                       reset='v = c; u += d; adapt += gain_adapt_pop; ratehat += 1/tau_rate_hat',
                       method='euler')
neurons.v = -65*mV + 5*mV*randn(2*N_ring)
neurons.u = b * neurons.v[:]
neurons.I = 0*mV/ms
neurons.g_nmda = 0*mV/ms
neurons.adapt = 0*mV/ms
neurons.ratehat = 0*Hz

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


syn_ee_A, syn_ei_A, syn_ie_A, syn_ii_A = costruisci_anello(A_E, A_I, seed=1000+SEED)
syn_ee_B, syn_ei_B, syn_ie_B, syn_ii_B = costruisci_anello(B_E, B_I, seed=2000+SEED)

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
gate : 1 (shared)
'''
stdp_on_pre = '''
v_post += w*gate
apre += A_ltp*gate
w = clip(w+apost*gate, w_min, w_max)
'''
stdp_on_post = '''
apost -= A_ltd*gate
w = clip(w+apre*gate, w_min, w_max)
'''

syn_ab = Synapses(A_E, B_E, model=stdp_model, on_pre=stdp_on_pre, on_post=stdp_on_post)
dmat_inter = dist_circ(xi[:, None] - xi[None, :])
prob_inter = p_inter_picco * np.exp(-dmat_inter**2 / (2*sigma_inter**2))
rng_inter = np.random.RandomState(3000+SEED)
mask_inter = rng_inter.rand(N_E, N_E) < prob_inter
i_inter, j_inter = np.nonzero(mask_inter)
syn_ab.connect(i=i_inter, j=j_inter)
syn_ab.w = w_init
syn_ab.gate = 1
print(f"Sinapsi inter-anello A->B create: {len(syn_ab)} (peso iniziale {w_init})")

syn_ba = Synapses(B_E, A_E, model=stdp_model, on_pre=stdp_on_pre, on_post=stdp_on_post)
mask_inter_ba = rng_inter.rand(N_E, N_E) < prob_inter
i_inter_ba, j_inter_ba = np.nonzero(mask_inter_ba)
syn_ba.connect(i=i_inter_ba, j=j_inter_ba)
syn_ba.w = w_init
syn_ba.gate = 0
print(f"Sinapsi inter-anello B->A create: {len(syn_ba)} (peso iniziale {w_init})")

# --- normalizzazione omeostatica: ogni 50ms, le sinapsi inter-anello in
# ingresso a ciascun neurone E vengono scalate per spingere il suo tasso
# di scarica recente (ratehat) verso target_rate. Fattore per aggiornamento
# limitato a +-3% per evitare oscillazioni, indipendente dal tetto w_max
# gia' esistente (che resta come sicurezza aggiuntiva) ---
target_rate = 2.5*Hz  # v2: bersaglio leggermente piu' conservativo (era 3Hz)
fattore_min, fattore_max = 0.99, 1.01  # v2: aggiornamenti piu' piccoli...


soglia_emergenza = 15*Hz  # v3: interruttore di emergenza -- molto sopra il
                            # bersaglio (2.5Hz) ma ben sotto un vero runaway
                            # (centinaia di Hz), per intercettare un'accensione
                            # che la correzione graduale non riesce a fermare


@network_operation(dt=20*ms)  # ...ma piu' frequenti (era 50ms) -- stesso
def omeostasi():              # margine di correzione nel tempo, reazione piu' rapida
    rate_b = np.asarray(B_E.ratehat/Hz)
    rate_a = np.asarray(A_E.ratehat/Hz)
    fattore_b = np.clip((target_rate/Hz) / (rate_b + 0.5), fattore_min, fattore_max)
    fattore_a = np.clip((target_rate/Hz) / (rate_a + 0.5), fattore_min, fattore_max)
    # interruttore di emergenza: azzeramento immediato (non graduale) delle
    # sinapsi in ingresso ai neuroni che superano la soglia critica
    emergenza_b = rate_b > (soglia_emergenza/Hz)
    emergenza_a = rate_a > (soglia_emergenza/Hz)
    if emergenza_b.any():
        fattore_b[emergenza_b] = 0.0
        B_E.adapt[emergenza_b] = np.maximum(np.asarray(B_E.adapt[emergenza_b]/(mV/ms)), 20.0) * (mV/ms)
        print(f"    [emergenza a t={defaultclock.t/ms:.0f}ms: {emergenza_b.sum()} neuroni di B azzerati+soppressi, ratehat max={rate_b.max():.1f}Hz]")
    if emergenza_a.any():
        fattore_a[emergenza_a] = 0.0
        A_E.adapt[emergenza_a] = np.maximum(np.asarray(A_E.adapt[emergenza_a]/(mV/ms)), 20.0) * (mV/ms)
        print(f"    [emergenza a t={defaultclock.t/ms:.0f}ms: {emergenza_a.sum()} neuroni di A azzerati+soppressi, ratehat max={rate_a.max():.1f}Hz]")
    j_ab = np.asarray(syn_ab.j[:])
    j_ba = np.asarray(syn_ba.j[:])
    syn_ab.w = np.clip(np.asarray(syn_ab.w/mV) * fattore_b[j_ab], w_min/mV, w_max/mV) * mV
    syn_ba.w = np.clip(np.asarray(syn_ba.w/mV) * fattore_a[j_ba], w_min/mV, w_max/mV) * mV


spikes = SpikeMonitor(neurons)
rate_mon = StateMonitor(neurons, 'ratehat', record=[0, N_ring], dt=50*ms)

MESSAGGI = [0.10, 0.40, 0.70]  # stesse posizioni del modello a rate a 3 messaggi
                                 # (dialogo_multi_repulsione_3msg.py, sezione 5.5.8)
larghezza_stimolo = 0.05
gruppi_stim = []
for m in MESSAGGI:
    peso_stim = np.exp(-dist_circ(xi - m)**2 / (2*larghezza_stimolo**2))
    gruppi_stim.append(np.where(peso_stim > 0.5)[0])
    print(f"Messaggio x={m}: {len(gruppi_stim[-1])} neuroni stimolati")

N_CICLI = 120  # 20 ripetizioni per ciascuna delle 6 combinazioni (3 messaggi x 2 direzioni)
T_STIM = 60*ms
T_PAUSA = 940*ms

print(f"\nTraining: {N_CICLI} turni, ciclando su {len(MESSAGGI)} messaggi x 2 "
      f"direzioni (A/B), interruttore di trasmissione + omeostasi + "
      f"soppressione di emergenza ({(T_STIM+T_PAUSA)*N_CICLI/ms:.0f}ms totali)")

turni = []  # (messaggio_idx, direzione, tempo_inizio_ms)
combinazioni = [(k, dirz) for k in range(len(MESSAGGI)) for dirz in ('A', 'B')]
t_corrente = 0*ms
for ciclo in range(N_CICLI):
    k, direzione = combinazioni[ciclo % len(combinazioni)]
    turni.append((k, direzione, t_corrente/ms))
    if direzione == 'A':
        syn_ab.gate = 1
        syn_ba.gate = 0
        A_E.I[gruppi_stim[k]] = 15*mV/ms
        run(T_STIM)
        A_E.I[:] = 0*mV/ms
    else:
        syn_ab.gate = 0
        syn_ba.gate = 1
        B_E.I[gruppi_stim[k]] = 15*mV/ms
        run(T_STIM)
        B_E.I[:] = 0*mV/ms
    run(T_PAUSA)
    t_corrente += T_STIM + T_PAUSA
    if ciclo % 10 == 0 or ciclo == N_CICLI-1:
        n_a_finora = int((spikes.i < N_E).sum())
        n_b_finora = int(((spikes.i >= N_ring) & (spikes.i < N_ring+N_E)).sum())
        print(f"  ciclo {ciclo:3d}/{N_CICLI} (msg x={MESSAGGI[k]}, turno {direzione}): "
              f"A->B peso medio={np.mean(syn_ab.w)/mV:.3f}mV max={np.max(syn_ab.w)/mV:.3f}mV  |  "
              f"B->A peso medio={np.mean(syn_ba.w)/mV:.3f}mV max={np.max(syn_ba.w)/mV:.3f}mV  |  "
              f"spike finora A={n_a_finora} B={n_b_finora}")

print("\nTraining completato.")

fig, axes = subplots(3, 1, figsize=(11, 9))

subplot(3, 1, 1)
mask_a = spikes.i < N_E
mask_b = (spikes.i >= N_ring) & (spikes.i < N_ring+N_E)
plot(spikes.t[mask_a]/ms, spikes.i[mask_a], '.', color='C0', markersize=1.5, label='A (E)')
plot(spikes.t[mask_b]/ms, spikes.i[mask_b]-N_ring, '.', color='C3', markersize=1.5, label='B (E)')
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title(f'Raster completo: {N_CICLI} turni, {len(MESSAGGI)} messaggi x=[{", ".join(str(m) for m in MESSAGGI)}]')
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
savefig('dialogo_spiking_bidirezionale_multi3.png', dpi=150)
print("\nGrafico salvato in dialogo_spiking_bidirezionale_multi3.png")

ciclo_ms = (T_STIM+T_PAUSA)/ms


def valuta_canale(k, direzione, mask_dest, n_turni_eval, escludi_diretto):
    subset = [t0 for (kk, d, t0) in turni if kk == k and d == direzione]
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
            if direzione == 'A':
                idx = idx - N_ring
            n_tot += m.sum()
            posizioni.extend(list(xi[idx % N_E]))
        centro = None
        if posizioni:
            centro = float(np.angle(np.mean(np.exp(2j*np.pi*np.array(posizioni)))) / (2*np.pi) % 1.0)
        return n_tot, centro

    return conta(primi), conta(ultimi)


for k, m in enumerate(MESSAGGI):
    print("\n" + "="*70)
    print(f"MESSAGGIO x={m}")
    print("="*70)
    (n_i, c_i), (n_f, c_f) = valuta_canale(k, 'A', mask_b, 5, escludi_diretto=False)
    err_i = dist_circ(c_i - m) if c_i is not None else None
    err_f = dist_circ(c_f - m) if c_f is not None else None
    print(f"  Turni A (B risponde): primi 5={n_i} spike (err={('N/A' if err_i is None else f'{err_i:.3f}')}), "
          f"ultimi 5={n_f} spike (err={('N/A' if err_f is None else f'{err_f:.3f}')})")
    (n_i, c_i), (n_f, c_f) = valuta_canale(k, 'B', mask_a, 5, escludi_diretto=True)
    err_i = dist_circ(c_i - m) if c_i is not None else None
    err_f = dist_circ(c_f - m) if c_f is not None else None
    print(f"  Turni B (A risponde): primi 5={n_i} spike (err={('N/A' if err_i is None else f'{err_i:.3f}')}), "
          f"ultimi 5={n_f} spike (err={('N/A' if err_f is None else f'{err_f:.3f}')})")

n_a_tot = int(mask_a.sum())
n_b_tot = int(mask_b.sum())
print(f"\nSpike totali sull'intera simulazione: A={n_a_tot}, B={n_b_tot} "
      f"(atteso: centinaia/migliaia se stabile; milioni se esplosione)")
if n_a_tot > 200000 or n_b_tot > 200000:
    print("  >>> ESPLOSIONE. <<<")
else:
    print("  >>> Stabile. <<<")

n_a_tot = int(mask_a.sum())
n_b_tot = int(mask_b.sum())
stabile = n_a_tot < 200000 and n_b_tot < 200000
righe_msg = []
for k, m in enumerate(MESSAGGI):
    (n_i_ab, c_i_ab), (n_f_ab, c_f_ab) = valuta_canale(k, 'A', mask_b, 5, escludi_diretto=False)
    (n_i_ba, c_i_ba), (n_f_ba, c_f_ba) = valuta_canale(k, 'B', mask_a, 5, escludi_diretto=True)
    err_ab = dist_circ(c_f_ab - m) if c_f_ab is not None else 1.0
    err_ba = dist_circ(c_f_ba - m) if c_f_ba is not None else 1.0
    righe_msg.append(f"m{k}(x={m}):ab_i={n_i_ab},ab_f={n_f_ab},err_ab={err_ab:.3f},ba_i={n_i_ba},ba_f={n_f_ba},err_ba={err_ba:.3f}")
print(f"MC_RIGA seed={SEED} stabile={int(stabile)} n_a_tot={n_a_tot} n_b_tot={n_b_tot} " + " ".join(righe_msg))
