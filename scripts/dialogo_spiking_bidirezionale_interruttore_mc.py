"""
Nono tentativo sulla bidirezionalita' spiking, e primo con un vero cambio
di architettura di circuito (non solo di parametri o di plasticita').

Diagnosi accumulata nelle sezioni 5.8-5.12: il canale B->A, quando
elettricamente sempre attivo insieme al canale A->B (entrambi sempre
"connessi", solo la PLASTICITA' era stata disattivata a turni in 5.9),
o non ha alcun effetto misurabile, o -- se reso abbastanza forte da
averne uno -- fa collassare l'intero sistema in un'esplosione di scarica
continua (5.12: ablazione statica e silenziamento attivo, entrambi
esplosi). L'ipotesi che emerge: il problema non e' la forza del canale
B->A in se', ma il fatto che A->B e B->A restano SEMPRE ENTRAMBI
elettricamente connessi -- due anelli NMDA-sostenuti mutuamente
eccitatori formano un anello di retroazione positiva che, sopra una
certa soglia di guadagno combinato, non ha un punto di equilibrio
stabile.

Soluzione provata qui: un vero INTERRUTTORE DI TRASMISSIONE (non solo di
plasticita'). Una variabile "gate" condivisa per ciascuna sinapsi
inter-anello, che azzera CONTEMPORANEAMENTE sia l'effetto post-sinaptico
(v_post += w*gate) sia l'aggiornamento STDP, non solo l'aggiornamento
STDP come nel gating della sezione 5.9. Durante i turni 'A': gate_ab=1
(A->B trasmette e impara normalmente, come in 5.7), gate_ba=0 (B->A e'
elettricamente SCONNESSA, non solo "non plastica" -- il ciclo di
retroazione e' matematicamente impossibile in questo turno). Durante i
turni 'B': gate_ab=0, gate_ba=1 (il verso opposto). I due anelli non sono
MAI entrambi accoppiati bidirezionalmente nello stesso istante -- ad ogni
dato momento il sistema e' una rete feedforward a senso unico, la stessa
architettura gia' verificata stabile in 5.7.
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

# interruttore di trasmissione: "gate" azzera SIA v_post += w*gate SIA
# l'aggiornamento STDP nello stesso istante -- quando gate=0 la sinapsi
# e' elettricamente sconnessa, non solo "non plastica"
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

spikes = SpikeMonitor(neurons)

MESSAGGIO = 0.3
larghezza_stimolo = 0.05
peso_stim = np.exp(-dist_circ(xi - MESSAGGIO)**2 / (2*larghezza_stimolo**2))
gruppo_stim = np.where(peso_stim > 0.5)[0]
print(f"Neuroni stimolati (vicino a x={MESSAGGIO}): {len(gruppo_stim)}")

N_CICLI = 60
T_STIM = 60*ms
T_PAUSA = 940*ms

print(f"\nTraining: {N_CICLI} turni alternati A/B/A/B/..., interruttore di "
      f"trasmissione: mai entrambi i canali connessi nello stesso istante "
      f"({(T_STIM+T_PAUSA)*N_CICLI/ms:.0f}ms totali)")

turni = []
t_corrente = 0*ms
for ciclo in range(N_CICLI):
    turno = 'A' if ciclo % 2 == 0 else 'B'
    turni.append((turno, t_corrente/ms))
    if turno == 'A':
        syn_ab.gate = 1
        syn_ba.gate = 0
        A_E.I[gruppo_stim] = 15*mV/ms
        run(T_STIM)
        A_E.I[:] = 0*mV/ms
    else:
        syn_ab.gate = 0
        syn_ba.gate = 1
        B_E.I[gruppo_stim] = 15*mV/ms
        run(T_STIM)
        B_E.I[:] = 0*mV/ms
    run(T_PAUSA)
    t_corrente += T_STIM + T_PAUSA
    if ciclo % 10 == 0 or ciclo == N_CICLI-1:
        n_a_finora = int((spikes.i < N_E).sum())
        n_b_finora = int(((spikes.i >= N_ring) & (spikes.i < N_ring+N_E)).sum())
        print(f"  ciclo {ciclo:3d}/{N_CICLI} (turno {turno}): "
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
title(f'Raster completo: {N_CICLI} turni, interruttore di trasmissione (x={MESSAGGIO})')
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
savefig('dialogo_spiking_bidirezionale_interruttore.png', dpi=150)
print("\nGrafico salvato in dialogo_spiking_bidirezionale_interruttore.png")

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
print("VALUTAZIONE: turni 'A' -- B risponde? (canale A->B connesso)")
print("="*70)
(n_i, c_i), (n_f, c_f) = valuta_canale('A', mask_b, 5, escludi_diretto=False)
print(f"  Primi 5 turni A: {n_i} spike di B, centro={('N/A' if c_i is None else f'{c_i:.3f}')}")
print(f"  Ultimi 5 turni A: {n_f} spike di B, centro={('N/A' if c_f is None else f'{c_f:.3f}')}")
print(f"  Peso medio A->B: {np.mean(syn_ab.w)/mV:.3f}mV (iniziale {w_init/mV:.3f}mV)")

print("\n" + "="*70)
print("VALUTAZIONE: turni 'B' -- A risponde? (canale B->A connesso, A->B "
      "completamente sconnesso in questo turno)")
print("="*70)
(n_i, c_i), (n_f, c_f) = valuta_canale('B', mask_a, 5, escludi_diretto=True)
print(f"  Primi 5 turni B: {n_i} spike di A (>20ms dall'inizio turno), centro={('N/A' if c_i is None else f'{c_i:.3f}')}")
print(f"  Ultimi 5 turni B: {n_f} spike di A (>20ms dall'inizio turno), centro={('N/A' if c_f is None else f'{c_f:.3f}')}")
print(f"  Peso medio B->A: {np.mean(syn_ba.w)/mV:.3f}mV (iniziale {w_init/mV:.3f}mV) max={np.max(syn_ba.w)/mV:.3f}mV")

n_a_tot = int(mask_a.sum())
n_b_tot = int(mask_b.sum())
print(f"\nSpike totali sull'intera simulazione: A={n_a_tot}, B={n_b_tot} "
      f"(atteso: centinaia/migliaia se stabile; milioni se esplosione)")
if n_a_tot > 200000 or n_b_tot > 200000:
    print("  >>> ESPLOSIONE: anche con l'interruttore di trasmissione il "
          "sistema e' instabile. <<<")
elif n_f > n_i and n_f > 2:
    print("  >>> Stabile E A mostra attivita' crescente nei turni B: "
          "possibile segnale causale reale, finalmente senza confondimenti "
          "ne' esplosioni. <<<")
elif n_i == 0 and n_f == 0:
    print("  >>> Stabile ma zero spike di A nei turni B: nessuna "
          "trasmissione causale rilevabile, anche eliminando ogni "
          "possibilita' di interferenza dall'altro canale. <<<")
else:
    print("  >>> Stabile, segnale debole o non chiaramente crescente. <<<")

stabile = n_a_tot < 200000 and n_b_tot < 200000
(n_ab_primi, _), (n_ab_ultimi, _) = valuta_canale('B', mask_a, 5, escludi_diretto=True)
print(f"MC_RIGA seed={SEED} stabile={int(stabile)} n_a_tot={n_a_tot} n_b_tot={n_b_tot} "
      f"a_da_b_primi={n_ab_primi} a_da_b_ultimi={n_ab_ultimi}")
