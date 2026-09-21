"""
Prima estensione della comunicazione spiking oltre la coppia: tre
popolazioni (A, B, C) invece di due. Mai tentato ne' nel modello a rate
ne' in quello spiking in questo progetto.

Per evitare di riaprire l'intera saga di stabilizzazione delle sezioni
5.8-5.19 (necessaria solo per la bidirezionalita' vera, dove due anelli
sono mutuamente eccitatori), questa prima rete a tre popolazioni usa
SOLO collegamenti unidirezionali con STDP vera -- esattamente il
meccanismo gia' collaudato e stabile di dialogo_spiking.py (sezione 5.7),
ripetuto tre volte per formare una STAFFETTA DIREZIONALE ad anello:
A->B->C->A. Nessun anello riceve MAI input da due sorgenti
contemporaneamente (ogni anello ha esattamente un collegamento in
ingresso, con STDP, e uno in uscita), quindi nessuna delle instabilita'
di retroazione mutua incontrate con la vera bidirezionalita' puo'
presentarsi qui per costruzione.

Domanda sperimentale: un messaggio iniettato SOLO in A (mai direttamente
in B o C) si propaga correttamente attraverso la catena, preservando
l'informazione di posizione dopo due salti (in C) e persino dopo il giro
completo (di nuovo in A, tramite C->A)?
"""

from brian2 import *

start_scope()

N_E = 100
N_I = 25
N_ring = N_E + N_I
N_POP = 3  # A, B, C

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

neurons = NeuronGroup(N_POP*N_ring, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV + 5*mV*randn(N_POP*N_ring)
neurons.u = b * neurons.v[:]
neurons.I = 0*mV/ms
neurons.g_nmda = 0*mV/ms

def sotto(pop_idx):
    base = pop_idx * N_ring
    return neurons[base:base+N_E], neurons[base+N_E:base+N_ring]

A_E, A_I = sotto(0)
B_E, B_I = sotto(1)
C_E, C_I = sotto(2)

xi = np.arange(N_E) / N_E
A_E.x = xi
B_E.x = xi
C_E.x = xi

background = PoissonGroup(N_POP*N_ring, rates=400*Hz)
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
syn_C = costruisci_anello(C_E, C_I, seed=2)

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


def collega(SRC_E, DST_E, seed_conn, w_iniziale):
    dmat_inter = dist_circ(xi[:, None] - xi[None, :])
    prob_inter = p_inter_picco * np.exp(-dmat_inter**2 / (2*sigma_inter**2))
    rng_inter = np.random.RandomState(seed_conn)
    mask_inter = rng_inter.rand(N_E, N_E) < prob_inter
    i_inter, j_inter = np.nonzero(mask_inter)
    syn = Synapses(SRC_E, DST_E, model=stdp_model, on_pre=stdp_on_pre, on_post=stdp_on_post)
    syn.connect(i=i_inter, j=j_inter)
    syn.w = w_iniziale
    return syn


# A->B parte dal valore originale di 5.7 (A ha uno stimolo esterno diretto
# forte, si innesca facilmente). B->C e C->A partono da un peso piu' alto:
# la loro sorgente (B, poi C) e' essa stessa un eco di secondo/terzo
# livello, con un segnale molto piu' debole e sparso di quello diretto su
# A -- lo stesso identico problema "uovo-gallina" gia' risolto in 5.7 per
# A->B (tentativo 1: w_init=0.3mV, B non si accendeva mai), qui riapplicato
# ai salti successivi della staffetta
syn_ab = collega(A_E, B_E, seed_conn=10, w_iniziale=w_init)
syn_bc = collega(B_E, C_E, seed_conn=11, w_iniziale=3.0*mV)
syn_ca = collega(C_E, A_E, seed_conn=12, w_iniziale=3.0*mV)
print(f"Staffetta A->B->C->A creata: A->B {len(syn_ab)} sinapsi, "
      f"B->C {len(syn_bc)} sinapsi, C->A {len(syn_ca)} sinapsi")

spikes = SpikeMonitor(neurons)

MESSAGGIO = 0.3
larghezza_stimolo = 0.05
peso_stim = np.exp(-dist_circ(xi - MESSAGGIO)**2 / (2*larghezza_stimolo**2))
gruppo_stim = np.where(peso_stim > 0.5)[0]
print(f"Neuroni stimolati SOLO in A (vicino a x={MESSAGGIO}): {len(gruppo_stim)}")

N_RIPETIZIONI_NORMALI = 80  # allena A->B->C come nella sezione 5.20
N_SONDE_SILENZIOSE = 10     # poi, cicli CONSECUTIVI senza alcuno stimolo
                              # su A -- se C->A e' causale, A dovrebbe
                              # mostrare attivita' tardiva stabile in
                              # OGNUNO di questi cicli (non un semplice
                              # decadimento della propria persistenza,
                              # che invece si esaurirebbe entro 1-2 cicli)
T_STIM = 60*ms
T_PAUSA = 440*ms  # piu' lunga di 5.7 (340ms) per dare tempo alla staffetta
                    # di propagare fino a C (due salti) prima del ciclo dopo
N_RIPETIZIONI = N_RIPETIZIONI_NORMALI + N_SONDE_SILENZIOSE

print(f"\nTraining: {N_RIPETIZIONI_NORMALI} ripetizioni normali (stimolo su A), "
      f"poi {N_SONDE_SILENZIOSE} cicli sonda SENZA alcuno stimolo esterno "
      f"({(T_STIM+T_PAUSA)*N_RIPETIZIONI/ms:.0f}ms totali)")

tipo_ciclo = []  # 'normale' o 'sonda', per ciclo
for rep in range(N_RIPETIZIONI):
    if rep < N_RIPETIZIONI_NORMALI:
        tipo_ciclo.append('normale')
        A_E.I[gruppo_stim] = 15*mV/ms
        run(T_STIM)
        A_E.I[:] = 0*mV/ms
    else:
        tipo_ciclo.append('sonda')
        run(T_STIM)  # nessuno stimolo esterno su nessun anello in questo ciclo
    run(T_PAUSA)
    if rep % 10 == 0 or rep == N_RIPETIZIONI-1 or tipo_ciclo[-1] == 'sonda':
        print(f"  ripetizione {rep:3d}/{N_RIPETIZIONI} ({tipo_ciclo[-1]}): "
              f"A->B={np.mean(syn_ab.w)/mV:.2f}mV  B->C={np.mean(syn_bc.w)/mV:.2f}mV  "
              f"C->A={np.mean(syn_ca.w)/mV:.2f}mV")

print("\nTraining completato.")

fig, axes = subplots(2, 1, figsize=(11, 7))
subplot(2, 1, 1)
colori = ['C0', 'C3', 'C2']
etichette = ['A (E)', 'B (E)', 'C (E)']
for k in range(N_POP):
    base = k*N_ring
    m = (spikes.i >= base) & (spikes.i < base+N_E)
    plot(spikes.t[m]/ms, spikes.i[m]-base, '.', color=colori[k], markersize=1.5, label=etichette[k])
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title(f'Raster completo: staffetta A->B->C->A, {N_RIPETIZIONI} ripetizioni (x={MESSAGGIO})')
legend(loc='upper right', fontsize=8, markerscale=6)

subplot(2, 1, 2)
t_tardi = (N_RIPETIZIONI-6) * (T_STIM+T_PAUSA)/ms
for k in range(N_POP):
    base = k*N_ring
    m = (spikes.i >= base) & (spikes.i < base+N_E) & (spikes.t/ms > t_tardi)
    plot(spikes.t[m]/ms, spikes.i[m]-base, '.', color=colori[k], markersize=3, label=etichette[k])
xlabel('Tempo (ms)')
ylabel('Indice neurone')
title('Dettaglio: ultime ripetizioni (fine training)')
legend(loc='upper right', fontsize=8, markerscale=4)

tight_layout()
savefig('dialogo_spiking_relay_3anelli.png', dpi=150)
print("\nGrafico salvato in dialogo_spiking_relay_3anelli.png")

ciclo_ms = (T_STIM+T_PAUSA)/ms


def valuta_pop(pop_idx, n_eval):
    base = pop_idx * N_ring
    mask_pop = (spikes.i >= base) & (spikes.i < base+N_E)

    def conta(rep_range):
        n_tot, posizioni = 0, []
        for rep in rep_range:
            t0 = rep * ciclo_ms
            fin = t0 + ciclo_ms
            m = mask_pop & (spikes.t/ms >= t0) & (spikes.t/ms < fin)
            idx = spikes.i[m] - base
            n_tot += m.sum()
            posizioni.extend(list(xi[idx]))
        centro = None
        if posizioni:
            centro = float(np.angle(np.mean(np.exp(2j*np.pi*np.array(posizioni)))) / (2*np.pi) % 1.0)
        return n_tot, centro

    primi = range(0, n_eval)
    ultimi = range(N_RIPETIZIONI-n_eval, N_RIPETIZIONI)
    return conta(primi), conta(ultimi)


print("\n" + "="*70)
print("VALUTAZIONE: l'informazione di posizione sopravvive alla staffetta?")
print("="*70)
for k, nome in enumerate(['A (sorgente diretta)', 'B (un salto: A->B)', 'C (due salti: A->B->C)']):
    (n_i, c_i), (n_f, c_f) = valuta_pop(k, 5)
    err_i = dist_circ(c_i - MESSAGGIO) if c_i is not None else None
    err_f = dist_circ(c_f - MESSAGGIO) if c_f is not None else None
    print(f"  {nome}:")
    print(f"    prime 5 rip.: {n_i} spike, centro={('N/A' if c_i is None else f'{c_i:.3f}')}, "
          f"errore={('N/A' if err_i is None else f'{err_i:.3f}')}")
    print(f"    ultime 5 rip.: {n_f} spike, centro={('N/A' if c_f is None else f'{c_f:.3f}')}, "
          f"errore={('N/A' if err_f is None else f'{err_f:.3f}')}")

print(f"\nPesi finali: A->B={np.mean(syn_ab.w)/mV:.3f}mV (iniziale {w_init/mV:.3f}mV), "
      f"B->C={np.mean(syn_bc.w)/mV:.3f}mV, C->A={np.mean(syn_ca.w)/mV:.3f}mV")

# --- Domanda centrale dell'esperimento: durante i cicli sonda (nessuno
# stimolo esterno su nessun anello), A mostra attivita' tardiva STABILE
# in OGNUNO dei 10 cicli (segno di un loop causale C->A che si
# autosostiene), o l'attivita' decade entro 1-2 cicli (segno che era
# solo persistenza residua di A stessa, non un vero segnale chiuso
# dall'anello)? ---
print("\n" + "="*70)
print("DETTAGLIO PER CICLO SONDA: attivita' di A senza alcuno stimolo esterno")
print("="*70)
base_A = 0
mask_A = (spikes.i >= base_A) & (spikes.i < base_A+N_E)
for rep in range(N_RIPETIZIONI_NORMALI, N_RIPETIZIONI):
    t0 = rep * ciclo_ms
    fin = t0 + ciclo_ms
    m = mask_A & (spikes.t/ms >= t0) & (spikes.t/ms < fin)
    idx = spikes.i[m] - base_A
    n_tot = m.sum()
    centro = None
    if n_tot > 0:
        centro = float(np.angle(np.mean(np.exp(2j*np.pi*np.array(xi[idx])))) / (2*np.pi) % 1.0)
    err = dist_circ(centro - MESSAGGIO) if centro is not None else None
    indice_sonda = rep - N_RIPETIZIONI_NORMALI
    print(f"  sonda {indice_sonda+1:2d}/{N_SONDE_SILENZIOSE}: {n_tot} spike in A, "
          f"centro={('N/A' if centro is None else f'{centro:.3f}')}, "
          f"errore={('N/A' if err is None else f'{err:.3f}')}")

# controllo di stabilita' (per costruzione non dovrebbe mai esplodere,
# ma verifichiamo comunque per coerenza con il resto del progetto)
n_tot_rete = len(spikes.i)
print(f"\nSpike totali sull'intera rete: {n_tot_rete} "
      f"(atteso migliaia se stabile; milioni se esplosione)")
print(">>> ESPLOSIONE <<<" if n_tot_rete > 300000 else ">>> Stabile <<<")
