"""
La sezione 5.20 ha dimostrato che l'informazione di posizione sopravvive
alla staffetta A->B->C->A per due salti (fino a C). Resta aperta la
domanda se il terzo salto, C->A, chiude davvero il cerchio -- cioe' se
C esercita un'influenza causale misurabile su A, oltre alla risposta di
A al proprio stimolo diretto.

Rischio di confondimento: in questo script (come in 5.20) i neuroni non
hanno adattamento di popolazione (le stesse equazioni semplici di 5.7),
quindi la risposta di A al proprio stimolo diretto potrebbe restare
attiva a lungo da sola (persistenza, sezione 5.6) -- indistinguibile,
guardando solo l'attivita' di A più tardi nel ciclo, da un'eventuale
influenza di C. Lo stesso identico problema gia' affrontato a fondo per
il canale di ritorno B->A nelle sezioni 5.8-5.12.

Soluzione, la stessa gia' validata in diagnosi_ablazione_ba.py (5.12):
ABLAZIONE vera. Due simulazioni IDENTICHE (stesso seed per rumore di
fondo e ogni connettivita' random), che differiscono SOLO per la
presenza o assenza delle sinapsi C->A (gia' allenate al peso finale
osservato in 5.20, per dare al canale ogni vantaggio possibile).
Qualunque differenza nell'attivita' di A durante la parte tarda del
ciclo (escludendo la finestra spiegabile dal solo stimolo diretto) e'
attribuibile ESCLUSIVAMENTE al canale C->A.
"""

from brian2 import *
import numpy as np

N_E = 100
N_I = 25
N_ring = N_E + N_I
N_POP = 3

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms
tau_nmda = 100*ms

sigma_exc = 0.05
CUTOFF_LOCALE = 3 * sigma_exc
w_ee_nmda = 0.325*mV/ms
p_conn = 0.15
w_exc = 1.0*mV
w_inh = 10.0*mV

MESSAGGIO = 0.3
larghezza_stimolo = 0.05
xi = np.arange(N_E) / N_E


def dist_circ(dx):
    dx = abs(dx)
    return np.minimum(dx, 1 - dx)


peso_stim = np.exp(-dist_circ(xi - MESSAGGIO)**2 / (2*larghezza_stimolo**2))
gruppo_stim = np.where(peso_stim > 0.5)[0]

sigma_inter = 0.08
p_inter_picco = 0.5
w_init_ab = 1.5*mV
w_ca_allenato = 2.89*mV  # peso medio finale osservato in 5.20 per C->A --
                           # si parte gia' li' per dare al canale ogni
                           # vantaggio possibile invece di ripartire da 1.5mV

N_RIPETIZIONI = 60
T_STIM = 60*ms
T_PAUSA = 440*ms


def esegui(con_ca, seed_base):
    start_scope()
    seed(seed_base)

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
    A_E.x = xi
    B_E.x = xi
    C_E.x = xi

    background = PoissonGroup(N_POP*N_ring, rates=400*Hz)
    bg_syn = Synapses(background, neurons, on_pre='v_post += 2.0*mV')
    bg_syn.connect(j='i')

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

    syn_ee_A, syn_ei_A, syn_ie_A, syn_ii_A = costruisci_anello(A_E, A_I, seed_locale=0)
    syn_ee_B, syn_ei_B, syn_ie_B, syn_ii_B = costruisci_anello(B_E, B_I, seed_locale=1)
    syn_ee_C, syn_ei_C, syn_ie_C, syn_ii_C = costruisci_anello(C_E, C_I, seed_locale=2)

    tau_stdp = 20*ms
    A_ltp = 0.15*mV
    A_ltd = 0.15*mV
    w_max = 5.0*mV
    w_min = 0*mV

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

    syn_ab = collega(A_E, B_E, seed_conn=10, w_iniziale=w_init_ab)
    syn_bc = collega(B_E, C_E, seed_conn=11, w_iniziale=3.0*mV)

    syn_ca = None
    if con_ca:
        syn_ca = Synapses(C_E, A_E, model='w : volt', on_pre='v_post += w')
        rng_inter = np.random.RandomState(12)
        dmat_inter = dist_circ(xi[:, None] - xi[None, :])
        prob_inter = p_inter_picco * np.exp(-dmat_inter**2 / (2*sigma_inter**2))
        mask_inter = rng_inter.rand(N_E, N_E) < prob_inter
        i_inter, j_inter = np.nonzero(mask_inter)
        syn_ca.connect(i=i_inter, j=j_inter)
        syn_ca.w = w_ca_allenato
    else:
        rng_inter = np.random.RandomState(12)
        _ = rng_inter.rand(N_E, N_E) < p_inter_picco  # consuma l'RNG allo stesso modo

    spikes = SpikeMonitor(neurons)

    for rep in range(N_RIPETIZIONI):
        A_E.I[gruppo_stim] = 15*mV/ms
        run(T_STIM)
        A_E.I[:] = 0*mV/ms
        run(T_PAUSA)

    ciclo_ms = (T_STIM+T_PAUSA)/ms
    mask_a = spikes.i < N_E
    n_tardi_tot, posizioni = 0, []
    for rep in range(N_RIPETIZIONI):
        t0 = rep * ciclo_ms
        fin = t0 + ciclo_ms
        m = mask_a & (spikes.t/ms >= t0 + 20) & (spikes.t/ms < fin) & (spikes.t/ms >= t0 + T_STIM/ms)
        idx = spikes.i[m]
        n_tardi_tot += m.sum()
        posizioni.extend(list(xi[idx]))
    centro = None
    if posizioni:
        centro = float(np.angle(np.mean(np.exp(2j*np.pi*np.array(posizioni)))) / (2*np.pi) % 1.0)
    n_tot_rete = len(spikes.i)
    return n_tardi_tot, centro, n_tot_rete


print("Esecuzione condizione 1/2: CON canale C->A (peso 2.89mV, gia' allenato)...")
n_con, centro_con, tot_con = esegui(con_ca=True, seed_base=42)
print(f"  spike di A dopo la fine dello stimolo diretto (tutti i cicli): {n_con}, centro={centro_con}")
print(f"  spike totali sulla rete: {tot_con}")

print("\nEsecuzione condizione 2/2: SENZA canale C->A (stesso seed, stesso tutto il resto)...")
n_senza, centro_senza, tot_senza = esegui(con_ca=False, seed_base=42)
print(f"  spike di A dopo la fine dello stimolo diretto (tutti i cicli): {n_senza}, centro={centro_senza}")
print(f"  spike totali sulla rete: {tot_senza}")

print("\n" + "="*70)
print("ABLAZIONE: differenza attribuibile ESCLUSIVAMENTE al canale C->A")
print("="*70)
print(f"  CON C->A:    {n_con} spike di A dopo lo stimolo diretto")
print(f"  SENZA C->A:  {n_senza} spike di A dopo lo stimolo diretto")
print(f"  differenza:  {n_con - n_senza:+d} spike")
if n_con > n_senza * 1.2:
    print("  >>> Il canale C->A produce un effetto misurabile: piu' attivita' "
          "di A dopo il proprio stimolo diretto quando il canale e' presente. <<<")
elif abs(n_con - n_senza) <= max(n_con, n_senza) * 0.2:
    print("  >>> Nessuna differenza sostanziale: l'attivita' di A dopo il "
          "proprio stimolo e' spiegata dalla sola persistenza propria, non "
          "dal canale C->A. Il cerchio NON si chiude in modo causale "
          "rilevabile con questa configurazione. <<<")
else:
    print("  >>> Differenza presente ma di segno/entita' ambigua. <<<")
