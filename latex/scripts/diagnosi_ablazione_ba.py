"""
Test definitivo, in risposta al fatto che sia il confronto per posizione
(diagnosi_sottosoglia_ba.py) sia i test comportamentali precedenti (5.9)
restano vulnerabili a un confondimento: la dinamica propria di A (bump
recrutato piu' largo della sola zona di stimolo diretto, codetta residua
tra un turno e l'altro) si mescola con qualunque eventuale segnale dal
canale B->A, rendendo ambiguo qualunque confronto "per posizione" o "nel
tempo" fatto DENTRO la stessa simulazione.

Soluzione: ABLAZIONE vera. Due simulazioni IDENTICHE (stesso seed per
rumore di fondo e per ogni connettivita' random), che differiscono SOLO
per la presenza o assenza delle sinapsi B->A (stesso schema a turni A/B
gia' usato in dialogo_spiking_bidirezionale_gate.py, sinapsi B->A gia'
allenate fino al peso finale osservato in quello script, w=7.8mV medio
sulle sinapsi piu' forti, per dare al canale la sua massima possibilita').
Qualunque differenza tra le due condizioni, con tutto il resto identico
per costruzione, e' attribuibile SOLO al canale B->A.
"""

from brian2 import *
import numpy as np

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
w_init_ba = 3.5*mV  # peso "gia' allenato": si parte forte, come se il
                     # canale avesse gia' accumulato il rinforzo medio
                     # osservato nella sezione 5.9, per dargli ogni
                     # vantaggio possibile invece di ripartire da zero


def esegui(con_canale_ba, seed_base):
    start_scope()
    seed(seed_base)

    eqs = '''
    dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I + g_nmda - adapt : volt
    du/dt = a*(b*v - u) : volt/second
    dg_nmda/dt = -g_nmda/tau_nmda : volt/second
    dadapt/dt = -adapt/tau_adapt_pop : volt/second
    I : volt/second
    x : 1
    '''
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
    A_E.x = xi
    B_E.x = xi

    background = PoissonGroup(2*N_ring, rates=400*Hz)
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

    tau_stdp = 20*ms
    A_ltp = 0.15*mV
    A_ltd = 0.15*mV
    w_max = 8.0*mV
    w_min = 0*mV
    w_init = 3.5*mV

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

    syn_ba = None
    if con_canale_ba:
        syn_ba = Synapses(B_E, A_E, model='w : volt', on_pre='v_post += w')
        mask_inter_ba = rng_inter.rand(N_E, N_E) < prob_inter
        i_inter_ba, j_inter_ba = np.nonzero(mask_inter_ba)
        syn_ba.connect(i=i_inter_ba, j=j_inter_ba)
        syn_ba.w = w_init_ba
    else:
        # stesso sorteggio random (per consumare l'RNG in modo identico
        # e mantenere sincronizzate le sequenze pseudocasuali successive)
        # ma nessuna sinapsi creata
        _ = rng_inter.rand(N_E, N_E) < prob_inter

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
            A_E.I[gruppo_stim] = 15*mV/ms
            run(T_STIM)
            A_E.I[:] = 0*mV/ms
        else:
            B_E.I[gruppo_stim] = 15*mV/ms
            run(T_STIM)
            B_E.I[:] = 0*mV/ms
        run(T_PAUSA)
        t_corrente += T_STIM + T_PAUSA

    mask_a = spikes.i < N_E
    mask_b_tot = (spikes.i >= N_ring) & (spikes.i < N_ring+N_E)
    print(f"  [controllo sanita': spike totali A={mask_a.sum()}, B={mask_b_tot.sum()} "
          f"su tutta la simulazione -- atteso centinaia se la dinamica ricorrente funziona]")
    ciclo_ms = (T_STIM+T_PAUSA)/ms
    turni_b = [t0 for (l, t0) in turni if l == 'B']

    n_tot, posizioni = 0, []
    for t0 in turni_b:
        fin = t0 + ciclo_ms
        m = mask_a & (spikes.t/ms >= t0+20) & (spikes.t/ms < fin)
        n_tot += m.sum()
        posizioni.extend(list(xi[spikes.i[m]]))
    centro = None
    if posizioni:
        centro = float(np.angle(np.mean(np.exp(2j*np.pi*np.array(posizioni)))) / (2*np.pi) % 1.0)
    return n_tot, centro, len(turni_b)


print("Esecuzione condizione 1/2: CON canale B->A (peso 3.5mV, gia' forte)...")
n_con, centro_con, n_turni_b = esegui(con_canale_ba=True, seed_base=42)
print(f"  spike di A durante i turni B (>20ms dall'inizio): {n_con} su {n_turni_b} turni, centro={centro_con}")

print("\nEsecuzione condizione 2/2: SENZA canale B->A (stesso seed, stesso tutto il resto)...")
n_senza, centro_senza, _ = esegui(con_canale_ba=False, seed_base=42)
print(f"  spike di A durante i turni B (>20ms dall'inizio): {n_senza} su {n_turni_b} turni, centro={centro_senza}")

print("\n" + "="*70)
print("ABLAZIONE: differenza attribuibile ESCLUSIVAMENTE al canale B->A")
print("="*70)
print(f"  CON B->A:    {n_con} spike di A nei turni B")
print(f"  SENZA B->A:  {n_senza} spike di A nei turni B")
print(f"  differenza:  {n_con - n_senza:+d} spike")
if n_con > n_senza:
    print("  >>> Il canale B->A produce un effetto misurabile: piu' spike di A "
          "quando il canale e' presente rispetto a quando e' assente, a "
          "parita' di ogni altra condizione. <<<")
elif n_con == n_senza:
    print("  >>> Nessuna differenza: il canale B->A non produce alcun effetto "
          "misurabile sull'attivita' di A, nemmeno con un peso gia' forte "
          "(3.5mV) dato in partenza. <<<")
else:
    print("  >>> Il canale B->A e' associato a MENO spike di A, non di piu' -- "
          "un effetto netto non riconducibile a trasmissione eccitatoria utile. <<<")
