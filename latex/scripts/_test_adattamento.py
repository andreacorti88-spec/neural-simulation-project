"""
Primo mattone verso il dominio spiking per la comunicazione bidirezionale
appresa (sezione 5.5): un vero ring attractor a impulsi, invece del modello
a rate usato finora (ring_attractor_rate.py e tutta la famiglia dialogo_*).
ring_attractor_rate.py aveva scelto deliberatamente il rate model "per
evitare i problemi di taratura della sincronia" -- questo script affronta
direttamente quel problema, verificando se un vero attrattore a impulsi
regge prima di tentare qualunque comunicazione tra due di essi.

Architettura: popolazione eccitatoria E disposta su un anello (la posizione
del neurone i-esimo e' i/N_E), con connettivita' E->E a caduta gaussiana
sulla distanza circolare (stessa larghezza sigma_exc=0.05 del modello a
rate) per dare localita' spaziale al bump. Popolazione inibitoria I
condivisa e non spazializzata (connettivita' casuale sparsa, stesso stile
di rete_300_neuroni.py) per l'inibizione globale che sopprime l'attivita'
lontano dal bump -- l'analogo spiking del termine "-global_inhib*r.mean()"
del modello a rate. Stesso modello di neurone (Izhikevich RS) usato in
tutto il resto del progetto.

Test: uno stimolo localizzato e transitorio innesca il bump; si verifica
se persiste dopo la rimozione dello stimolo (stesso test di
ring_attractor_rate.py) e se resta spazialmente localizzato invece di
diffondersi o spegnersi.
"""

from brian2 import *

start_scope()

N_E = 100
N_I = 25
N = N_E + N_I

a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

# tentativi 1-5: sinapsi E->E istantanee (v_post += peso) -- o non bastavano
# a formare un bump, o lo facevano propagare a tutto l'anello come un'onda
# viaggiante invece di restare fermo. Qui la componente E->E diventa LENTA
# (costante di tempo tipo NMDA, ~100ms, Compte et al. 2000 / Wang 2002):
# integra molti spike nel tempo invece di dare uno scatto istantaneo,
# smorzando le instabilita' veloci che innescano le onde
tau_nmda = 100*ms

# adattamento di popolazione, analogo alla variabile aggiunta nel modello a
# rate (dialogo_stdp_adattamento.py): una corrente lenta che cresce ad ogni
# spike e si sottrae dalla spinta della cellula, cosi' che un bump sostenuto
# a lungo si affievolisca da solo -- la variabile "u" di Izhikevich gia'
# presente non basta, e' troppo debole a livello di rete per far spegnere
# un bump sostenuto dalla ricorrenza (verificato: persistenza per 900ms in
# ring_attractor_spiking.py senza alcun segno di affievolimento)
tau_adapt_pop = 200*ms
gain_adapt_pop = 0.15*mV/ms  # tentativo 1: 0.4 spegneva il bump quasi subito
                              # (16 spike appena dopo lo stimolo, poi nulla)

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I + g_nmda - adapt : volt
du/dt = a*(b*v - u) : volt/second
dg_nmda/dt = -g_nmda/tau_nmda : volt/second
dadapt/dt = -adapt/tau_adapt_pop : volt/second
I : volt/second
x : 1
'''

neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d; adapt += gain_adapt_pop',
                       method='euler')
neurons.v = -65*mV + 5*mV*randn(N)
neurons.u = b * neurons.v[:]
neurons.I = 0*mV/ms
neurons.g_nmda = 0*mV/ms
neurons.adapt = 0*mV/ms

E = neurons[:N_E]
I_pop = neurons[N_E:]

# posizione sull'anello, solo per E (I non e' spazializzata)
E.x = np.arange(N_E) / N_E

# rumore di fondo: 1500Hz/3mV (tentativo 1) dominava la dinamica; 400Hz/2mV
# (tentativi 2-4) era troppo debole per sostenere un'attivita' di base --
# via di mezzo qui, cosi' l'eccitazione ricorrente locale amplifica
# un'attivita' di base gia' presente invece di dover accendere tutto da zero
background = PoissonGroup(N, rates=400*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 2.0*mV')
bg_syn.connect(j='i')

sigma_exc = 0.05
CUTOFF_LOCALE = 3 * sigma_exc  # nessuna connessione E->E oltre questa distanza
                                # circolare -- i tentativi precedenti usavano
                                # una probabilita' gaussiana SENZA taglio netto,
                                # permettendo connessioni sporadiche a lungo
                                # raggio che facevano propagare l'attivita'
                                # a tutto l'anello invece di restare locale
w_ee_peak = 4.0*mV  # tentativo 1: 1.2mV (nessun bump); tentativo 2: 4.0mV senza
                     # taglio (localizzato ma non autosostenuto); tentativo 3:
                     # 7.0mV (propagazione a tutto l'anello, poi silenzio
                     # totale); tentativo 4: 5.0mV+inibizione piu' forte
                     # (ancora propagazione, solo piu' debole); tentativo 5:
                     # con taglio netto ma sinapsi istantanee (onde viaggianti
                     # persistenti, ne' ferme ne' spente)
# tentativo 6: la stessa w_ee_peak, ma come incremento a g_nmda (corrente
# lenta) invece che salto istantaneo di v -- calibrato piu' basso perche'
# si accumula nel tempo invece di agire in un colpo solo
w_ee_nmda = 0.325*mV/ms  # tentativo 8: 0.35 cresceva senza limite; tentativo
                          # 9: 0.30 si spegneva tra le due finestre (ma bump
                          # stretto e ben posizionato finche' dura) -- via di
                          # mezzo qui, il punto di equilibrio dovrebbe stare
                          # tra i due
p_conn = 0.15  # densita' I aumentata rispetto a 0.1, per un'inibizione piu' pronta
w_exc = 1.0*mV
w_inh = 10.0*mV  # tentativo 6: 6.0mV (diffusione); tentativo 7: 12.0mV
                  # (bump localizzato che resta fermo, ma si affievolisce
                  # gradualmente invece di restare stabile)


def dist_circ(dx):
    dx = abs(dx)
    return np.minimum(dx, 1 - dx)


# E->E: caduta gaussiana sulla distanza circolare, connessione esplicita
# (non tramite la sintassi probabilistica standard di Brian2, per poter
# calcolare il peso in funzione della distanza -- stesso approccio usato
# per la connettivita' spaziale in rete_smallworld.py)
xi = np.arange(N_E) / N_E
dmat = dist_circ(xi[:, None] - xi[None, :])
prob_ee = np.exp(-dmat**2 / (2*sigma_exc**2))
prob_ee[dmat > CUTOFF_LOCALE] = 0.0  # taglio netto, nessuna coda a lungo raggio
np.fill_diagonal(prob_ee, 0)
rng = np.random.RandomState(0)
mask_ee = rng.rand(N_E, N_E) < prob_ee
i_idx, j_idx = np.nonzero(mask_ee)

syn_ee = Synapses(E, E, on_pre='g_nmda_post += w_ee_nmda')
syn_ee.connect(i=i_idx, j=j_idx)
print(f"Sinapsi E->E create: {len(i_idx)} (media {len(i_idx)/N_E:.1f} per neurone)")

syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
syn_ei.connect(p=p_conn)
syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
syn_ie.connect(p=p_conn)
syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
syn_ii.connect(condition='i!=j', p=p_conn)

spikes = SpikeMonitor(neurons)
rate_mon = PopulationRateMonitor(E)

run(200*ms)
print("Baseline (200ms) completata.")

# stimolo: impulso di corrente su un gruppo localizzato di neuroni E
# intorno alla posizione 0.3 sull'anello, per 60ms
centro_stimolo = 0.3
larghezza_stimolo = 0.05
peso_stim = np.exp(-dist_circ(xi - centro_stimolo)**2 / (2*larghezza_stimolo**2))
gruppo_stim = np.where(peso_stim > 0.5)[0]
print(f"Neuroni stimolati (vicino a x={centro_stimolo}): {len(gruppo_stim)}")

E.I[gruppo_stim] = 15*mV/ms
run(60*ms)
E.I[:] = 0*mV/ms
print("Stimolo rimosso, verifica persistenza per altri 800ms (finestra piu' "
      "lunga per distinguere un vero equilibrio da una lenta deriva)...")

run(800*ms)

fig, axes = subplots(2, 1, figsize=(10, 7))

subplot(2, 1, 1)
mask_e = spikes.i < N_E
plot(spikes.t[mask_e]/ms, spikes.i[mask_e], '.', color='C0', markersize=2, label=f'E (n={N_E})')
plot(spikes.t[~mask_e]/ms, spikes.i[~mask_e]-N_E, '.', color='C3', markersize=2, label=f'I (n={N_I}, indice traslato)')
axvspan(200, 260, color='red', alpha=0.15, label='stimolo')
xlabel('Tempo (ms)')
ylabel('Indice neurone (E) / I traslato')
title(f'Raster: ring attractor spiking, stimolo a x={centro_stimolo} rimosso dopo 60ms')
legend(loc='upper right', fontsize=8, markerscale=4)

subplot(2, 1, 2)
plot(rate_mon.t/ms, rate_mon.smooth_rate(window='flat', width=5*ms)/Hz, color='black')
axvspan(200, 260, color='red', alpha=0.15)
xlabel('Tempo (ms)')
ylabel('Frequenza popolazione E (Hz)')
title('Attivita\' collettiva di E nel tempo')
grid(alpha=0.3)

tight_layout()
savefig('ring_attractor_spiking.png', dpi=150)
print("\nGrafico salvato in ring_attractor_spiking.png")

# verifica quantitativa di persistenza, localizzazione, E STABILITA' (due
# finestre distanziate: se un vero equilibrio, dovrebbero assomigliarsi; se
# deriva/crescita lenta, la seconda finestra sara' piu' attiva e piu' larga)
spikes_e = spikes.i[mask_e]
tempi_e = spikes.t[mask_e]


def analizza_finestra(t0, t1, etichetta):
    m = (tempi_e > t0*ms) & (tempi_e < t1*ms)
    n = m.sum()
    if n == 0:
        print(f"  {etichetta} ({t0}-{t1}ms): NESSUNO spike -- il bump si e' spento.")
        return None, 0, None
    posizioni_attive = xi[spikes_e[m]]
    centro = np.angle(np.mean(np.exp(2j*np.pi*posizioni_attive))) / (2*np.pi) % 1.0
    larghezza = np.std([dist_circ(p - centro) for p in posizioni_attive])
    print(f"  {etichetta} ({t0}-{t1}ms): {n} spike, centro x={centro:.3f}, "
          f"larghezza (dev.std. posizioni) = {larghezza:.3f}")
    return centro, n, larghezza


print()
c0, n0, l0 = analizza_finestra(260, 360, "Subito dopo stimolo")
c1, n1, l1 = analizza_finestra(460, 660, "Finestra precoce")
c2, n2, l2 = analizza_finestra(900, 1100, "Finestra tardiva")

if n1 > 20 and n2 > 20:
    crescita = (n2 - n1) / n1 * 100
    print(f"\n  Variazione di attivita' tra le due finestre: {crescita:+.0f}%")
    if l1 is not None and l2 is not None:
        print(f"  Variazione di larghezza: {l1:.3f} -> {l2:.3f}")
    if abs(crescita) < 30 and l2 < l1 * 1.5:
        print("  >>> Il bump sembra vicino a un vero equilibrio: localizzato, "
              "persistente, ne' in crescita ne' in espansione marcata. <<<")
    else:
        print("  >>> Il bump persiste ma non e' ancora stabile: sta ancora "
              "crescendo o allargandosi tra le due finestre. <<<")
elif n1 > 20 and n2 <= 20:
    print("\n  >>> Il bump si e' spento tra le due finestre: persistenza insufficiente. <<<")
elif n1 <= 20:
    print("\n  >>> Attivita' gia' troppo scarsa nella finestra precoce. <<<")
