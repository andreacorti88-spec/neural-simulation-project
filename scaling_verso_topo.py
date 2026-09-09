"""
SCALING VERSO LA SCALA DI UN CERVELLO DI TOPO (~70 milioni di neuroni)
====================================================================
Un cervello di topo ha circa 70 milioni di neuroni. Simularlo per
intero, con dettaglio biofisico, richiede risorse di calcolo che
nessun laptop possiede (serve un cluster o un supercomputer).

Questo script NON pretende di arrivarci. Fa una cosa piu' onesta e
piu' utile: misura, passo dopo passo, qual e' il limite pratico REALE
del tuo Mac -- quanti neuroni riesci a simulare prima che diventi
troppo lento o richieda troppa memoria.

DIFFERENZA IMPORTANTE rispetto agli script precedenti:
Finora usavamo connettivita' a PROBABILITA' fissa (ogni coppia di
neuroni ha una probabilita' p di essere connessa). Questo pero' fa
crescere il numero di sinapsi come N^2 -- con un milione di neuroni,
sarebbe un numero di sinapsi ingestibile.
Qui usiamo connettivita' a GRADO FISSO: ogni neurone si connette
esattamente a K altri neuroni scelti a caso (K=15 di default). Il
numero di sinapsi cresce cosi' in modo LINEARE con N, non quadratico
-- e' cio' che rende possibile scalare a milioni di neuroni.

USO: esegui questo script. Prova' automaticamente dimensioni via via
piu' grandi, stampando tempo e memoria usati. Fermati (Ctrl+C) quando
i tempi diventano troppo lunghi per i tuoi scopi, oppure lascialo
finire: si fermera' da solo se una dimensione fallisce per mancanza
di memoria.
"""

from brian2 import *
import time as pytime
import numpy as np

try:
    import resource
    HAS_RESOURCE = True
except ImportError:
    HAS_RESOURCE = False  # su alcuni sistemi (es. Windows) non e' disponibile


def run_network(N, K_local=15, sim_time=200*ms):
    """Rete E/I bilanciata (80% eccitatori, 20% inibitori) con
    connettivita' a grado fisso K. Ritorna tempo, numero di sinapsi,
    spike totali e (se disponibile) memoria usata."""
    start_scope()

    N_E = int(N * 0.8)
    N_I = N - N_E

    a = 0.02/ms
    b = 0.2/ms
    c = -65*mV
    d = 6*mV/ms

    eqs = '''
    dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
    du/dt = a*(b*v - u) : volt/second
    I : volt/second
    '''

    neurons = NeuronGroup(N, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                           method='euler')
    neurons.v = -65*mV + 15*mV*rand(N) - 5*mV
    neurons.u = b * neurons.v[:] + 2*mV/ms*randn(N)
    neurons.I = 0*mV/ms

    E = neurons[:N_E]
    I_pop = neurons[N_E:]

    # rumore di fondo, come negli script precedenti
    background = PoissonGroup(N, rates=1500*Hz)
    bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
    bg_syn.connect(j='i')

    w_exc = 1.0*mV
    w_inh = 4.0*mV

    K_ei = min(K_local, N_I)
    K_ii = min(K_local, N_I)

    syn_ee = Synapses(E, E, on_pre='v_post += w_exc')
    syn_ee.connect(j=f'k for k in sample(N_E, size={K_local})', skip_if_invalid=True)
    syn_ei = Synapses(E, I_pop, on_pre='v_post += w_exc')
    syn_ei.connect(j=f'k for k in sample(N_I, size={K_ei})', skip_if_invalid=True)
    syn_ie = Synapses(I_pop, E, on_pre='v_post -= w_inh')
    syn_ie.connect(j=f'k for k in sample(N_E, size={K_local})', skip_if_invalid=True)
    syn_ii = Synapses(I_pop, I_pop, on_pre='v_post -= w_inh')
    syn_ii.connect(j=f'k for k in sample(N_I, size={K_ii})', skip_if_invalid=True)

    # record=False: contiamo solo il NUMERO di spike, non salviamo i
    # tempi di ognuno -- risparmia moltissima memoria a grandi scale
    spikes = SpikeMonitor(neurons, record=False)

    t0 = pytime.time()
    run(sim_time)
    elapsed = pytime.time() - t0

    n_synapses = len(syn_ee) + len(syn_ei) + len(syn_ie) + len(syn_ii) + len(bg_syn)
    n_spikes = int(sum(spikes.count))

    mem_mb = None
    if HAS_RESOURCE:
        # su macOS ru_maxrss e' in BYTES (non KB come su Linux)
        mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024*1024)

    return elapsed, n_synapses, n_spikes, mem_mb


# ---------------------------------------------------------------
# SCALING PROGRESSIVO: prova dimensioni via via piu' grandi
# ---------------------------------------------------------------
MOUSE_BRAIN_NEURONS = 70_000_000

sizes_to_test = [1_000, 10_000, 100_000, 500_000,
                  1_000_000, 3_000_000, 5_000_000,
                  10_000_000, 20_000_000]

print(f"{'N neuroni':>12} | {'N sinapsi':>14} | {'tempo (s)':>10} | "
      f"{'RAM (MB)':>10} | {'% cervello topo':>16}")
print("-" * 80)

results = []
for N in sizes_to_test:
    try:
        elapsed, n_syn, n_spk, mem = run_network(N)
        pct_mouse = 100 * N / MOUSE_BRAIN_NEURONS
        mem_str = f"{mem:.0f}" if mem is not None else "n/d"
        print(f"{N:>12,} | {n_syn:>14,} | {elapsed:>10.2f} | "
              f"{mem_str:>10} | {pct_mouse:>15.3f}%")
        results.append((N, elapsed, mem))

        # fermati automaticamente se un singolo run supera i 5 minuti
        if elapsed > 300:
            print("\nTempo superiore a 5 minuti: mi fermo qui per non "
                  "bloccarti il Mac troppo a lungo.")
            break
    except MemoryError:
        print(f"{N:>12,} | MEMORIA ESAURITA -- questo e' il tuo limite pratico.")
        break
    except Exception as e:
        print(f"{N:>12,} | FALLITO ({type(e).__name__}): {e}")
        break

if results:
    max_n_reached = results[-1][0]
    pct = 100 * max_n_reached / MOUSE_BRAIN_NEURONS
    print(f"\n{'='*55}")
    print(f"Dimensione massima raggiunta: {max_n_reached:,} neuroni")
    print(f"Corrisponde al {pct:.3f}% della scala di un cervello di topo "
          f"({MOUSE_BRAIN_NEURONS:,} neuroni)")
    print(f"{'='*55}")
