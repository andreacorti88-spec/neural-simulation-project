"""
SCALING FASE 3: oltre i 25 milioni, con stima di memoria ricalibrata
====================================================================
Nella fase 2 hai raggiunto 25 milioni di neuroni (35.7% della scala
di un cervello di topo) usando solo 11.2 GB di RAM -- molto meno dei
18.7 GB che la stima precedente prevedeva. La mia stima era troppo
prudente. Qui la ricalibro usando ANCHE questo nuovo dato reale, cosi'
la previsione per le dimensioni successive e' piu' accurata.

Due cambiamenti rispetto alla fase 2:
  1. La regressione della memoria include ora il punto a 25M appena
     misurato (6 punti invece di 5)
  2. Il margine di sicurezza scende dal 30% al 15% (i dati reali
     mostrano che la crescita della memoria e' piu' lenta del previsto)
  3. Il limite di tempo per singolo tentativo sale da 5 a 25 minuti
     (hai detto che sei disposto ad aspettare di piu')
"""

from brian2 import *
import time as pytime
import numpy as np

try:
    import resource
    HAS_RESOURCE = True
except ImportError:
    HAS_RESOURCE = False


def run_network(N, K_local=15, sim_time=200*ms):
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

    spikes = SpikeMonitor(neurons, record=False)

    t0 = pytime.time()
    run(sim_time)
    elapsed = pytime.time() - t0

    n_synapses = len(syn_ee) + len(syn_ei) + len(syn_ie) + len(syn_ii) + len(bg_syn)
    n_spikes = int(sum(spikes.count))

    mem_mb = None
    if HAS_RESOURCE:
        mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024*1024)

    return elapsed, n_synapses, n_spikes, mem_mb


# ---------------------------------------------------------------
# DATI OSSERVATI FINORA (fasi 1 e 2, sul tuo Mac) — 6 punti reali
# ---------------------------------------------------------------
N_observed = np.array([1_000_000, 3_000_000, 5_000_000, 10_000_000,
                        20_000_000, 25_000_000])
mem_observed = np.array([1720, 4152, 6254, 9838, 11147, 11495])

slope_mb_per_neuron, intercept = np.polyfit(N_observed, mem_observed, 1)
SAFETY_FACTOR = 1.15   # margine ridotto, visto che i dati reali sono
                        # risultati piu' bassi della stima precedente

def predict_memory_mb(N):
    return SAFETY_FACTOR * (intercept + slope_mb_per_neuron * N)

MOUSE_BRAIN_NEURONS = 70_000_000
TOTAL_RAM_MB = 24 * 1024
SAFETY_MARGIN_MB = 2 * 1024
MAX_SAFE_MEM_MB = TOTAL_RAM_MB - SAFETY_MARGIN_MB
TIME_LIMIT_S = 25 * 60   # 25 minuti per singolo tentativo

sizes_to_test = [30_000_000, 35_000_000, 40_000_000,
                  50_000_000, 60_000_000, 70_000_000]

print(f"{'N neuroni':>12} | {'mem stimata':>12} | {'mem reale':>12} | "
      f"{'tempo (s)':>10} | {'% cervello topo':>16}")
print("-" * 80)

results = []
for N in sizes_to_test:
    predicted = predict_memory_mb(N)
    if predicted > MAX_SAFE_MEM_MB:
        print(f"{N:>12,} | {predicted:>11.0f}MB | -- SALTATO: supererebbe "
              f"la soglia di sicurezza ({MAX_SAFE_MEM_MB:.0f}MB) --")
        print("\nMi fermo qui per non rischiare di bloccare il Mac.")
        break

    try:
        elapsed, n_syn, n_spk, mem = run_network(N)
        pct_mouse = 100 * N / MOUSE_BRAIN_NEURONS
        mem_str = f"{mem:.0f}" if mem is not None else "n/d"
        print(f"{N:>12,} | {predicted:>11.0f}MB | {mem_str:>10}MB | "
              f"{elapsed:>10.2f} | {pct_mouse:>15.3f}%")
        results.append((N, elapsed, mem))

        if elapsed > TIME_LIMIT_S:
            print(f"\nTempo superiore a {TIME_LIMIT_S//60} minuti: mi fermo qui.")
            break
    except MemoryError:
        print(f"{N:>12,} | MEMORIA ESAURITA -- questo e' il tuo limite pratico.")
        break
    except Exception as e:
        print(f"{N:>12,} | FALLITO ({type(e).__name__}): {e}")
        break

if results:
    max_n = results[-1][0]
    pct = 100 * max_n / MOUSE_BRAIN_NEURONS
    print(f"\n{'='*60}")
    print(f"MASSIMO RAGGIUNTO (fase 3): {max_n:,} neuroni")
    print(f"Percentuale della scala di un cervello di topo: {pct:.2f}%")
    print(f"{'='*60}")
else:
    print(f"\nNessuna dimensione testata con successo in questa fase — "
          f"il massimo verificato resta 25 milioni (35.7%).")
