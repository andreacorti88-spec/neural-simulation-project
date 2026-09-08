# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
CLASSIFICAZIONE CON STDP: un neurone che impara a riconoscere un pattern
==========================================================================
Questo e' il passo concettualmente piu' grande finora: un sistema che,
tramite ripetizione e la regola STDP (vista nello script precedente),
impara da solo a distinguere due pattern di input diversi.

STRUTTURA:
  - 6 neuroni di input, divisi in due gruppi:
      Gruppo A = neuroni 0,1,2 (rappresentano il "Pattern A")
      Gruppo B = neuroni 3,4,5 (rappresentano il "Pattern B")
  - Un neurone di OUTPUT (O) riceve da tutti e 6 tramite sinapsi
    plastiche (STDP), pesi iniziali tutti uguali e bassi

FASE DI TRAINING (10 ripetizioni):
  - Ogni volta che presentiamo il Pattern A, diamo anche un piccolo
    aiuto esterno ("maestro") che fa sparare O subito dopo — come un
    insegnante che dice "quando vedi A, rispondi!"
  - Quando presentiamo il Pattern B, NESSUN aiuto: O deve arrangiarsi
  - La regola STDP, ripetizione dopo ripetizione, rafforza le sinapsi
    del Gruppo A (che "predicono" correttamente lo spike di O) e lascia
    indietro (o indebolisce leggermente) quelle del Gruppo B

FASE DI TEST (nessun aiuto esterno):
  - Presentiamo il Pattern A da solo -> se l'apprendimento ha
    funzionato, O deve sparare DA SOLO, senza maestro
  - Presentiamo il Pattern B da solo -> O NON deve sparare

Se questo succede, il sistema ha davvero "imparato" a classificare,
nel senso piu' semplice possibile ma concettualmente reale: ha
sviluppato una risposta selettiva a un pattern specifico, basandosi
solo su ripetizione ed esperienza, non su una regola scritta a mano.
"""

from brian2 import *

# ---------------------------------------------------------------
# PARAMETRI DEL NEURONE (Izhikevich, come negli script precedenti)
# ---------------------------------------------------------------
a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
du/dt = a*(b*v - u) : volt/second
I : volt/second
'''

# ---------------------------------------------------------------
# PARAMETRI STDP
# ---------------------------------------------------------------
tau_stdp = 20*ms
A_ltp = 0.6*mV       # quanto rafforza (potenziamento)
A_ltd = 0.6*mV       # quanto indebolisce (depressione)
w_max = 15*mV        # peso massimo (saturazione)
w_min = 0*mV         # peso minimo

start_scope()

# ---------------------------------------------------------------
# IL NEURONE DI OUTPUT
# ---------------------------------------------------------------
out = NeuronGroup(1, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                   method='euler')
out.v = -65*mV
out.u = b * out.v[:]
out.I = 0*mV/ms

# ---------------------------------------------------------------
# COSTRUZIONE DEL PROTOCOLLO DI STIMOLAZIONE
# 10 ripetizioni di training (A con maestro, B senza), poi test finale
# ---------------------------------------------------------------
indices, times, teacher_times = [], [], []

t = 20
N_TRAIN_REPS = 10
for rep in range(N_TRAIN_REPS):
    for k in [0, 1, 2]:                 # Pattern A
        indices.append(k); times.append(t)
    teacher_times.append(t + 2)         # il "maestro" fa sparare O
    t += 40
    for k in [3, 4, 5]:                 # Pattern B (nessun aiuto)
        indices.append(k); times.append(t)
    t += 40

train_end = t

# Fase di TEST: pattern A da solo, poi pattern B da solo, senza maestro
test_A_time = train_end + 20
test_B_time = train_end + 80
for k in [0, 1, 2]:
    indices.append(k); times.append(test_A_time)
for k in [3, 4, 5]:
    indices.append(k); times.append(test_B_time)

total_duration = test_B_time + 30

driver = SpikeGeneratorGroup(6, indices=indices, times=[x*ms for x in times])
teacher = SpikeGeneratorGroup(1, indices=[0]*len(teacher_times),
                                times=[x*ms for x in teacher_times])

# Il maestro forza uno spike diretto in O (solo durante il training di A)
teacher_syn = Synapses(teacher, out, on_pre='v_post += 50*mV')
teacher_syn.connect()

# ---------------------------------------------------------------
# SINAPSI PLASTICHE INPUT -> OUTPUT (STDP)
# ---------------------------------------------------------------
