"""
CLASSIFICAZIONE CON STDP: un neurone che impara a riconoscere un pattern
==========================================================================
Questo e' il passo concettualmente piu' grande finora: un sistema che,
tramite ripetizione e la regola STDP (vista nello script precedente),
impara da solo a distinguere due pattern di input diversi.

STRUTTURA:
  - 6 neuroni di input, divisi in due gruppi:
      Gruppo A = neuroni 0,1,2  (rappresentano il "Pattern A")
      Gruppo B = neuroni 3,4,5  (rappresentano il "Pattern B")
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
stdp_eqs = '''
w : volt
dapre/dt = -apre / tau_stdp : volt (event-driven)
dapost/dt = -apost / tau_stdp : volt (event-driven)
'''
syn = Synapses(driver, out, model=stdp_eqs,
                on_pre='''
                v_post += w
                apre += A_ltp
                w = clip(w + apost, w_min, w_max)
                ''',
                on_post='''
                apost -= A_ltd
                w = clip(w + apre, w_min, w_max)
                ''',
                method='euler')
syn.connect(i=range(6), j=0)
syn.w = 6*mV   # tutti i pesi partono uguali e bassi

# ---------------------------------------------------------------
# MONITORAGGIO E SIMULAZIONE
# ---------------------------------------------------------------
w_mon = StateMonitor(syn, 'w', record=range(6))
spikes = SpikeMonitor(out)
driver_spikes = SpikeMonitor(driver)

run(total_duration*ms)

# ---------------------------------------------------------------
# RISULTATI
# ---------------------------------------------------------------
print("=" * 55)
print("PESI SINAPTICI FINALI (dopo il training)")
print("=" * 55)
for k in range(6):
    grp = "A" if k < 3 else "B"
    print(f"  input {k} (gruppo {grp}): {syn.w[k]/mV:.2f} mV")

spikes_near_A = sum((spikes.t/ms > test_A_time) & (spikes.t/ms < test_A_time+15))
spikes_near_B = sum((spikes.t/ms > test_B_time) & (spikes.t/ms < test_B_time+15))

print("\n" + "=" * 55)
print("RISULTATO DEL TEST (nessun aiuto esterno)")
print("=" * 55)
print(f"Pattern A presentato da solo a t={test_A_time}ms -> "
      f"O ha sparato: {'SI' if spikes_near_A > 0 else 'NO'}")
print(f"Pattern B presentato da solo a t={test_B_time}ms -> "
      f"O ha sparato: {'SI' if spikes_near_B > 0 else 'NO'}")

if spikes_near_A > 0 and spikes_near_B == 0:
    print("\n>>> Il neurone ha imparato a riconoscere selettivamente il Pattern A. <<<")

# ---------------------------------------------------------------
# VISUALIZZAZIONE
# ---------------------------------------------------------------
figure(figsize=(11, 7))

subplot(3, 1, 1)
colors = ['C0', 'C0', 'C0', 'C1', 'C1', 'C1']
for k in range(6):
    plot(w_mon.t/ms, w_mon.w[k]/mV, color=colors[k],
         label=f'Gruppo A' if k == 0 else (f'Gruppo B' if k == 3 else None))
axvline(train_end, color='gray', linestyle=':', linewidth=1)
axvspan(train_end, total_duration, color='yellow', alpha=0.1)
ylabel('Peso sinaptico (mV)')
title('Evoluzione dei pesi durante il training (STDP)')
legend(loc='center left', fontsize=8)

subplot(3, 1, 2)
plot(driver_spikes.t[driver_spikes.i < 3]/ms,
     driver_spikes.i[driver_spikes.i < 3], 'o', color='C0', markersize=4, label='Gruppo A')
plot(driver_spikes.t[driver_spikes.i >= 3]/ms,
     driver_spikes.i[driver_spikes.i >= 3], 'o', color='C1', markersize=4, label='Gruppo B')
axvline(train_end, color='gray', linestyle=':', linewidth=1)
axvspan(train_end, total_duration, color='yellow', alpha=0.1)
yticks(range(6), [f'in{i}' for i in range(6)])
ylabel('Neurone input')
legend(loc='center left', fontsize=8)
title('Pattern presentati (zona gialla = fase di TEST, senza maestro)')

subplot(3, 1, 3)
plot(spikes.t/ms, [0]*len(spikes.t), 'ko', markersize=8)
axvline(train_end, color='gray', linestyle=':', linewidth=1)
axvspan(train_end, total_duration, color='yellow', alpha=0.1)
axvline(test_A_time, color='C0', linestyle='--', linewidth=1, label='test A')
axvline(test_B_time, color='C1', linestyle='--', linewidth=1, label='test B')
yticks([])
xlabel('Tempo (ms)')
title('Spike del neurone di OUTPUT O')
legend(loc='upper right', fontsize=8)

tight_layout()
savefig('classificazione_stdp.png', dpi=150)
print("\nGrafico salvato in classificazione_stdp.png")
