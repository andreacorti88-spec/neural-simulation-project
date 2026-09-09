"""
STDP: Spike-Timing-Dependent Plasticity
====================================================================
Questa e' la regola di apprendimento sinaptico piu' usata in
neuroscienza computazionale, ed e' molto piu' realistica della
semplice "facilitazione" vista prima.

Principio (regola di Hebb, versione temporale):
  - Se A spara POCO PRIMA di B (es. 5ms prima)  -> la sinapsi A->B
    si RAFFORZA (potenziamento, "LTP": Long-Term Potentiation)
  - Se A spara POCO DOPO B (es. 5ms dopo)        -> la sinapsi A->B
    si INDEBOLISCE (depressione, "LTD": Long-Term Depression)
  - Piu' i due spike sono vicini nel tempo, piu' l'effetto e' forte
    (decade esponenzialmente con la distanza temporale)

Questo e' esattamente il meccanismo biologico con cui il cervello
"impara" relazioni causa-effetto tra eventi/neuroni: se A predice
sistematicamente B, la connessione si rafforza nel tempo.

In questo script:
  1. Stimoliamo DUE VOLTE la coppia A-B, sempre con A che spara
     circa 5ms prima di B (relazione causale ripetuta)
  2. Osserviamo il peso sinaptico salire ad ogni ripetizione
  3. Poi facciamo l'esperimento opposto: B prima di A -> il peso scende
"""

from brian2 import *

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
tau_stdp = 20*ms          # finestra temporale entro cui l'STDP ha effetto
A_ltp = 1.0*mV            # quanto si rafforza se A precede B (potenziamento)
A_ltd = 1.0*mV            # quanto si indebolisce se B precede A (depressione)
w_max = 40*mV             # peso sinaptico massimo (satura, come nella realta')
w_min = 0*mV              # peso sinaptico minimo (non puo' diventare negativo)

# ---------------------------------------------------------------
# COSTRUZIONE DEL MODELLO
# Usiamo due variabili di traccia (apre, apost) che si accumulano
# ad ogni spike e decadono nel tempo: e' il modo standard di
# implementare STDP in Brian2 (traccia esponenziale pre e post-sinaptica)
# ---------------------------------------------------------------
start_scope()

neurons = NeuronGroup(2, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV
neurons.u = b * neurons.v
neurons.I = 0*mV/ms

# Usiamo un SpikeGeneratorGroup per controllare esattamente QUANDO
# sparano A e B nei due esperimenti, cosi' vediamo l'effetto pulito
# senza dipendere dalla dinamica naturale del neurone.
# Esperimento: 3 ripetizioni di "A a t, B a t+5ms" (A causa B)
pattern_causale = []
pattern_indices = []
for rep in range(5):
    t0 = 20 + rep*30
    pattern_causale += [t0, t0+5]
    pattern_indices += [0, 1]   # 0 = A, 1 = B

driver = SpikeGeneratorGroup(2, indices=pattern_indices,
                               times=[t*ms for t in pattern_causale])

# driver spara direttamente A e B ai tempi voluti (forziamo lo spike
# dando un impulso fortissimo, così partiamo da tempi certi)
force_syn = Synapses(driver, neurons, on_pre='v_post += 50*mV')
force_syn.connect(j='i')

# ---------------------------------------------------------------
# LA SINAPSI PLASTICA A -> B, CON REGOLA STDP
# ---------------------------------------------------------------
stdp_eqs = '''
w : volt
dapre/dt = -apre / tau_stdp : volt (event-driven)
dapost/dt = -apost / tau_stdp : volt (event-driven)
'''

syn = Synapses(neurons, neurons, model=stdp_eqs,
                on_pre='''
                apre += A_ltp
                w = clip(w + apost, w_min, w_max)
                ''',
                on_post='''
                apost -= A_ltd
                w = clip(w + apre, w_min, w_max)
                ''',
                method='euler')
syn.connect(i=0, j=1)
syn.w = 5*mV   # peso iniziale, deliberatamente basso
initial_w = syn.w[0]   # salviamo il valore iniziale PRIMA di far girare la simulazione

w_mon = StateMonitor(syn, 'w', record=0)
spikes = SpikeMonitor(neurons)

run(200*ms)

# ---------------------------------------------------------------
# VISUALIZZAZIONE
# ---------------------------------------------------------------
figure(figsize=(10, 6))

subplot(2, 1, 1)
plot(w_mon.t/ms, w_mon.w[0]/mV, color='purple', linewidth=2)
ylabel('Peso sinaptico w (mV)')
title('STDP: il peso cresce ad ogni ripetizione "A prima di B" (causale)')
grid(alpha=0.3)

subplot(2, 1, 2)
plot(spikes.t[spikes.i == 0]/ms, [0]*sum(spikes.i == 0), 'o', color='C0', label='A')
plot(spikes.t[spikes.i == 1]/ms, [1]*sum(spikes.i == 1), 'o', color='C1', label='B')
yticks([0, 1], ['A', 'B'])
xlabel('Tempo (ms)')
title('Pattern di stimolazione: A spara sempre 5ms prima di B')
legend(loc='upper right', fontsize=8)

tight_layout()
savefig('stdp.png', dpi=150)

print("Simulazione completata. Grafico salvato in stdp.png")
print(f"Peso sinaptico iniziale: {initial_w/mV:.2f} mV")
print(f"Peso sinaptico finale:   {w_mon.w[0][-1]/mV:.2f} mV")
print(f"Rafforzamento totale dopo 5 ripetizioni 'A causa B': "
      f"+{(w_mon.w[0][-1] - initial_w)/mV:.2f} mV")
print("\nProva ora a invertire l'ordine nel pattern (B prima di A)")
print("nel codice sopra, per vedere il peso SCENDERE invece di salire.")
