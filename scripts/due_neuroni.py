"""
DUE NEURONI COLLEGATI DA UNA SINAPSI
=====================================
Modello: Izhikevich (2003) — scelto perché con poche equazioni riproduce
il comportamento realistico di un neurone (soglia, spike, refrattarietà)
senza la complessità completa di Hodgkin-Huxley. Ottimo punto di partenza.

Cosa succede in questo script:
1. Creiamo il Neurone A e il Neurone B con le stesse proprietà biofisiche.
2. Diamo al Neurone A una corrente di stimolo esterna (lo "svegliamo" noi).
3. Colleghiamo A -> B con una sinapsi eccitatoria: quando A spara,
   B riceve uno stimolo (con un piccolo ritardo, il "delay sinaptico").
4. Registriamo l'attività di entrambi e la plottiamo.

Se B inizia a sparare DOPO che A ha sparato, la comunicazione funziona.
"""

from brian2 import *

# ---------------------------------------------------------------
# 1. PARAMETRI DEL MODELLO IZHIKEVICH
# ---------------------------------------------------------------
# Questi 4 parametri (a, b, c, d) definiscono il "tipo" di neurone.
# I valori qui sotto corrispondono a un neurone "regular spiking" (RS),
# il tipo più comune nella corteccia cerebrale.
a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 6*mV/ms

# Equazioni del modello: v = potenziale di membrana, u = variabile di recovery
# (rappresenta l'attivazione dei canali K+ e l'inattivazione dei canali Na+)
# Nota: qui I è espressa direttamente come termine di "spinta" in mV/ms
# (semplificazione standard usata negli esempi Brian2 per evitare la
# gestione esplicita di capacità/resistenza di membrana).
eqs = '''
dv/dt = (0.04/mV*v**2 + 5*v + 140*mV)/ms - u + I : volt
du/dt = a*(b*v - u) : volt/second
I : volt/second
'''

# ---------------------------------------------------------------
# 2. CREAZIONE DEI DUE NEURONI (un unico gruppo di 2)
# ---------------------------------------------------------------
# threshold='v > 30*mV' : quando il potenziale supera 30 mV, è uno spike
# reset='v = c; u += d'  : dopo lo spike, il potenziale torna a riposo
neurons = NeuronGroup(2, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                       method='euler')
neurons.v = -65*mV   # entrambi partono a riposo
neurons.u = b * neurons.v

# Solo il Neurone A (indice 0) riceve stimolo esterno.
# Il Neurone B (indice 1) NON riceve nulla di suo: si attiva SOLO
# se la sinapsi da A lo eccita abbastanza.
neurons.I = [15*mV/ms, 0*mV/ms]

# ---------------------------------------------------------------
# 3. LA SINAPSI: A -> B
# ---------------------------------------------------------------
# on_pre = cosa succede a B quando A spara ("pre" = il neurone che manda)
# Aggiungiamo un impulso di corrente a B, con un ritardo sinaptico
# realistico di 1-2 ms (il tempo che ci mette il neurotrasmettitore
# a diffondere e aprire i canali nel neurone ricevente).
syn = Synapses(neurons, neurons, on_pre='v_post += 20*mV')
syn.connect(i=0, j=1)   # collega neurone 0 (A) -> neurone 1 (B)
syn.delay = 1.5*ms

# ---------------------------------------------------------------
# 4. MONITORAGGIO E SIMULAZIONE
# ---------------------------------------------------------------
mon = StateMonitor(neurons, 'v', record=True)
spikes = SpikeMonitor(neurons)

run(200*ms)

# ---------------------------------------------------------------
# 5. VISUALIZZAZIONE
# ---------------------------------------------------------------
figure(figsize=(10, 5))

subplot(2, 1, 1)
plot(mon.t/ms, mon.v[0]/mV, label='Neurone A (stimolato)')
plot(mon.t/ms, mon.v[1]/mV, label='Neurone B (riceve da A)')
ylabel('Potenziale (mV)')
legend(loc='upper right')
title('Potenziale di membrana: comunicazione A -> B')

subplot(2, 1, 2)
plot(spikes.t/ms, spikes.i, 'o', markersize=8)
yticks([0, 1], ['A', 'B'])
xlabel('Tempo (ms)')
ylabel('Neurone')
title('Raster degli spike (i puntini = quando ogni neurone spara)')

tight_layout()
savefig('due_neuroni_output.png', dpi=150)
print("Simulazione completata. Grafico salvato in due_neuroni_output.png")
print(f"Numero di spike del Neurone A: {sum(spikes.i == 0)}")
print(f"Numero di spike del Neurone B: {sum(spikes.i == 1)}")
