"""
TOPOLOGIA CONVERGENTE: sommazione di input
====================================================================
Finora avevamo catene lineari (uno dopo l'altro). Qui invece
costruiamo una struttura "molti a uno": 3 neuroni indipendenti
(N0, N1, N2) che convergono TUTTI sullo stesso neurone target (N3).

Questo introduce un fenomeno che con le catene non si vedeva:
la SOMMAZIONE DI INPUT. N3 non riceve stimolo diretto: si attiva
solo se riceve abbastanza spinta combinata dagli altri tre.

Facciamo un confronto diretto:
  - Caso A: i tre neuroni sparano IN FASE (sincronizzati) -> i loro
    effetti su N3 si sommano nello stesso istante -> facile superare soglia
  - Caso B: i tre neuroni sparano SFASATI (in tempi diversi) -> gli
    effetti su N3 arrivano separati -> più difficile superare soglia

E' lo stesso principio per cui, biologicamente, un neurone "conta"
quanti input arrivano vicini nel tempo, non solo quanti in totale.
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

def run_convergent(offsets_ms, synapse_strength_mV, label):
    """N0, N1, N2 sparano una volta ciascuno, con un offset temporale
    diverso specificato in offsets_ms. Tutti convergono su N3."""
    start_scope()

    neurons = NeuronGroup(4, eqs, threshold='v > 30*mV', reset='v = c; u += d',
                           method='euler')
    neurons.v = -65*mV
    neurons.u = b * neurons.v
    neurons.I = 0*mV/ms  # nessuno riceve corrente continua

    # Usiamo SpikeGeneratorGroup per far sparare N0,N1,N2 esattamente
    # ai tempi che vogliamo noi (controllo preciso, utile per confrontare
    # sincrono vs sfasato in modo pulito)
    driver = SpikeGeneratorGroup(3, indices=[0, 1, 2],
                                   times=[t*ms for t in offsets_ms])

    # driver -> i primi 3 neuroni del gruppo "neurons" (fanno da tramite)
    # e poi i primi 3 -> N3 (indice 3), il neurone target
    drive_syn = Synapses(driver, neurons, on_pre='v_post += 35*mV')
    drive_syn.connect(i=[0, 1, 2], j=[0, 1, 2])

    converge_syn = Synapses(neurons, neurons,
                              on_pre=f'v_post += {synapse_strength_mV}*mV')
    converge_syn.connect(i=[0, 1, 2], j=[3, 3, 3])
    converge_syn.delay = 1.5*ms

    spikes = SpikeMonitor(neurons)
    mon = StateMonitor(neurons, 'v', record=3)  # solo N3, il target

    run(50*ms)

    n_target_spikes = sum(spikes.i == 3)
    print(f"{label}: N3 ha sparato {n_target_spikes} volte "
          f"(input a t={offsets_ms} ms)")

    return mon.t/ms, mon.v[0]/mV, spikes.t/ms, spikes.i


# ---------------------------------------------------------------
# CASO A: input sincroni (tutti a t=10ms)
# ---------------------------------------------------------------
t_sync, v_sync, spk_t_sync, spk_i_sync = run_convergent(
    [10, 10, 10], 10, "SINCRONO")

# ---------------------------------------------------------------
# CASO B: input sfasati (a t=10, 15, 20ms)
# ---------------------------------------------------------------
t_async, v_async, spk_t_async, spk_i_async = run_convergent(
    [10, 15, 20], 10, "SFASATO")

# ---------------------------------------------------------------
# VISUALIZZAZIONE COMPARATIVA
# ---------------------------------------------------------------
figure(figsize=(10, 7))

subplot(2, 2, 1)
plot(t_sync, v_sync, color='C0')
axhline(30, color='red', linestyle='--', linewidth=0.8, label='soglia')
title('N3: input SINCRONI (t=10,10,10 ms)')
ylabel('Potenziale (mV)')
legend(fontsize=8)

subplot(2, 2, 2)
plot(t_async, v_async, color='C1')
axhline(30, color='red', linestyle='--', linewidth=0.8, label='soglia')
title('N3: input SFASATI (t=10,15,20 ms)')
legend(fontsize=8)

subplot(2, 2, 3)
for idx in range(4):
    mask = spk_i_sync == idx
    plot(spk_t_sync[mask], spk_i_sync[mask], 'o', markersize=8)
yticks(range(4), ['N0', 'N1', 'N2', 'N3 (target)'])
xlabel('Tempo (ms)')
title('Raster: caso sincrono')

subplot(2, 2, 4)
for idx in range(4):
    mask = spk_i_async == idx
    plot(spk_t_async[mask], spk_i_async[mask], 'o', markersize=8)
yticks(range(4), ['N0', 'N1', 'N2', 'N3 (target)'])
xlabel('Tempo (ms)')
title('Raster: caso sfasato')

suptitle('Sommazione di input: la sincronia conta')
tight_layout()
savefig('rete_convergente.png', dpi=150)
print("\nGrafico salvato: rete_convergente.png")
