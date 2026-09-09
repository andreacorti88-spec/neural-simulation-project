"""
PROVE RIPETUTE E MEDIA (PSTH): la propagazione e' reale o e' rumore?
====================================================================
Nello script precedente, una singola prova sembrava mostrare una
chiara risposta del Modulo 2 allo stimolo del Modulo 1. Ma guardando
il grafico con attenzione, quel rialzo era paragonabile alle normali
fluttuazioni spontanee della rete — non potevamo essere sicuri che
fosse un vero effetto o solo una coincidenza casuale.

LA SOLUZIONE STANDARD in neuroscienza: ripetere l'esperimento molte
volte, allineare ogni prova rispetto al momento dello stimolo, e
fare la MEDIA. Se c'e' un vero effetto, nella media emerge sopra il
rumore (che tende a cancellarsi mediando su piu' prove). Se non c'e'
nessun vero effetto, la media resta piatta. Questa tecnica si chiama
PSTH (Peristimulus Time Histogram) ed e' lo standard per analizzare
dati neurali reali.

In questo script:
  1. Ripetiamo lo stimolo al Modulo 1 per 25 volte, spaziate nel tempo
  2. Per ogni ripetizione, registriamo l'attivita' del Modulo 2 nella
     finestra [-50, +150] ms attorno allo stimolo
  3. Facciamo la media (e calcoliamo l'errore standard) su tutte le
     25 ripetizioni
  4. Confrontiamo il risultato mediato con la variabilita' naturale
     di una singola prova, per capire se l'effetto e' reale
"""

from brian2 import *
import numpy as np

start_scope()

# ---------------------------------------------------------------
# STRUTTURA DELLA RETE (identica allo script dei due moduli)
# ---------------------------------------------------------------
N_E_mod = 120
N_I_mod = 30
N_mod = N_E_mod + N_I_mod
N = N_mod * 2

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

M1_E = neurons[0:120]; M1_I = neurons[120:150]
M2_E = neurons[150:270]; M2_I = neurons[270:300]

background = PoissonGroup(N, rates=1500*Hz)
bg_syn = Synapses(background, neurons, on_pre='v_post += 3*mV')
bg_syn.connect(j='i')

p_local = 0.1
w_exc = 1.0*mV
w_inh = 4.0*mV

def connect_module(E, I):
    see = Synapses(E, E, on_pre='v_post += w_exc')
    see.connect(condition='i!=j', p=p_local)
    sei = Synapses(E, I, on_pre='v_post += w_exc')
    sei.connect(p=p_local)
    sie = Synapses(I, E, on_pre='v_post -= w_inh')
    sie.connect(p=p_local)
    sii = Synapses(I, I, on_pre='v_post -= w_inh')
    sii.connect(condition='i!=j', p=p_local)
    return see, sei, sie, sii

syn1 = connect_module(M1_E, M1_I)
syn2 = connect_module(M2_E, M2_I)

p_inter = 0.02
w_inter = 1.0*mV
syn_12 = Synapses(M1_E, M2_E, on_pre='v_post += w_inter')
syn_12.connect(p=p_inter)
syn_21 = Synapses(M2_E, M1_E, on_pre='v_post += w_inter')
syn_21.connect(p=p_inter)

# ---------------------------------------------------------------
# BURN-IN
# ---------------------------------------------------------------
run(200*ms)
print("Burn-in completato.")

spikes = SpikeMonitor(neurons)

# ---------------------------------------------------------------
# 25 RIPETIZIONI DELLO STIMOLO, spaziate di 300ms l'una dall'altra
# (abbastanza tempo perche' la rete torni alla normalita' tra una
# ripetizione e la successiva)
# ---------------------------------------------------------------
N_TRIALS = 25
ISI = 300*ms          # tempo tra uno stimolo e il successivo
stim_dur = 20*ms
stim_group = M1_E[:20]

stim_times = []
run(50*ms)
for trial in range(N_TRIALS):
    stim_times.append(float(defaultclock.t/ms))
    stim_group.I = 40*mV/ms
    run(stim_dur)
    stim_group.I = 0*mV/ms
    run(ISI - stim_dur)

print(f"Completate {N_TRIALS} ripetizioni dello stimolo.")

# ---------------------------------------------------------------
# COSTRUZIONE DEL PSTH (media sulle prove)
# ---------------------------------------------------------------
t_spikes = np.array(spikes.t/ms)
i_spikes = np.array(spikes.i)
mod2_mask = i_spikes >= 150
mod1_mask = i_spikes < 150

window_pre, window_post, bin_size = 50, 150, 5
bins = np.arange(-window_pre, window_post, bin_size)
bin_centers = bins[:-1] + bin_size/2

def compute_psth(mask, n_neurons):
    all_trials = np.zeros((N_TRIALS, len(bins)-1))
    for k, st in enumerate(stim_times):
        rel_t = t_spikes[mask] - st
        sel = (rel_t >= -window_pre) & (rel_t < window_post)
        hist, _ = np.histogram(rel_t[sel], bins=bins)
        all_trials[k] = hist / (bin_size/1000) / n_neurons  # Hz
    return all_trials

psth1 = compute_psth(mod1_mask, 150)
psth2 = compute_psth(mod2_mask, 150)

mean1, sem1 = psth1.mean(axis=0), psth1.std(axis=0)/np.sqrt(N_TRIALS)
mean2, sem2 = psth2.mean(axis=0), psth2.std(axis=0)/np.sqrt(N_TRIALS)

baseline2 = mean2[bin_centers < 0].mean()
peak2 = mean2[(bin_centers >= 0) & (bin_centers < 50)].max()
trial_std2 = psth2[:, bin_centers < 0].std()  # variabilita' di una singola prova

print("\n" + "=" * 60)
print("RISULTATO (media su 25 prove, con errore standard)")
print("=" * 60)
print(f"Modulo 2 -- baseline media: {baseline2:.2f} Hz")
print(f"Modulo 2 -- picco medio dopo stimolo: {peak2:.2f} Hz")
print(f"Modulo 2 -- effetto medio: +{peak2-baseline2:.2f} Hz")
print(f"Variabilita' naturale di UNA SINGOLA prova (deviazione standard): "
      f"{trial_std2:.2f} Hz")
print(f"\nConfronto: l'effetto medio (+{peak2-baseline2:.2f} Hz) e' "
      f"{'chiaramente sopra' if (peak2-baseline2) > trial_std2 else 'paragonabile a'} "
      f"il rumore di una singola prova ({trial_std2:.2f} Hz)")
print("Questo e' il motivo per cui una sola prova non basta a "
      "concludere se la propagazione e' reale.")

# ---------------------------------------------------------------
# VISUALIZZAZIONE
# ---------------------------------------------------------------
figure(figsize=(11, 7))

subplot(2, 1, 1)
plot(bin_centers, mean1, color='C0', label='Modulo 1 (stimolato)')
fill_between(bin_centers, mean1-sem1, mean1+sem1, color='C0', alpha=0.3)
axvspan(0, 20, color='red', alpha=0.15)
ylabel('Frequenza media (Hz)')
title(f'PSTH Modulo 1 -- media su {N_TRIALS} prove (banda = errore standard)')
legend(loc='upper right', fontsize=8)

subplot(2, 1, 2)
plot(bin_centers, mean2, color='C2', label='Modulo 2 (non stimolato)')
fill_between(bin_centers, mean2-sem2, mean2+sem2, color='C2', alpha=0.3)
axhline(baseline2, color='gray', linestyle='--', linewidth=1, label='baseline')
axvspan(0, 20, color='red', alpha=0.15, label='stimolo')
xlabel('Tempo relativo allo stimolo (ms)')
ylabel('Frequenza media (Hz)')
title(f'PSTH Modulo 2 -- media su {N_TRIALS} prove (banda = errore standard)')
legend(loc='upper right', fontsize=8)

tight_layout()
savefig('psth_moduli.png', dpi=150)
print("\nGrafico salvato in psth_moduli.png")
