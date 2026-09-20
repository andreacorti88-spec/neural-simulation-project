"""
Ring attractor a tasso di scarica (Amari 1977 / Ben-Yishai 1995):
tau*dr/dt = -r + F(J*r + I_ext), kernel a cappello messicano su anello.
Rate model invece di spiking per evitare i problemi di taratura della
sincronia. Verifica di persistenza del bump dopo rimozione dello stimolo.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(0)

N = 200
x = np.arange(N) / N
tau = 10.0
dt = 0.5
T_total = 600
n_steps = int(T_total/dt)

# kernel a cappello messicano: eccitazione stretta, inibizione piu' ampia
def circ_dist(d):
    return 0.5 - np.abs(np.abs(d) - 0.5)

d_matrix = circ_dist(x[:, None] - x[None, :])  # NxN distanze circolari

sigma_exc = 0.05
sigma_inh = 0.15
A_exc = 3.2
A_inh = 1.0

J = A_exc*np.exp(-d_matrix**2/(2*sigma_exc**2)) - A_inh*np.exp(-d_matrix**2/(2*sigma_inh**2))

def F(u):
    return np.maximum(u, 0)

# stimolo localizzato, poi rimosso -- verifica persistenza del bump
def run_sim(stim_center, stim_duration_steps, stim_strength=2.0, r0=None):
    r = np.zeros(N) if r0 is None else r0.copy()
    history = np.zeros((n_steps, N))
    stim_profile = stim_strength * np.exp(-circ_dist(x-stim_center)**2/(2*0.03**2))

    for step in range(n_steps):
        I_ext = stim_profile if step < stim_duration_steps else 0.0
        rec_input = J @ r / N   # input ricorrente (normalizzato per N)
        dr = (-r + F(rec_input + I_ext)) / tau
        r = r + dt*dr
        history[step] = r
    return history

stim_duration_steps = int(80/dt)  # 80ms, poi rimosso
history = run_sim(stim_center=0.25, stim_duration_steps=stim_duration_steps)

t_after_stim = stim_duration_steps + int(50/dt)  # 50ms dopo fine stimolo
r_after = history[t_after_stim]
r_final = history[-1]
peak_after = r_after.max()
peak_final = r_final.max()
print(f"Picco del bump 50ms dopo la rimozione dello stimolo: {peak_after:.3f}")
print(f"Picco del bump a fine simulazione (t={T_total}ms): {peak_final:.3f}")
print(f"Posizione del picco finale: x={x[np.argmax(r_final)]:.3f} (stimolo era a x=0.25)")

if peak_final > 0.1:
    print(">>> IL BUMP PERSISTE senza input esterno. Attrattore confermato. <<<")
else:
    print(">>> Il bump si e' spento. Serve ritarare i parametri. <<<")

np.save('ring_history.npy', history)
np.save('ring_x.npy', x)
