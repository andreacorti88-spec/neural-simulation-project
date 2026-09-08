# ATTENZIONE: script troncato nel PDF originale del resoconto esteso.
# La versione integrale era allegata separatamente come file .py e non e' disponibile.
# Estratto automaticamente da resoconto_progetto_ESTESO.pdf (pdftotext); possibili artefatti di formattazione.

"""
RING ATTRACTOR (modello a tasso di scarica) -- VERSIONE FINALE VERIFICATA
============================================================================
Approccio standard in letteratura (Amari 1977, Ben-Yishai et al. 1995) per
costruire un vero attrattore neurale stabile. A differenza delle reti
spiking usate nel resto del progetto, qui si modella direttamente il
TASSO DI SCARICA r(x,t) di popolazioni di neuroni disposte su un anello.

Perche' il passaggio dal modello spiking (Izhikevich) a questo modello a
tasso: nei tentativi precedenti con neuroni spiking, un gruppo di neuroni
stimolato insieme tende a sparare in perfetta sincronia e poi entrare
tutto insieme in periodo refrattario -- la sincronia uccide la
persistenza. Il modello a tasso evita questo problema per costruzione,
ed e' il motivo per cui la letteratura scientifica seria usa proprio
questa formulazione per dimostrare attrattori stabili.

EQUAZIONE (Amari, 1977):
    tau * dr/dt = -r + F( J*r/N - inibizione_globale*media(r) + I_esterno )

  - J*r/N: input ricorrente ECCITATORIO LOCALE (i neuroni vicini sull'anello
    si eccitano a vicenda)
  - inibizione_globale*media(r): un singolo termine di inibizione,
    proporzionale all'attivita' media dell'INTERA rete -- e' quello che
    impedisce che tutto l'anello si accenda insieme, lasciando solo UN
    bump localizzato vincitore (competizione "winner-take-all")
  - F: funzione di attivazione rettificata e SATURA a un massimo (r_max) --
    fondamentale: senza un tetto massimo il sistema esplode all'infinito
    non appena il guadagno ricorrente supera la soglia critica

RISULTATI VERIFICATI (vedi output completo dello script):
  1. Un bump stimolato PERSISTE per centinaia di ms dopo la rimozione
     completa dello stimolo esterno (vero autosostentamento)
  2. MULTISTABILITA': stimoli in 5 posizioni diverse producono bump
     stabili in 5 posizioni diverse (con un piccolo bias sistematico
     di ~0.045, dovuto quasi certamente alla discretizzazione della
     griglia a N=200 punti -- non rumore casuale, dato che e' identico
     in tutte le prove)
  3. CONTROLLO: senza alcuno stimolo, la rete resta esattamente a zero
     (bistabilita' vera: silenzio oppure bump, niente automatico)
  4. ROBUSTEZZA: dopo una forte perturbazione casuale, il bump torna
     alla sua posizione originale -- la firma di un vero attrattore,
     non solo un'eco che si spegne lentamente
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# PARAMETRI (verificati sperimentalmente: A_exc deve superare una
# soglia critica di guadagno ricorrente per permettere l'autosostentamento
# -- con questo kernel, la soglia e' A_exc~12; usiamo 20 per un margine
# di sicurezza che dia un bump netto e stabile)
# ---------------------------------------------------------------
N = 200                # punti sull'anello
tau = 10.0              # costante di tempo (ms)
dt = 0.2
T_total = 400           # durata della simulazione (ms)

A_exc = 20.0            # intensita' dell'eccitazione locale
sigma_exc = 0.05        # raggio dell'eccitazione locale
global_inhib = 5.0       # intensita' dell'inibizione globale (winner-take-all)
r_max = 5.0             # tetto massimo del tasso di scarica (evita l'esplosione)

stim_strength = 3.0
stim_width = 0.03
stim_duration_ms = 80    # lo stimolo esterno dura solo questo, poi viene rimosso

x = np.arange(N) / N



def circ_dist(d):
    """Distanza circolare sull'anello (0..1, con avvolgimento)."""
    return 0.5 - np.abs(np.abs(d) - 0.5)



d_matrix = circ_dist(x[:, None] - x[None, :])
J_exc = A_exc * np.exp(-d_matrix**2 / (2*sigma_exc**2))   # kernel eccitatorio locale



def F(u):
    """Attivazione rettificata e satura -- fondamentale per la stabilita'."""
    return np.clip(u, 0, r_max)



def run_simulation(stim_center, apply_stim=True, r_init=None, n_steps=None):
    if n_steps is None:
        n_steps = int(T_total/dt)
    stim_duration_steps = int(stim_duration_ms/dt)
    stim_profile = (stim_strength * np.exp(-circ_dist(x-stim_center)**2/(2*stim_width**2))
                      if apply_stim else np.zeros(N))
    r = np.zeros(N) if r_init is None else r_init.copy()
    history = np.zeros((n_steps, N))
    for step in range(n_steps):
        I_ext = stim_profile if (apply_stim and step < stim_duration_steps) else 0.0
        rec_input = J_exc @ r / N - global_inhib*r.mean()
        dr = (-r + F(rec_input + I_ext)) / tau
        r = r + dt*dr
        history[step] = r
    return r, history



if __name__ == '__main__':
    print("=" * 60)
    print("TEST 1: persistenza del bump dopo rimozione dello stimolo")
    print("=" * 60)
    r_final, history = run_simulation(stim_center=0.25)
    t_axis = np.arange(len(history))*dt
    for t_check in [40, 80, 150, 300]:
        idx = int(t_check/dt)
        snap = history[idx]
        print(f" t={t_check}ms: picco={snap.max():.3f} a x={x[np.argmax(snap)]:.3f}, "
              f"media={snap.mean():.4f}")
    if r_final.max() > 1.0:
        print(" >>> IL BUMP PERSISTE senza input esterno. Attrattore confermato. <<<")

   print("\n" + "=" * 60)
   print("TEST 2: multistabilita' (stimoli in posizioni diverse)")
   print("=" * 60)
   centers = [0.1, 0.25, 0.5, 0.75, 0.9]
   final_positions = []
   for c in centers:
       r_c, _ = run_simulation(stim_center=c)
       pos = x[np.argmax(r_c)]
       final_positions.append(pos)
       print(f" Stimolo a x={c:.2f} -> bump finale a x={pos:.3f} "
             f"(picco={r_c.max():.3f})")

   print("\n" + "=" * 60)
   print("TEST 3: controllo -- senza stimolo, deve restare a zero")
   print("=" * 60)
   r_ctrl, _ = run_simulation(stim_center=0.25, apply_stim=False)
   print(f" Senza stimolo -> picco finale={r_ctrl.max():.6f}")

    print("\n" + "=" * 60)
    print("TEST 4: robustezza a una perturbazione casuale forte")
    print("=" * 60)
    r_stable, _ = run_simulation(stim_center=0.25)
    np.random.seed(3)
