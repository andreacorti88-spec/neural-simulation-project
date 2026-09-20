"""
Verifica diretta dell'ipotesi proposta in dialogo_multi_repulsione.py
(sezione 5.5.6 del resoconto): l'interferenza tra due messaggi osservata
in dialogo_multi_messaggio.py non nasceva da una vera competizione per
risorse condivise, ma dal fatto che la ricerca esplorativa casuale
riscopriva per puro caso la STESSA zona "facile" dell'anello per entrambi
-- lo stesso bias sistematico di posizione gia' documentato (mai spiegato
fino in fondo) per il ring attractor nella sezione 9 del resoconto.

Se l'ipotesi e' corretta, allenando molti messaggi diversi in modo
COMPLETAMENTE INDIPENDENTE (una rete fresca per ciascuno, nessuna
condivisione di pesi, nessuna repulsione) le posizioni di risposta finali
dovrebbero mostrare una tendenza a raggrupparsi in poche zone ricorrenti
invece di distribuirsi in modo uniforme sull'anello. Se invece l'ipotesi
e' sbagliata (l'interferenza era davvero competizione per risorse
condivise, non un artefatto della ricerca), le posizioni indipendenti
dovrebbero risultare sparse senza pattern riconoscibile.

Nessuna sinapsi condivisa qui: ogni messaggio allena la propria coppia
W_AB/W_BA da zero, con il proprio seed casuale indipendente -- l'unica
cosa in comune tra le run e' la fisica del ring attractor stesso.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from dialogo_stdp_adattamento import (
    init_weak_bias, run_episode, eval_batch, circ_dist,
    OFFSET_AB_INNATO, OFFSET_BA_INNATO, W_INIT_SCALE, DECAY_PER_EPISODE,
)

MESSAGGI_TEST = [0.05, 0.20, 0.35, 0.50, 0.65, 0.80]
N_EPISODES = 1800  # sufficiente perche' A->B raggiunga l'autosufficienza
                    # (verificato nei test di regressione delle sezioni precedenti)

if __name__ == '__main__':
    print("=" * 70)
    print(f"DIAGNOSI: {len(MESSAGGI_TEST)} messaggi allenati in modo INDIPENDENTE "
          f"(reti separate, nessun peso condiviso)")
    print("=" * 70)

    posizioni_target = []
    for k, messaggio in enumerate(MESSAGGI_TEST):
        W_AB = init_weak_bias(OFFSET_AB_INNATO, W_INIT_SCALE)
        W_BA = init_weak_bias(OFFSET_BA_INNATO, W_INIT_SCALE)
        rng = np.random.RandomState(1000 + k)  # seed diverso per ciascun messaggio

        for ep in range(N_EPISODES):
            run_episode(W_AB, W_BA, messaggio, rng, train=True, allena_BA=False)
            W_AB *= DECAY_PER_EPISODE

        r = eval_batch(W_AB, W_BA, messaggio, n_trial=10, seed0=5000 + k)
        target = r['media_B']
        offset = circ_dist(target - messaggio) if target is not None else None
        posizioni_target.append(target)
        print(f"  messaggio x={messaggio:.2f}  ->  target B={target if target is None else f'{target:.3f}'}"
              f"  (offset={offset if offset is None else f'{offset:.3f}'})  "
              f"tasso B(t1)={r['tasso_B1']*100:.0f}%")

    print("\n" + "=" * 70)
    print("VERIFICA CLUSTERING: le posizioni target sono raggruppate o sparse?")
    print("=" * 70)
    validi = [p for p in posizioni_target if p is not None]
    if len(validi) >= 2:
        distanze = []
        for i in range(len(validi)):
            for j in range(i+1, len(validi)):
                distanze.append(abs(circ_dist(validi[i]-validi[j])))
        distanza_media = np.mean(distanze)
        distanza_attesa_uniforme = 0.26  # valore atteso per punti uniformi su un anello [0,1)
        print(f"  Posizioni target: {[f'{p:.3f}' for p in validi]}")
        print(f"  Distanza circolare media tra coppie: {distanza_media:.3f}")
        print(f"  Distanza attesa se le posizioni fossero uniformi e indipendenti: "
              f"~{distanza_attesa_uniforme:.2f}")
        if distanza_media < distanza_attesa_uniforme * 0.6:
            print("  >>> Le posizioni sono PIU' VICINE tra loro di quanto ci si aspetterebbe "
                  "per caso: indizio di un bias sistematico verso poche zone ricorrenti "
                  "dell'anello, coerente con l'ipotesi. <<<")
        elif distanza_media > distanza_attesa_uniforme * 1.4:
            print("  >>> Le posizioni sono PIU' SPARSE di quanto ci si aspetterebbe per caso. <<<")
        else:
            print("  >>> Le posizioni sono compatibili con una distribuzione uniforme -- "
                  "nessuna evidenza chiara di clustering in questo campione. <<<")

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={'projection': 'polar'})
    for messaggio, target in zip(MESSAGGI_TEST, posizioni_target):
        ang_m = messaggio * 2*np.pi
        ax.plot([ang_m], [1.0], 'o', color='C0', markersize=10)
        ax.text(ang_m, 1.12, f'{messaggio:.2f}', ha='center', fontsize=8, color='C0')
        if target is not None:
            ang_t = target * 2*np.pi
            ax.plot([ang_t], [0.6], 's', color='C3', markersize=10)
            ax.plot([ang_m, ang_t], [1.0, 0.6], '-', color='gray', alpha=0.4, linewidth=1)
    ax.plot([], [], 'o', color='C0', label='messaggio (A)')
    ax.plot([], [], 's', color='C3', label='target appreso (B)')
    ax.set_yticklabels([])
    ax.set_title("Messaggi indipendenti: posizione del messaggio (blu)\ne del target appreso (rosso) sull'anello")
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=8)

    plt.tight_layout()
    plt.savefig('diagnosi_bias_posizione.png', dpi=150)
    print("\nGrafico salvato in diagnosi_bias_posizione.png")
