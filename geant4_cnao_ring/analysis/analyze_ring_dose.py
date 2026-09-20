"""
Analizza i CSV prodotti da cnaoRingImpact: mappa di dose/ionizzazioni
per posizione (xi) sull'anello di 100 "neuroni", per ciascuna
configurazione (particella, energia) del CNAO simulata.

Collega esplicitamente ai risultati del progetto neuroni-progetto:
sovrappone le posizioni usate per i messaggi nella comunicazione
spiking multi-messaggio (x=0.15, 0.40, 0.65, 0.70 -- sezioni 5.10,
5.18, 5.19) alla mappa di dose fisica, per vedere se il fascio
colpisce in modo diseguale le posizioni usate funzionalmente dalla
rete o le attraversa in modo sostanzialmente uniforme.

Uso:
    python3 analyze_ring_dose.py results/
"""

import sys
import csv
import glob
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

POSIZIONI_MESSAGGI = [0.10, 0.15, 0.40, 0.65, 0.70]  # tutte le posizioni usate
                                                        # nei vari esperimenti spiking


def carica_csv(path):
    xi, sum_edep, mean_edep, hit_count, hit_frac = [], [], [], [], []
    n_events = None
    with open(path) as f:
        for riga in f:
            if riga.startswith('# n_events='):
                n_events = int(riga.strip().split('=')[1])
                continue
            if riga.startswith('neuron_index'):
                continue
            parti = riga.strip().split(',')
            if len(parti) != 6:
                continue
            _, x, se, me, hc, hf = parti
            xi.append(float(x))
            sum_edep.append(float(se))
            mean_edep.append(float(me))
            hit_count.append(int(hc))
            hit_frac.append(float(hf))
    return {
        'n_events': n_events,
        'xi': np.array(xi),
        'sum_edep_MeV': np.array(sum_edep),
        'mean_edep_keV': np.array(mean_edep),
        'hit_count': np.array(hit_count),
        'hit_frac': np.array(hit_frac),
    }


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else 'results'
    files = sorted(glob.glob(os.path.join(results_dir, 'dose_*.csv')))
    if not files:
        print(f"Nessun file dose_*.csv trovato in {results_dir}")
        return

    dati = {}
    for f in files:
        nome = os.path.basename(f).replace('dose_', '').replace('.csv', '')
        dati[nome] = carica_csv(f)
        d = dati[nome]
        print(f"{nome}: n_events={d['n_events']}, "
              f"dose totale sull'anello={d['sum_edep_MeV'].sum():.2f} MeV, "
              f"neuroni con almeno un urto diretto del primario={int((d['hit_count']>0).sum())}/100, "
              f"neuroni con dose>0={int((d['sum_edep_MeV']>0).sum())}/100")

    # --- Figura 1: mappa di dose media per posizione, tutte le configurazioni ---
    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

    ax = axes[0]
    for nome, d in dati.items():
        colore = 'C0' if 'proton' in nome else 'C3'
        stile = '-' if '070' in nome or '120' in nome else ('--' if '150' in nome or '260' in nome else ':')
        ax.plot(d['xi'], d['mean_edep_keV'], stile, color=colore, alpha=0.8, label=nome)
    for xm in POSIZIONI_MESSAGGI:
        ax.axvline(xm, color='gray', linestyle=':', linewidth=0.7)
    ax.set_ylabel('Dose media per evento (keV) per neurone')
    ax.set_title("Mappa di dose per posizione sull'anello -- protoni (blu) vs carbonio-12 (rosso)")
    ax.legend(fontsize=7, ncol=2)

    ax = axes[1]
    for nome, d in dati.items():
        colore = 'C0' if 'proton' in nome else 'C3'
        stile = '-' if '070' in nome or '120' in nome else ('--' if '150' in nome or '260' in nome else ':')
        ax.plot(d['xi'], d['hit_frac'], stile, color=colore, alpha=0.8, label=nome)
    for xm in POSIZIONI_MESSAGGI:
        ax.axvline(xm, color='gray', linestyle=':', linewidth=0.7)
    ax.set_xlabel('Posizione sull\'anello (xi)')
    ax.set_ylabel('Frazione di eventi con urto diretto del primario')
    ax.set_title("Frazione di attraversamenti diretti per posizione")
    ax.legend(fontsize=7, ncol=2)

    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, 'mappa_dose_anello.png'), dpi=150)
    print(f"\nGrafico salvato in {os.path.join(results_dir, 'mappa_dose_anello.png')}")

    # --- Tabella riassuntiva: dose totale e picco/media, confronto
    #     protoni vs carbonio a energia comparabile ---
    print("\n" + "="*78)
    print("RIASSUNTO: dose totale sull'anello e rapporto picco/media per neurone")
    print("="*78)
    print(f"{'configurazione':<22}{'dose_tot_MeV':>14}{'dose_media_keV':>16}{'picco/media':>14}")
    for nome, d in dati.items():
        tot = d['sum_edep_MeV'].sum()
        media = d['mean_edep_keV'].mean()
        picco = d['mean_edep_keV'].max()
        rapporto = picco/media if media > 0 else float('nan')
        print(f"{nome:<22}{tot:>14.2f}{media:>16.4f}{rapporto:>14.2f}")

    # --- Confronto specifico sulle posizioni usate per i messaggi ---
    print("\n" + "="*78)
    print("DOSE MEDIA NELLE POSIZIONI USATE DAI MESSAGGI vs RESTO DELL'ANELLO")
    print("="*78)
    for nome, d in dati.items():
        idx_msg = np.array([np.argmin(np.abs(d['xi'] - xm)) for xm in POSIZIONI_MESSAGGI])
        mask_msg = np.zeros(len(d['xi']), dtype=bool)
        mask_msg[idx_msg] = True
        dose_msg = d['mean_edep_keV'][mask_msg].mean()
        dose_resto = d['mean_edep_keV'][~mask_msg].mean()
        print(f"  {nome}: posizioni messaggio={dose_msg:.4f} keV, resto anello={dose_resto:.4f} keV, "
              f"rapporto={dose_msg/dose_resto if dose_resto>0 else float('nan'):.2f}")


if __name__ == '__main__':
    main()
