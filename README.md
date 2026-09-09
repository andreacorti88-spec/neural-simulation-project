# Dalla comunicazione tra due neuroni a una politica decisionale appresa

Un percorso completo di simulazione neurale computazionale — dalla scala di un singolo neurone fino a 40 milioni, dalla propagazione di un segnale fino a un rudimento di apprendimento per rinforzo.

**Autore:** Andrea Corti — Università di Pavia (TLB / Fisica)
**Periodo:** Settembre 2026
**Ambiente:** Python 3.14, Brian2 2.10.1, NumPy, SciPy, Matplotlib — MacBook Air M5, 24 GB RAM

![Rete a 118.000 neuroni](figures/rete_118000_neuroni_t200ms.png)

## Cosa contiene questo repository

Sedici esperimenti progressivi, ciascuno costruito sul precedente, con codice verificato eseguendolo realmente (non solo scritto):

- **Comunicazione di base:** soglia sinaptica, plasticità a breve termine, propagazione, sommazione temporale
- **Apprendimento:** STDP, classificazione di pattern
- **Scala:** rete bilanciata E/I fino a 40 milioni di neuroni (57% della scala di un cervello di topo)
- **Struttura spaziale:** connettività small-world, forma 3D a cervello, visualizzazioni interattive
- **Attrattori:** ring attractor verificato (persistenza, multistabilità, robustezza), catene di attrattori collegati
- **Decisione:** competizione tra alternative, apprendimento per rinforzo, politica contestuale
- **Sintesi finale:** struttura spaziale + sinapsi biofisicamente realistiche + densità sinaptica pari al modello più avanzato al mondo (Kuriyama, Akira et al. 2025)

Il resoconto completo (91 pagine, con ogni esperimento spiegato in dettaglio, codice sorgente incluso, ed equazioni) è in [`report/resoconto_progetto_ESTESO.pdf`](report/resoconto_progetto_ESTESO.pdf) (anche in [`.docx`](report/resoconto_progetto_ESTESO.docx)). Una descrizione più breve è in [`README.pdf`](README.pdf).

## Struttura del repository

| Cartella | Contenuto |
|---|---|
| [`/scripts`](scripts) | Tutti gli script Python, uno per esperimento |
| [`/figures`](figures) | I grafici prodotti da ogni esperimento |
| [`/visualizations`](visualizations) | Le visualizzazioni 3D interattive (HTML, autonome) |
| [`/report`](report) | Il resoconto completo (PDF/DOCX) |

## Come eseguire

Ogni script è autonomo. Per gli esperimenti a grande scala, attesi tempi di esecuzione:

- **40 milioni di neuroni:** ~10 minuti (singolo blocco) fino a diverse ore (esperimento completo multi-fase, throttling termico incluso)
- **118.000 neuroni a densità realistica:** ~20-25 minuti

```bash
pip install brian2 numpy scipy matplotlib --break-system-packages
python3 scripts/due_neuroni.py   # il punto di partenza più semplice
```

## Nota metodologica

Il progetto include deliberatamente i tentativi falliti (un attrattore che non ha funzionato, un'instabilità scoperta a piena scala e poi corretta) oltre ai risultati positivi. Un risultato negativo verificato rigorosamente vale quanto un risultato positivo — è la disciplina che ha guidato l'intero lavoro.

## Licenza

Progetto personale, a scopo di apprendimento e portfolio. Codice liberamente riutilizzabile citando la fonte.
