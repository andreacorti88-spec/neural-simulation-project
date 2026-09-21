# Dalla comunicazione tra due neuroni a una politica decisionale appresa

Un percorso completo di simulazione neurale computazionale — dalla scala di un singolo neurone fino a 40 milioni, dalla propagazione di un segnale fino a un rudimento di apprendimento per rinforzo.

**Autore:** Andrea Corti — Università di Pavia (TLB / Fisica)
**Periodo:** Settembre 2026
**Ambiente:** Python 3.14, Brian2 2.10.1, NumPy, SciPy, Matplotlib — MacBook Air M5, 24 GB RAM

![Rete a 118.000 neuroni](figures/rete_realistica_118K.png)

*([English version](README_EN.md) available — shorter Italian sections, full English project summary in [`report/PROJECT_SUMMARY_EN.md`](report/PROJECT_SUMMARY_EN.md).)*

## Cosa contiene questo repository

Ventitré esperimenti progressivi, ciascuno costruito sul precedente, con codice verificato eseguendolo realmente (non solo scritto):

- **Comunicazione di base:** soglia sinaptica, plasticità a breve termine, propagazione, sommazione temporale
- **Apprendimento:** STDP, classificazione di pattern
- **Scala:** rete bilanciata E/I fino a 40 milioni di neuroni (57% della scala di un cervello di topo)
- **Struttura spaziale:** connettività small-world, forma 3D a cervello, visualizzazioni interattive
- **Attrattori:** ring attractor verificato (persistenza, multistabilità, robustezza), catene di attrattori collegati
- **Decisione:** competizione tra alternative, apprendimento per rinforzo, politica contestuale
- **Sintesi finale:** struttura spaziale + sinapsi biofisicamente realistiche + densità sinaptica pari al modello più avanzato al mondo (Kuriyama, Akira et al. 2025)
- **Comunicazione bidirezionale appresa (dominio spiking):** dal fallimento totale alla soluzione completa (interruttore di trasmissione + normalizzazione omeostatica + soppressione di emergenza), verificata su 8 seed indipendenti, estesa a più messaggi e a una rete a tre popolazioni
- **Fisica delle radiazioni (Geant4), su due scale collegate:** impatto di fasci ionici clinici del CNAO (protoni, carbonio-12) sulla stessa architettura ad anello, dal plateau al vero picco di Bragg — collegato direttamente alla rete spiking tramite un esperimento di knockout funzionale, con un gradiente dose-risposta verificato, un recupero funzionale, un punto di non ritorno, un modello RBE clinico (Kase/NIRS-Chiba) validato, e infine collegato a livello nanodosimetrico (Geant4-DNA) agli elettroni secondari realmente prodotti nei neuroni colpiti — la catena fisica macroscopica-a-nanoscopica verificata end-to-end

Il resoconto completo (quasi 400 pagine, con ogni esperimento spiegato in dettaglio, codice sorgente incluso, ed equazioni) è in [`report/resoconto_progetto_ESTESO.pdf`](report/resoconto_progetto_ESTESO.pdf) (anche in [`.docx`](report/resoconto_progetto_ESTESO.docx) e in [LaTeX](latex/), compilato in [`latex/main.pdf`](latex/main.pdf)). Una descrizione più breve è in [`README.pdf`](README.pdf).

## Struttura del repository

| Cartella | Contenuto |
|---|---|
| [`/scripts`](scripts) | Tutti gli script Python, uno per esperimento |
| [`/figures`](figures) | I grafici prodotti da ogni esperimento |
| [`/visualizations`](visualizations) | Le visualizzazioni 3D interattive (HTML, autonome) |
| [`/report`](report) | Il resoconto completo (PDF/DOCX) |
| [`/geant4_cnao_ring`](geant4_cnao_ring) | Progetto Geant4 (C++): impatto fisico dei fasci CNAO sull'architettura neurale — vedi il README dedicato |
| [`/geant4dna_icsd`](geant4dna_icsd) | Progetto Geant4-DNA (C++): nanodosimetria a livello di elettroni secondari, collegato a `geant4_cnao_ring` (sezioni 5.39-5.41) — vedi il README dedicato |
| [`/geant4dna_ssbdsb`](geant4dna_ssbdsb) | Progetto Geant4-DNA (C++, esempio `moleculardna`): danno reale al DNA (SSB/DSB) da elettroni secondari CNAO (sezione 5.42) — vedi il README dedicato |

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
