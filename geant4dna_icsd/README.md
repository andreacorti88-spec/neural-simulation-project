# nanoICSD — Geant4-DNA reproduction of the paper's ICSD test case

**Parte integrante di `neuroni-progetto`**, come `geant4_cnao_ring`:
vive in `neuroni-progetto/geant4dna_icsd/`, stesso repository Git.
Dalla sezione 5.39 e' collegato direttamente a `geant4_cnao_ring`: gli
elettroni secondari realmente prodotti da un fascio clinico CNAO
dentro un neurone (tracciati li', con il taglio di produzione corretto
per renderli espliciti) vengono usati QUI come energie primarie reali,
al posto dei soli valori monoenergetici del paper -- vedi
"Collegamento a geant4_cnao_ring" piu' sotto.

Progetto Geant4-DNA per riprodurre la geometria e la fisica usate in:

> Villagrasa C, Baiocco G, Chaoui Z-E-A, et al. "Evaluation of the
> uncertainty in calculating nanodosimetric quantities due to the use
> of different interaction cross sections in Monte Carlo track
> structure codes." PLOS One 2026;21(1):e0340500.
> https://doi.org/10.1371/journal.pone.0340500

Fa parte del "Livello 2" descritto nel documento `Analisi_critica_nanodosimetria.docx`:
verificare in modo indipendente, con un codice pubblicamente disponibile
(Geant4-DNA), l'effetto delle sezioni d'urto sul danno nanodosimetrico
che nel paper originale è testato con un solo codice non pubblico
(PARTRAC).

## Cosa fa, davvero, questo codice

Simula un elettrone monoenergetico che parte dal centro di una sfera
d'acqua nanometrica (8 nm per energie fino a 1 keV, 100 nm per 5 e 10
keV, come nel paper) e conta, per ogni evento, quante interazioni di
ionizzazione avvengono dentro quella sfera — dalla particella primaria
e da tutte le sue secondarie. La sequenza di questi conteggi, evento
per evento, **è** la Ionization Cluster Size Distribution (ICSD) usata
nel paper. Da lì si calcolano M1 (media), F2 e F3 (probabilità di 2+ o
3+ ionizzazioni) esattamente come nelle Tabelle 3–7 dell'articolo.

`analysis/analyze_icsd.py` confronta automaticamente i risultati con i
valori M1/F2 pubblicati per Geant4-DNA opzione 2 (sezioni d'urto
originali, Tabelle 3 e 4) — è il controllo di sanità più immediato:
se il codice è impostato correttamente, i numeri dovrebbero avvicinarsi
a quelli del paper entro l'incertezza statistica.

## Cosa NON fa (ancora)

**Non usa il dataset di sezioni d'urto "comuni"** definito nel paper —
quello resta il vero passo successivo, ed è più impegnativo di quanto
sembri:

1. I valori numerici delle sezioni d'urto comuni sono nel materiale
   supplementare del paper (Figure S1–S3 e i dati sorgente collegati),
   non riprodotti qui — vanno estratti da lì.
2. Iniettarli in Geant4-DNA non è questione di un parametro da
   cambiare: richiede scrivere una sottoclasse di `G4VEmModel` (o
   `G4VDNAModel`, a seconda della versione di Geant4-DNA) che carichi
   quella tabella al posto del modello nativo per ionizzazione,
   eccitazione e scattering elastico, e registrarla nel `PhysicsList`
   al posto di `G4EmDNAPhysics_option2` di default.
3. Il paper stesso segnala che l'operazione, fatta dai gruppi
   originali, ha richiesto accesso al codice sorgente e "una
   conoscenza approfondita" di come ciascun codice gestisce le sezioni
   d'urto internamente — non è un'estensione da weekend.

**Non riproduce il confronto sul danno al DNA (SSB/DSB)** fatto con
PARTRAC: quello richiede un modello geometrico di nucleo cellulare con
cromatina multi-scala, che in Geant4-DNA esiste (vedi gli esempi
ufficiali `dnadamage1`/`molecularDNA` nel repository Geant4-DNA) ma va
integrato a parte, dopo aver validato questo primo livello (la ICSD).
È il passo naturale successivo, non incluso qui.

## Prerequisiti

- Geant4 **11.x** compilato con i dataset Geant4-DNA (`G4EMLOW`)
  installati e le variabili d'ambiente dati (`G4LEDATA`, `G4EMLOW`,
  ecc.) impostate — le stesse richieste da qualunque esempio ufficiale
  Geant4-DNA.
- CMake ≥ 3.16, un compilatore C++17.
- **Compilato ed eseguito con successo** (sezione 5.40, sulla stessa
  installazione Geant4 11.4.0 usata per `geant4_cnao_ring`) — la nota
  storica sotto resta per contesto, ma non e' piu' vera:
  ~~Non è stato possibile compilare né eseguire questo progetto
  nell'ambiente in cui è stato scritto (nessuna installazione Geant4
  disponibile, una build da sorgente richiede tipicamente 1–3 ore anche
  su hardware potente).~~

## Build

```bash
mkdir build && cd build
cmake -DGeant4_DIR=/percorso/a/geant4/lib/cmake/Geant4 ..
make -j$(nproc)
```

## Esecuzione — un singolo punto

```bash
cd build
./nanoICSD macros/run.mac 2 8      # opzione 2, target 8 nm, energia da macro
```

Modifica `/gun/energy` in `macros/run.mac` per cambiare l'energia
dell'elettrone.

## Esecuzione — l'intera griglia di energie del paper

```bash
cd build
chmod +x macros/scan_energies.sh
./macros/scan_energies.sh 2 100000   # opzione 2, 1e5 eventi per punto
python3 ../analysis/analyze_icsd.py results_opt2/
```

Con 1e5 eventi per punto (la statistica usata nel paper) il tempo di
calcolo totale, su una singola CPU, è dell'ordine di alcune ore per
l'intera griglia — parallelizzabile per energia se hai accesso a più
core o a un cluster.

## Prima di fidarti dei numeri

- Verifica il nome del processo di ionizzazione
  (`"e-_G4DNAIonisation"` in `SteppingAction.cc`) contro la tua
  versione di Geant4 — è cambiato tra alcune release major.
- Comincia con pochi eventi (1e3–1e4) per controllare che il codice
  giri end-to-end prima di lanciare run da 1e5+ eventi.
- Se i tuoi M1/F2 divergono sistematicamente (non solo per rumore
  statistico) da quelli pubblicati per l'opzione 2, il problema è quasi
  certamente nella geometria o nel process name, non nella fisica.

## Collegamento a geant4_cnao_ring (sezioni 5.39-5.41)

Il "livello nanodosimetrico locale" elencato come sviluppo futuro nel
README di `geant4_cnao_ring` e' stato completato:

1. **Sezione 5.39**: aggiunto a `geant4_cnao_ring` un tracking dedicato
   dell'energia cinetica di ogni elettrone secondario nato dentro un
   volume "Neuron". Primo tentativo: zero elettroni registrati (il
   taglio di produzione di default di QBBC, ~1mm, e' troppo grande
   rispetto alla scala del neurone, 10um -- quasi nessun raggio delta
   ha un range sufficiente per essere generato come track esplicito).
   Corretto con una `G4Region` dedicata (taglio 100nm). Risultato:
   18.021 elettroni secondari su 500.000 eventi (carbonio-12 al picco
   di Bragg), 100% sotto 1 MeV (mediana 1.86 keV, media 3.57 keV,
   massimo 96.2 keV) -- dentro il dominio di validita' di Geant4-DNA,
   confermato empiricamente e non solo assunto.
2. **Sezione 5.40** (qui): usate le energie reali di quello spettro
   (mediana 1.86 keV, 90-esimo percentile 7.33 keV) come primari in
   `nanoICSD`, al posto dei valori monoenergetici del paper. Risultati
   (`results/icsd_1860eV_median_real.csv`,
   `results/icsd_7332eV_p90_real.csv`):

   | Energia reale | Diametro | M1 | F2 | F3 |
   |---|---|---|---|---|
   | 1860 eV (mediana) | 8nm | 1.090 | 0.257 | 0.142 |
   | 7332 eV (p90) | 100nm | 6.320 | 0.801 | 0.670 |

   Entrambi i punti si inseriscono coerentemente tra i valori gia'
   pubblicati nel paper alla stessa dimensione di bersaglio (8nm:
   600eV->2.90, 1000eV->1.83, il nostro 1860eV->1.09, stessa tendenza
   decrescente; 100nm: 5000eV->7.10, il nostro 7332eV->6.32, 10000eV
   ->4.02, esattamente tra i due) -- nessuna sorpresa, il collegamento
   macro-a-nano e' verificato end-to-end.
3. **Sezione 5.41**: nanoICSD usava lo stesso pattern di seed basato
   sull'orologio di cnaoRingImpact prima della correzione 5.32 --
   verificato se M1=1.090 fosse un dato solido o un run rumoroso.
   Aggiunto un seed esplicito (`./nanoICSD [macro] [opzione] [diametro]
   [seed]`) e ripetuta la stessa configurazione su 4 seed indipendenti:
   M1 = 1.0949, 1.0892, 1.0954, 1.0883 -- coefficiente di variazione
   solo 0.34%. A DIFFERENZA della mappa di dose per-neurone (5.32), qui
   non c'e' un problema di rumore statistico: ogni evento contribuisce
   direttamente alla statistica aggregata, senza la suddivisione in
   100 bin sparsi che rendeva rumorosa la mappa di dose. Il valore
   della 5.40 e' confermato solido, non un artefatto di un seed
   fortunato.

## Struttura del progetto

```
nanoICSD.cc                    main: legge opzione fisica e diametro target da riga di comando
include/ src/                  DetectorConstruction, PhysicsList, Primary/Run/Event/SteppingAction
macros/run.mac                 macro di esempio per un singolo punto
macros/scan_energies.sh        script per scansionare l'intera griglia di energie
analysis/analyze_icsd.py       calcola M1/F2/F3 e confronta con i valori pubblicati (opzione 2)
```
