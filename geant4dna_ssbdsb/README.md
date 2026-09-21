# geant4dna_ssbdsb — danno reale al DNA (SSB/DSB) da elettroni secondari CNAO

**Parte integrante di `neuroni-progetto`**, come `geant4_cnao_ring` e
`geant4dna_icsd`: vive in `neuroni-progetto/geant4dna_ssbdsb/`, stesso
repository Git. Chiude, alla sezione 5.42 del resoconto, il punto
lasciato aperto nel capitolo 9 ("sviluppi futuri"): estendere la
catena fisica gia' verificata (dose macroscopica CNAO -> spettro di
elettroni secondari reali -> ICSD nanodosimetrico, sezioni 5.39-5.41)
fino al danno diretto al DNA -- rotture singole e doppie di filamento
(SSB/DSB).

## Cosa contiene

Copia locale, con una piccola modifica, dell'esempio ufficiale
Geant4-DNA `moleculardna` (Geant4-DNA collaboration,
http://moleculardna.org): simula fisica, fisico-chimica e chimica in
geometrie di DNA molecolare per predire il danno precoce. Usa la
geometria di riferimento "cylinders" (200.000 segmenti di DNA di 216
coppie di basi in un volume di collocazione di 100x30x100nm, pensata
per studi parametrici), scaricata automaticamente da CERN al primo
`cmake` (nessuno strumento esterno necessario).

**Modifica rispetto all'originale**: aggiunta un'opzione `-s <seed>`
a `molecular.cc` (il codice originale non chiamava mai `setTheSeed`,
quindi ogni run senza `-s` riproduce la stessa identica sequenza --
utile per la riproducibilita', ma non permette una verifica Monte
Carlo multi-seed senza questa aggiunta).

## Perche' non dnadamage1

Il primo tentativo, con l'esempio ufficiale piu' semplice `dnadamage1`,
e' stato bloccato da una dipendenza reale: richiede un file di
geometria del DNA (`VoxelStraight.fab2g4dna`) generato da uno
strumento esterno (DnaFabric, https://bitbucket.org/sylMeylan/opendnafabric)
non disponibile e non installabile rapidamente in questo ambiente.
Verificato con una ricerca del file (assente ovunque sul sistema)
prima di abbandonare la strada, non solo ipotizzato.

## Risultati (sezione 5.42)

Eseguito con le energie reali dello spettro di elettroni secondari
nati dentro un neurone colpito da un fascio CNAO (carbonio-12 al
picco di Bragg, sezione 5.39 di `geant4_cnao_ring`): mediana 1860 eV
e 90-esimo percentile 7332 eV.

| Energia reale | N primari | SB totali | DSB totali | SB/primario | DSB/primario | Frazione DSB |
|---|---|---|---|---|---|---|
| 1860 eV (mediana) | 10000 | 1055 | 152 | 0.1055 | 0.0152 | 12.6% |
| 7332 eV (p90) | 5000 | 2688 | 215 | 0.5376 | 0.0430 | 7.4% |

L'elettrone a energia piu' alta produce piu' danno assoluto (traccia
piu' lunga), ma una frazione relativamente MINORE di rotture doppie
-- coerente con un LET piu' basso a energia piu' alta in questo
intervallo (il danno si concentra meno).

**Verifica Monte Carlo** (3 seed indipendenti alla mediana, 5000
eventi ciascuno, `results_mc_check/`): SB stabile (CV 1.8%: 522, 514,
504), DSB piu' rumoroso (CV 12.4%: 60, 77, 68 -- atteso per conteggi
molto piu' piccoli). Media dei 3 seed coerente entro il rumore col run
singolo originale a 10.000 eventi.

## Limiti dichiarati

- La geometria "cylinders" e' un volume di studio parametrico, non un
  intero nucleo cellulare -- questi sono rese (yield) relative per
  primario in quel volume definito, non una previsione assoluta di
  danno a livello di cellula intera.
- Solo due punti energetici testati (mediana e 90-esimo percentile
  dello spettro reale), non l'intero spettro continuo pesato per la
  sua vera distribuzione.
- Il modello di danno usa parametri di default dell'esempio (soglia
  energia diretta 17.5 eV, probabilita' di reazione indiretta 65% sul
  filamento) -- non calibrati specificamente su dati sperimentali per
  questo scenario.

## Build ed esecuzione

```bash
mkdir build && cd build
cmake -DGeant4_DIR=/percorso/a/geant4-install/lib/cmake/Geant4 ..   # scarica la geometria al primo run
make -j$(nproc)

# energia reale mediana (1860 eV), 5000 eventi, seed esplicito
./molecular -m ../cnao_median_1860eV.mac -t 4 -p 2 -s 301

# analisi con ROOT (TTree tuples/classification: colonne SSB, SSBp, 2SSB, DSB, DSBp, DSBpp)
root -l molecular-dna.root
```
