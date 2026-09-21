# cnaoRingImpact — impatto fisico di fasci ionici CNAO su un'architettura a ring attractor

*([English version](README_EN.md))*

**Parte integrante di `neuroni-progetto`** (non un progetto satellite):
vive in `neuroni-progetto/geant4_cnao_ring/`, stesso repository Git dei
codici Brian2. Due linguaggi (C++/Geant4 per la fisica, Python/Brian2 per
la rete spiking) e due toolchain di build separate per necessita'
tecnica, ma un solo progetto, una sola cronologia, un solo resoconto.

Calcola la dose e il numero di attraversamenti diretti depositati da
fasci ionici clinici del CNAO (protoni e carbonio-12) su una geometria
fisica che rappresenta il ring attractor spiking dello stesso progetto
(sezioni 5.6-5.20 del resoconto): un anello di 100 sfere tessuto-
equivalenti, una per ciascun neurone eccitatorio della rete, posizionate
esattamente alle stesse coordinate angolari (`xi = i/100`) usate nel
modello Brian2. Il collegamento e' concreto, non solo tematico: la
sezione "Collegamento a Brian2" piu' sotto usa direttamente le mappe di
dose calcolate qui per spegnere neuroni specifici in
`dialogo_spiking_cnao_knockout.py`, nella stessa directory principale
del progetto.

Progetto gemello di `geant4dna_icsd` (nanodosimetria a livello di
elettroni, sezione dedicata alla verifica del paper Villagrasa/Baiocco
2026), ma con uno scopo e una fisica completamente diversi — vedi sotto.

## Perché QBBC e non Geant4-DNA

`geant4dna_icsd` usa la fisica Geant4-DNA (opzioni 2/4/6), valida **solo
fino a ~1 MeV e solo per elettroni** in acqua liquida (vedi
`PhysicsList.cc` di quel progetto, `SetMaxEnergy(1*MeV)`). I fasci
clinici del CNAO sono protoni da 60-250 MeV e ioni carbonio-12 da
120-400 MeV/u: decine-centinaia di MeV per nucleone, un dominio
energetico e un tipo di trasporto (comprese le interazioni nucleari di
frammentazione del carbonio) per cui Geant4-DNA non ha un modello
validato. Qui si usa **QBBC**, la physics list di riferimento di Geant4
esplicitamente raccomandata (ed usata nell'esempio ufficiale
"Hadrontherapy") per protoni e ioni a energie cliniche.

Le due fisiche non sono intercambiabili né componibili direttamente in
un solo progetto senza un lavoro multi-scala serio (vedi "Sviluppi
futuri" sotto) — è la stessa distinzione discussa nel documento
`Analisi_critica_nanodosimetria.pdf`, dove la nanodosimetria (Geant4-DNA)
e il trasporto macroscopico del fascio primario sono esplicitamente due
livelli separati della stessa catena fisica.

## Geometria

- **Mondo**: blocco d'acqua (tessuto-equivalente) 4mm x 4mm x 4mm.
- **Anello**: 100 sfere ("Neuron"), raggio 10um (soma corticale
  tipico, ~20um di diametro), disposte su un cerchio di raggio 500um
  (anello di ~1mm di diametro) nel piano Z=0, alla posizione angolare
  `xi = i/100`, `i=0..99` — stessa indicizzazione delle posizioni
  usate per i messaggi nei dialoghi spiking (es. x=0.15, x=0.65).
- **Fascio**: pencil beam lungo +Z, posizione laterale (x,y) campionata
  uniformemente su un disco di raggio 550um che copre tutto l'anello.
  Un vero spot CNAO clinico ha sigma di alcuni mm, molto più largo
  dell'anello: su questa scala microscopica il fascio è quindi
  effettivamente uniforme, e qualunque differenza di dose tra posizioni
  osservata nei risultati è vera stocasticità microdosimetrica
  (fluttuazioni statistiche reali nel numero e nel tipo di interazioni),
  non un artefatto della forma del fascio.

**Limite dichiarato**: la scala fisica (anello 1mm, soma 20um) è una
scelta di letteratura per un microcircuito corticale tipico, non una
misura della rete simulata in Brian2 (che è un modello topologico, non
metrico — non ha mai avuto una scala fisica propria). Il fascio entra a
~2mm di profondità nel blocco, un punto qualunque della fase di plateau
della curva di Bragg per queste energie: questo progetto NON simula il
picco di Bragg né la dipendenza dalla profondità reale di trattamento,
solo l'impatto locale a una profondità fissa.

## Cosa viene registrato

Per ciascun neurone, su tutta la durata di un run (accumulato su tutti
gli eventi primari simulati):
- **energia depositata totale** (da tutte le tracce, primaria e
  secondarie, incluse le interazioni nucleari del carbonio)
- **numero di eventi in cui la traccia PRIMARIA** (protone o ione,
  `ParentID==0`) attraversa fisicamente quel neurone

Output: un CSV per run (`neuron_dose_run0.csv`), rinominato per
configurazione da `analysis/analyze_ring_dose.py`.

## Risultati (50.000 primari per configurazione)

| Configurazione | Dose totale anello (MeV) | Dose media/evento (keV) |
|---|---|---|
| protone 70 MeV | 22.5 | 0.0045 |
| protone 150 MeV | 12.4 | 0.0025 |
| protone 225 MeV | 9.5 | 0.0019 |
| carbonio-12 120 MeV/u | 509.4 | 0.102 |
| carbonio-12 260 MeV/u | 288.5 | 0.058 |
| carbonio-12 400 MeV/u | 239.9 | 0.048 |

Due osservazioni fisicamente coerenti, entrambe attese:

1. **Il carbonio-12 deposita 10-25 volte più dose per neurone del
   protone** a energie per nucleone comparabili — coerente con il LET
   molto più alto di uno ione pesante (scala circa con Z², qui Z=6
   contro Z=1).
2. **La dose diminuisce con l'energia** per entrambe le particelle —
   a energia più alta la particella è più lontana dal proprio picco di
   Bragg a questa profondità fissa, quindi ha un potere frenante (LET)
   minore.

**Confronto con le posizioni funzionali della rete**: nessuna delle sei
configurazioni mostra una differenza sistematica tra la dose media
sulle posizioni usate per i messaggi nei dialoghi spiking (x=0.10,
0.15, 0.40, 0.65, 0.70) e il resto dell'anello (rapporto tra 0.97 e
1.20, compatibile con puro rumore statistico) — conferma diretta che
il fascio, su questa scala, non privilegia alcuna posizione: qualunque
neurone della rete corre lo stesso rischio dosimetrico fisico,
indipendentemente dal suo ruolo funzionale nella rete.

⚠️ **Questa tabella usa la geometria originale (disco di campionamento
550um), poi scoperta problematica in profondita' (sezione "Limiti"
sotto e sezione 5.25 del resoconto)**. Le conclusioni qualitative
restano valide (tutte le run sono alla stessa profondita' fissa, dove
l'artefatto colpisce le posizioni in modo uniforme), ma i numeri
assoluti di dose sono superati dalla tabella seguente.

### Risultati corretti: plateau vs vero picco di Bragg (`results_picco/`)

Con la geometria corretta (disco 5mm, sigma di uno spot clinico reale)
e l'anello posizionato ESATTAMENTE al picco di Bragg per ciascuna
energia (trovato con `braggProfile/`, validato entro l'1% contro la
letteratura):

| Configurazione | Profondità picco | Dose plateau (MeV) | Dose al picco (MeV) | Rapporto |
|---|---|---|---|---|
| protone 70 MeV | 39.5mm | 0.118 | 0.953 | 8.11x |
| protone 150 MeV | 156.5mm | 0.099 | 0.499 | 5.06x |
| protone 225 MeV | 315.5mm | 0.104 | 0.996 | 9.56x |
| carbonio-12 120 MeV/u | 34.5mm | 7.068 | 20.996 | 2.97x |
| carbonio-12 260 MeV/u | 135.5mm | 3.504 | 11.258 | 3.21x |
| carbonio-12 400 MeV/u | 275.5mm | 3.254 | 4.412 | 1.36x |

Tutte e sei le configurazioni mostrano il picco fisicamente più alto
del plateau, come atteso — il segno distintivo del vantaggio clinico
dell'adroterapia (concentrare la dose sul bersaglio, risparmiare il
tessuto sano attraversato prima). Nota di rigore: il primo tentativo
per il protone 225 MeV (il caso più profondo e col picco più stretto)
dava 50.000 eventi un rapporto anomalo di 0.24x (picco più BASSO del
plateau) — un artefatto di pura statistica Monte Carlo insufficiente
su un picco molto stretto a grande profondità, non un errore
sistematico: risolto con 500.000 eventi (rapporto 9.56x, coerente con
gli altri cinque casi). Lezione tenuta a mente per qualunque estensione
futura a profondità ancora maggiori.

## Limiti di questo primo livello

- **Statistica per neurone bassa**: con probabilità geometrica di
  attraversamento diretto ~0.03% per neurone per evento, 50.000 primari
  danno ~15-25 attraversamenti diretti per neurone — sufficiente per una
  dose media stabile ma con rumore visibile nella mappa punto-per-punto
  (Fig. `mappa_dose_anello.png`). Per una mappa punto-per-punto meno
  rumorosa servirebbero 10-100x più eventi.
- **Nessun modello di danno biologico**: questo progetto calcola solo
  energia depositata (dose fisica), non un endpoint biologico (rotture
  del DNA, morte cellulare) — lo stesso limite dichiarato esplicitamente
  nel paper Villagrasa/Baiocco per il passaggio dose->danno->sopravvivenza
  (vedi `Analisi_critica_nanodosimetria.pdf`, sezione 5).
- **Profondità del picco di Bragg — RISOLTO**: i risultati sopra sono
  tutti a profondità fissa (~2mm, plateau). `braggProfile/` (sotto-
  progetto separato) trova il vero picco per una data energia
  scansionando la dose in profondità in un fantoccio a fette (validato
  entro l'1% contro valori di riferimento noti: 150 MeV protoni →
  picco a 156.5mm, letteratura ~157.7mm). `cnaoRingImpact` ora accetta
  la profondità dell'anello come parametro. Durante il confronto
  plateau-vs-picco è emerso un artefatto reale: il disco di
  campionamento laterale originale (550um) è troppo stretto a
  profondità maggiori, dove la diffusione multipla del fascio diluisce
  la fluenza primaria in quella zona — dava dose PIÙ BASSA al picco
  che in plateau, il contrario della fisica attesa. Diagnosticato
  confrontando con i dati indipendenti di `braggProfile` (stessa
  energia, fetta larga 5x5cm, insensibile all'artefatto) e corretto
  allargando il disco a 5mm (sigma di uno spot clinico CNAO reale,
  sempre molto più largo della diffusione multipla anche a diversi cm
  di profondità). Dopo la correzione: rapporto picco/plateau 2.81x per
  protoni 70 MeV, nella direzione fisica corretta.
- **Nessun collegamento diretto (accoppiato) con la simulazione Brian2**:
  questo e' un output fisico puro (CSV di dose per posizione), analizzato
  separatamente rispetto alla funzione della rete -- non una
  co-simulazione fisica-neurale in tempo reale (che richiederebbe
  tradurre dose/ionizzazioni in un modello di danno/inattivazione
  neuronale con soglie che non hanno, ad oggi, un valore consolidato in
  letteratura per un singolo neurone).

## Collegamento a Brian2: dal dato fisico al knockout funzionale

`neuroni-progetto/dialogo_spiking_cnao_knockout.py` prende i 10 neuroni
a dose piu' alta di una mappa di questo progetto e li rende "spenti"
(sinapsi in uscita rimosse) nella rete spiking multi-messaggio (sezione
5.10 del resoconto). Risultato, verificato su 8 seed e con un controllo
negativo (knockout casuale non basta): il messaggio che usa i neuroni
spenti crolla, quello che non li usa resta intatto -- un effetto
causale, non un artefatto, dimostrato ribaltandolo deliberatamente: con
il set di neuroni corretto dopo la scoperta dell'artefatto sopra (che
per caso include piu' neuroni della zona dell'altro messaggio),
l'effetto si ribalta esattamente come previsto. Dettagli completi:
resoconto, sezioni 5.22-5.26.

## Sviluppi futuri

1. ~~Statistica maggiore~~ **TESTATO** (sezione 5.37): 10x eventi
   (5.000.000 invece di 500.000) migliora la sovrapposizione top-10 solo
   da 1.0/10 (puro rumore) a 2/10 -- non basta. Estrapolando lo stesso
   scaling servirebbe un ulteriore fattore 100-1000, computazionalmente
   proibitivo nei tempi di questo progetto. Confermato che mediare la
   dose su piu' seed a statistica modesta (5.33) e' la strada giusta,
   non aumentare brutalmente gli eventi per run.
2. ~~Livello nanodosimetrico locale~~ **FATTO** (sezioni 5.39-5.40):
   collegato a `geant4dna_icsd` (ora anch'esso dentro `neuroni-progetto/`).
   Aggiunto il tracking dello spettro degli elettroni secondari nati
   dentro un neurone (scoperto e corretto un bug: il taglio di
   produzione di default era troppo grosso, zero secondari registrati
   finche' non e' stata aggiunta una G4Region dedicata con taglio
   100nm). Spettro verificato: 100% sotto 1 MeV (mediana 1.86 keV,
   coerente col dominio di Geant4-DNA). Usate quelle energie reali come
   primari in `nanoICSD` (mai compilato/eseguito prima d'ora): gli ICSD
   risultanti si inseriscono coerentemente tra i valori gia' pubblicati
   nel paper Villagrasa/Baiocco -- il collegamento macro-a-nano e'
   verificato end-to-end. Vedi il README di `geant4dna_icsd` per i
   dettagli completi.
3. ~~Soglia di danno neuronale assoluta~~ **TESTATO** (sezione 5.38): il
   knockout resta un criterio relativo (una soglia assoluta di danno
   biologico non e' comunque consolidata in letteratura per un singolo
   neurone, lo stesso limite dichiarato nel paper Villagrasa/Baiocco),
   ma la sensibilita' alla DIMENSIONE del knockout e' stata testata:
   TOP-5 (overlap 0-vs-0 con le zone di codifica) -> effetto nullo come
   previsto; TOP-20 (overlap 3-vs-1, segno OPPOSTO al TOP-10 gia'
   testato) -> l'effetto si inverte nettamente (1.74x, 7/8 seed) esattamente
   come previsto. Il meccanismo non dipende dalla scelta specifica di
   "esattamente 10 neuroni".
4. ~~Knockout al vero picco di Bragg~~ **FATTO** (sezione 5.27 del
   resoconto): con la mappa di dose al picco (non piu' in plateau), il
   set dei 10 neuroni piu' colpiti cambia ancora, con una sovrapposizione
   quasi bilanciata tra i due messaggi (~1 neurone ciascuno) -- e
   l'effetto del knockout, come previsto dal meccanismo gia' verificato,
   diventa nullo (rapporto 1.06x, 4/8 seed). Messo in fila con i due
   test precedenti (2-vs-0 -> 4.7x; 3-vs-1 -> 1.81x; ~1-vs-1 -> 1.06x),
   emerge un gradiente dose-risposta pulito: la forza dell'effetto scala
   con quanto sono sbilanciati i neuroni spenti tra i due messaggi, non
   un fenomeno di tutto-o-niente.

## Aggiornamenti (sezioni 5.28-5.42)

- **5.28 -- Gradiente confermato con i protoni**: mappa di dose al picco
  di Bragg per protone 70 MeV, sbilanciamento 4-vs-0 (il piu' netto
  finora) -> rapporto 7.07x, 8/8 seed, l'effetto piu' forte e piu'
  consistente della serie.
- **5.29 -- Recupero funzionale dopo il danno**: training esteso 5x (300
  ripetizioni) dopo il knockout piu' severo (5.28): il messaggio
  danneggiato recupera progressivamente (+294% tra prima e seconda meta'
  del training, verificato su 4 seed), senza mai eguagliare il messaggio
  sano -- plasticita' STDP residua che compensa, non ripara.
- **5.30 -- Punto di non ritorno**: spegnendo TUTTI gli 11 neuroni della
  zona di codifica di un messaggio (non piu' un sotto-insieme a dose
  fisica), il recupero della 5.29 sparisce del tutto, confermato su 5
  run indipendenti (1 singolo + 4 seed MC) -- la plasticita' compensativa
  richiede che sopravviva almeno una via d'uscita funzionante.
- **5.31 -- Modello RBE-LET (Kanai/NIRS-Chiba)**: aggiunto il tracking
  del LET dose-mediato per neurone (stesso formalismo usato a NIRS/HIMAC
  per il calcolo clinico della dose biologica). Il divario RBE tra
  particelle e' netto (carbonio ~3.3, protone ~1.15 al rispettivo picco,
  coerente col riferimento clinico RBE=3.0 a 80 keV/um), ma entro la
  stessa profondita' la classifica dei neuroni cambia pochissimo (9/10
  di sovrapposizione) -- il gradiente 5.22-5.28 regge al raffinamento.
- **5.32 -- Limite scoperto: la mappa di dose per-neurone non e'
  riproducibile**: aggiunto un seed esplicito a `cnaoRingImpact` (prima
  basato sull'orologio, diverso a ogni run). Ripetendo la stessa
  configurazione fisica su 4 seed indipendenti, la sovrapposizione media
  tra le classifiche top-10 e' solo 1.0/10 -- dominata dal rumore di
  Poisson (~1-2 colpi diretti attesi per neurone su 500k eventi). La
  dose TOTALE sull'anello resta stabile (+-15%); e' solo la sua
  distribuzione tra i 100 neuroni a essere rumorosa.
- **5.33 -- Il gradiente su basi statistiche solide**: ricostruito il
  set di knockout dalla dose MEDIATA sui 4 seed della 5.32 (overlap
  1-vs-0 con le zone di codifica) -> rapporto 1.29x, 6/8 seed,
  statisticamente indistinguibile dal controllo negativo 0-vs-0 (1.39x,
  6/8, sezione 5.24). Emerge una soglia: serve uno sbilanciamento di
  ALMENO 2 neuroni perche' l'effetto superi il rumore di base della
  rete -- il gradiente si raffina, non si invalida.
- **5.34 -- Il vero modello MKM di Kase/NIRS**: sostituita l'approssimazione
  lineare RBE(LET) della 5.31 con il modello MKM (Microdosimetric
  Kinetic Model) saturazione-corretto usato clinicamente a NIRS/HIMAC
  (Kase et al. 2011, J Radiat Res 52:59-68 -- lo stesso paper che
  descrive il TEPC usato di routine per la QA della dose biologica a
  HIMAC), con i parametri reali delle cellule HSG. Validato PRIMA di
  pubblicare: applicato ai valori misurati e pubblicati da Kase et al.,
  riproduce i loro RBE10 entro il 10-15% (es. al punto di riferimento
  clinico NIRS, LET_d=80 keV/um -> RBE10 calcolato 2.67 contro il
  riferimento clinico 3.0). Applicato ai nostri dati: RBE10 fino a 3.44
  per il carbonio al picco (fisicamente sensato), ~1.0 per il protone
  (come atteso) -- e la conclusione della 5.31 (il gradiente 5.22-5.28
  regge al raffinamento RBE) si conferma anche col modello clinico
  vero, non era un artefatto dell'approssimazione lineare. Script
  riutilizzabile: `analysis/mkm_rbe.py`.
- **5.35 -- Il modello MKM su tutte e sei le configurazioni CNAO**:
  completate le quattro configurazioni mancanti con tracking del LET.
  Sovrapposizione fisica/biologica 9-10/10 in TUTTE e sei -- conferma
  generalizzata. Due osservazioni non banali: overkill del carbonio
  (RBE10 max si stabilizza ~3.44 anche se il LET massimo cresce) e un
  RBE anomalo per il protone a 150 MeV (fino a 1.46).
- **5.36 -- Verificato, non ipotizzato**: diagnostica per-step aggiunta
  per identificare la causa dell'anomalia della 5.35. Non e' un
  frammento secondario come ipotizzato: e' il protone primario stesso,
  quasi fermo (range straggling) -- lo stesso meccanismo dietro il
  dibattito reale sull'RBE elevato al bordo distale in protonterapia.
- **5.37 -- Quantificato il costo di una statistica maggiore**: 10x
  eventi (5.000.000 invece di 500.000) migliora la sovrapposizione
  top-10 solo da 1.0/10 (puro rumore) a 2/10 -- non basta. Confermato
  che mediare su seed (5.33) e' molto piu' efficiente che aumentare gli
  eventi per run.
- **5.38 -- La dimensione del knockout non e' magica**: testato TOP-5
  (overlap 0-vs-0 con le zone di codifica -> effetto nullo come
  previsto) e TOP-20 (overlap 3-vs-1, segno OPPOSTO al TOP-10 gia'
  testato -> l'effetto si inverte nettamente, 1.74x, 7/8 seed,
  esattamente come previsto). Il meccanismo generalizza oltre la
  scelta specifica di "esattamente 10 neuroni" -- quarta previsione
  deliberatamente falsificabile confermata in questo progetto.
- **5.39-5.41 -- Collegamento nanodosimetrico completato e verificato**:
  tracking dello spettro di elettroni secondari nati dentro un neurone
  (scoperto e corretto un bug: taglio di produzione troppo grosso, zero
  secondari finche' non e' stata aggiunta una G4Region dedicata).
  Spettro verificato: 100% sotto 1 MeV, dentro il dominio di Geant4-DNA.
  Quelle energie reali usate come primari in `nanoICSD` (mai eseguito
  prima): gli ICSD risultanti si inseriscono coerentemente tra i valori
  gia' pubblicati nel paper Villagrasa/Baiocco. Verifica MC su 4 seed
  (5.41): M1=1.090 confermato stabile (CV 0.34%) -- a differenza della
  mappa di dose per-neurone (5.32), qui non c'e' un problema di rumore
  statistico. La catena fisica macro-a-nano e' ora verificata
  end-to-end. Dettagli completi nel README di
  `geant4dna_icsd` (ora anch'esso dentro `neuroni-progetto/`).
- **5.42 -- Oltre gli ICSD: danno reale al DNA**: chiuso l'ultimo
  punto aperto (capitolo 9 del resoconto). Primo tentativo con
  `dnadamage1` bloccato da una dipendenza esterna mancante (verificato,
  non ipotizzato); usato invece `moleculardna` (Geant4-DNA
  collaboration, geometria scaricata automaticamente, nessuno strumento
  esterno necessario). Eseguito con le energie reali della 5.39
  (mediana 1860 eV, p90 7332 eV): frazione di rotture doppie (DSB) piu'
  alta a bassa energia (12.6% contro 7.4%), coerente col LET piu' alto
  a energia piu' bassa. Verificato su 3 seed (aggiunta un'opzione -s,
  il codice originale non ne aveva una): SB stabile (CV 1.8%), DSB piu'
  rumoroso (CV 12.4%, atteso per conteggi piccoli). Dettagli completi
  nel README di `geant4dna_ssbdsb` (nuovo progetto dentro
  `neuroni-progetto/`).

## Build ed esecuzione

Stessi prerequisiti di `nanoICSD` (Geant4 11.x), ma qui bastano i
dataset standard (non serve G4EMLOW/Geant4-DNA):

```bash
source /path/a/geant4-install/bin/geant4.sh
mkdir build && cd build
cmake -DGeant4_DIR=/path/a/geant4-install/lib/cmake/Geant4 ..
make -j$(nproc)

# singola configurazione: particella (0=protone,1=carbonio-12), energia MeV totale,
# profondita' anello mm, [seed esplicito -- default: orologio di sistema, diverso a ogni run]
./cnaoRingImpact macros/run.mac 0 150            # protone 150 MeV, plateau (2mm)
./cnaoRingImpact macros/run.mac 1 3120           # carbonio-12 260 MeV/u (260*12=3120 MeV totali)
./cnaoRingImpact macros/run.mac 0 70 39.5 12345  # protone 70 MeV al picco (39.5mm), seed fisso
                                                  # -- per verifica MC multi-seed, vedi sezione 5.32

python3 ../analysis/analyze_ring_dose.py results/
```
