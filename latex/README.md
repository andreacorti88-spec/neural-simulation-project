# Versione LaTeX del resoconto

Conversione automatica di `resoconto_progetto_ESTESO.docx` (la fonte
originale, Italiano) in LaTeX, generata con pandoc e poi ripulita da
uno script (`build_latex.py`) che sostituisce ogni blocco di codice
sorgente inserito nel docx come tabella con un vero
`\lstinputlisting` che punta al file sorgente reale — molto piu'
pulito della trascrizione automatica di pandoc, e sempre sincronizzato
con il codice vero.

## Compilare

Richiede un motore LaTeX in grado di gestire Unicode e font moderni
(XeLaTeX o LuaLaTeX). Su Overleaf: caricare l'intera cartella
`latex/` come progetto, impostare il motore su XeLaTeX, compilare
`main.tex`.

In locale, con [tectonic](https://tectonic-typesetting.github.io/)
(motore LaTeX autosufficiente, scarica i pacchetti al volo):

```bash
brew install tectonic   # macOS
tectonic main.tex
```

Oppure con una distribuzione TeX Live/MacTeX completa:

```bash
xelatex main.tex
xelatex main.tex   # una seconda volta per indice/riferimenti
```

Compilato e verificato: 397 pagine, nessun errore, solo warning
cosmetici di *overfull hbox* su alcune righe di codice molto lunghe
(non troncano nulla, solo un'estetica non perfetta).

## Rigenerare da zero

Se il resoconto (`.docx`) viene aggiornato con nuove sezioni:

```bash
# 1. aggiornare le copie dei sorgenti usate dai \lstinputlisting
cp ../*.py scripts/
cp ../geant4_cnao_ring/src/*.cc ../geant4_cnao_ring/include/*.hh \
   ../geant4_cnao_ring/*.cc ../geant4_cnao_ring/analysis/*.py geant4_src/

# 2. rigenerare main.tex dal docx aggiornato
python3 build_latex.py

# 3. ricompilare
tectonic main.tex
```

Richiede `pandoc` (`brew install pandoc`).

## Struttura

| File/cartella | Contenuto |
|---|---|
| `main.tex` | Documento LaTeX completo, generato — non modificare a mano (verrebbe sovrascritto al prossimo `build_latex.py`) |
| `main.pdf` | PDF compilato |
| `build_latex.py` | Script di conversione/pulizia (docx → main.tex) |
| `scripts/` | Copie dei sorgenti Python, referenziate dai `\lstinputlisting` |
| `geant4_src/` | Copie dei sorgenti C++/Python del progetto Geant4 |
| `media/` | Immagini estratte dal docx |

## Nota

Per il resto (README in inglese, sintesi del progetto in inglese)
vedi [`../README_EN.md`](../README_EN.md) e
[`../report/PROJECT_SUMMARY_EN.md`](../report/PROJECT_SUMMARY_EN.md) —
questa cartella e' solo la trasposizione LaTeX del resoconto italiano
completo, non una sua traduzione.
