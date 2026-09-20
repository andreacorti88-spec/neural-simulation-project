#!/usr/bin/env python3
"""
Rigenera main.tex a partire da resoconto_progetto_ESTESO.docx (nella
directory del progetto, un livello sopra).

Pipeline:
1. pandoc converte il docx in LaTeX standalone (preambolo completo,
   indice, immagini estratte in media/).
2. Questo script sostituisce ogni blocco di codice sorgente (inserito
   nel docx come tabella con una riga per linea, in font Consolas --
   pandoc lo converte in prosa illeggibile, con caratteri LaTeX speciali
   mal escapati) con un vero \\lstinputlisting che punta al file
   sorgente reale in scripts/ o geant4_src/ -- molto piu' pulito e
   sempre sincronizzato con il codice vero, non con una sua trascrizione.
3. Corregge alcuni simboli Unicode (≥, →, µ, ecc.) che non hanno un
   glifo nel font di default sotto XeTeX, sostituendoli con l'equivalente
   in modalita' matematica LaTeX.
4. Aggiunge la configurazione del pacchetto listings (syntax highlighting
   di base per Python/C++) subito prima di \\begin{document}.

Uso:
    python3 build_latex.py            # rigenera main.tex
    tectonic main.tex                 # compila in main.pdf (o pdflatex/xelatex)

Prima di eseguire: aggiornare scripts/ e geant4_src/ con le versioni
correnti dei file sorgente (cp ../*.py scripts/; cp ../geant4_cnao_ring/src/*.cc
../geant4_cnao_ring/include/*.hh ../geant4_cnao_ring/*.cc
../geant4_cnao_ring/analysis/*.py geant4_src/), altrimenti i listing
punteranno a versioni vecchie.
"""

import re
import os
import subprocess
import sys
import tempfile

BASE = os.path.dirname(os.path.abspath(__file__))
DOCX = os.path.join(BASE, '..', 'resoconto_progetto_ESTESO.docx')
OUT_TEX = os.path.join(BASE, 'main.tex')

SYMBOL_SUBS = [
    ('≥', '$\\geq$'),
    ('≤', '$\\leq$'),
    ('≈', '$\\approx$'),
    ('→', '$\\rightarrow$'),
    ('←', '$\\leftarrow$'),
    ('×', '$\\times$'),
    ('·', '$\\cdot$'),
    ('−', '-'),
    ('²', '\\textsuperscript{2}'),
    ('³', '\\textsuperscript{3}'),
    ('¹', '\\textsuperscript{1}'),
    ('⁴', '\\textsuperscript{4}'),
    ('⁵', '\\textsuperscript{5}'),
    ('⁷', '\\textsuperscript{7}'),
    ('⁸', '\\textsuperscript{8}'),
    ('⁻', '-'),
    ('µ', '$\\mu$'),
    ('ξ', '$\\xi$'),
    ('τ', '$\\tau$'),
]

LISTINGS_SETUP = r'''
\usepackage{listings}
\lstdefinelanguage{Python}{
  keywords={def,class,import,from,as,return,if,elif,else,for,while,in,not,and,or,is,None,True,False,with,try,except,finally,lambda,yield,pass,break,continue,global,nonlocal,assert,del,raise},
  keywordstyle=\color{blue}\bfseries,
  sensitive=true,
  comment=[l]{\#},
  commentstyle=\color{gray}\itshape,
  stringstyle=\color{teal},
  morestring=[b]',
  morestring=[b]"
}
\lstset{
  basicstyle=\ttfamily\scriptsize,
  breaklines=true,
  breakatwhitespace=false,
  frame=single,
  numbers=left,
  numberstyle=\tiny\color{gray},
  showstringspaces=false,
  tabsize=4,
  columns=fullflexible,
  captionpos=b
}

\begin{document}'''


def run_pandoc(tmp_tex):
    cmd = [
        'pandoc', DOCX, '-o', tmp_tex,
        '--extract-media=media', '--wrap=preserve', '--standalone',
        '--toc', '--toc-depth=2',
        '-M', 'title=Dalla comunicazione tra due neuroni a una politica decisionale appresa',
        "-M", "author=Andrea Corti --- Universita' di Pavia (TLB / Fisica)",
        '-M', 'date=Settembre 2026',
        '-M', 'lang=it',
        '--pdf-engine=xelatex',
    ]
    subprocess.run(cmd, cwd=BASE, check=True)


def main():
    with tempfile.NamedTemporaryFile(suffix='.tex', dir=BASE, delete=False) as tf:
        tmp_tex = tf.name
    try:
        run_pandoc(os.path.basename(tmp_tex))
        with open(tmp_tex, encoding='utf-8') as f:
            text = f.read()
    finally:
        os.remove(tmp_tex)

    # media/ viene estratto gia' relativo a BASE (niente prefisso da correggere
    # quando si esegue pandoc con cwd=BASE)

    for ch, repl in SYMBOL_SUBS:
        text = text.replace(ch, repl)

    # Refuso pre-esistente nel docx sorgente (glitch di digitazione: lettere
    # cirilliche omoglife al posto di quelle latine in "costruire")
    text = text.replace('костruire', 'costruire')

    name_to_path = {}
    for d, relprefix in [(os.path.join(BASE, 'scripts'), 'scripts'),
                          (os.path.join(BASE, 'geant4_src'), 'geant4_src')]:
        for fn in os.listdir(d):
            name_to_path[fn] = f'{relprefix}/{fn}'

    def unescape_for_lookup(raw):
        return raw.replace('\\_', '_').replace('\\.', '.')

    def lang_for(fn):
        if fn.endswith('.py'):
            return 'Python'
        if fn.endswith(('.cc', '.hh', '.cpp', '.h')):
            return 'C++'
        return 'Python'

    pattern = re.compile(r'\\emph\{\\textbf\{Codice sorgente: (.*?)\}\}')
    matches = list(pattern.finditer(text))

    out = []
    cursor = 0
    n_listings = 0
    n_missing = 0
    for idx, m in enumerate(matches):
        out.append(text[cursor:m.start()])
        raw_caption_escaped = m.group(1)
        lookup_caption = unescape_for_lookup(raw_caption_escaped)
        filenames = re.findall(r'[\w][\w\-]*\.(?:py|cc|hh|cpp|h)', lookup_caption)

        search_start = m.end()
        next_marker_pos = matches[idx + 1].start() if idx + 1 < len(matches) else len(text)
        next_heading = re.search(r'\\(sub)?section\{', text[search_start:next_marker_pos])
        block_end = search_start + next_heading.start() if next_heading else next_marker_pos

        out.append(f'\\textit{{\\textbf{{Codice sorgente: {raw_caption_escaped}}}}}\n\n')
        for fn in filenames:
            rel = name_to_path.get(fn)
            if rel is None:
                out.append(f'% ATTENZIONE: file non trovato per la didascalia: {fn}\n')
                n_missing += 1
                continue
            out.append(f'\\lstinputlisting[language={lang_for(fn)}]{{{rel}}}\n')
            n_listings += 1
        cursor = block_end

    out.append(text[cursor:])
    result = ''.join(out)
    result = result.replace('\n\\begin{document}', LISTINGS_SETUP, 1)

    with open(OUT_TEX, 'w', encoding='utf-8') as f:
        f.write(result)

    print(f'{len(matches)} marker di codice trovati, {n_listings} listing inseriti, '
          f'{n_missing} file non trovati (controllare)')
    print(f'Scritto: {OUT_TEX}')
    if n_missing:
        sys.exit(1)


if __name__ == '__main__':
    main()
