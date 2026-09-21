#!/usr/bin/env python3
"""
Regenerates main_en.tex from resoconto_progetto_ESTESO_EN.docx (English
translation, one directory above), mirroring build_latex.py's pipeline
for the Italian version.

Pipeline:
1. pandoc converts the docx to standalone LaTeX (full preamble, table of
   contents, images extracted to media/).
2. This script replaces every source-code block (inserted in the docx as
   a one-row-per-line table in Consolas font -- pandoc turns it into
   unreadable prose with badly escaped LaTeX special characters) with a
   real \\lstinputlisting pointing at the actual source file in scripts/
   or geant4_src/ -- much cleaner and always in sync with the real code,
   not a transcription of it. Source files are shared with the Italian
   build (code and file names are identical in both languages).
3. Fixes a handful of Unicode symbols (>=, ->, µ, etc.) that have no
   glyph in the default font under XeTeX, replacing them with the LaTeX
   math-mode equivalent.
4. Adds the listings package configuration (basic Python/C++ syntax
   highlighting) right before \\begin{document}.

Usage:
    python3 build_latex_en.py         # regenerates main_en.tex
    tectonic main_en.tex              # compiles to main_en.pdf

Run build_latex.py first (or keep scripts/ and geant4_src/ in sync some
other way) since this script reuses the same source-file directories.
"""

import re
import os
import subprocess
import sys
import tempfile

BASE = os.path.dirname(os.path.abspath(__file__))
DOCX = os.path.join(BASE, '..', 'resoconto_progetto_ESTESO_EN.docx')
OUT_TEX = os.path.join(BASE, 'main_en.tex')

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
        '-M', 'title=From Communication Between Two Neurons to a Learned Decision Policy',
        '-M', "author=Andrea Corti --- University of Pavia (TLB / Physics)",
        '-M', 'date=September 2026',
        '-M', 'lang=en',
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

    # media/ is extracted already relative to BASE (shared with the Italian
    # build -- no prefix to fix since pandoc runs with cwd=BASE)

    for ch, repl in SYMBOL_SUBS:
        text = text.replace(ch, repl)

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

    pattern = re.compile(r'\\emph\{\\textbf\{Source code: (.*?)\}\}')
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
        # Bug fix (same as build_latex.py): \subsubsection{ needs TWO "sub"
        # prefixes, so (sub)? missed every \subsubsection{ heading -- (sub)*
        # matches any nesting depth.
        next_heading = re.search(r'\\(?:sub)*section\{', text[search_start:next_marker_pos])
        block_end = search_start + next_heading.start() if next_heading else next_marker_pos

        out.append(f'\\textit{{\\textbf{{Source code: {raw_caption_escaped}}}}}\n\n')
        for fn in filenames:
            rel = name_to_path.get(fn)
            if rel is None:
                out.append(f'% WARNING: file not found for caption: {fn}\n')
                n_missing += 1
                continue
            out.append(f'\\lstinputlisting[language={lang_for(fn)}]{{{rel}}}\n')
            n_listings += 1
        cursor = block_end

    out.append(text[cursor:])
    result = ''.join(out)
    result = result.replace('\n\\begin{document}', LISTINGS_SETUP, 1)

    # Same underscore-line-break fix as the Italian build (build_latex.py):
    # long filenames in the final File/Section/Content index table need an
    # explicit break point after each escaped underscore, or they overflow
    # into the neighbouring column. \lstinputlisting paths use plain
    # (non-escaped) underscores and are therefore untouched by this.
    result = result.replace('\\_', '\\_\\allowbreak{}')

    with open(OUT_TEX, 'w', encoding='utf-8') as f:
        f.write(result)

    print(f'{len(matches)} code markers found, {n_listings} listings inserted, '
          f'{n_missing} files not found (check these)')
    print(f'Written: {OUT_TEX}')
    if n_missing:
        sys.exit(1)


if __name__ == '__main__':
    main()
