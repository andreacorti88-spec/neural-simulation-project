# cnaoRingImpact — physical impact of CNAO ion beams on a ring-attractor architecture

*(Italian version: [`README.md`](README.md).)*

**Integral part of `neuroni-progetto`** (not a satellite project): lives
in `neuroni-progetto/geant4_cnao_ring/`, same Git repository as the
Brian2 code. Two languages (C++/Geant4 for the physics, Python/Brian2
for the spiking network) and two separate build toolchains for
technical necessity, but one project, one history, one report.

Computes the dose and direct-crossing count deposited by clinical CNAO
ion beams (protons and carbon-12) on a physical geometry representing
the spiking ring attractor of the same project (report sections
5.6-5.20): a ring of 100 tissue-equivalent spheres, one per excitatory
neuron of the network, positioned at exactly the same angular
coordinates (`xi = i/100`) used in the Brian2 model. The link is
concrete, not just thematic: the "Link to Brian2" section below uses
the dose maps computed here directly to silence specific neurons in
`dialogo_spiking_cnao_knockout.py`, in the project's main directory.

Sibling project to `geant4dna_icsd` (electron-scale nanodosimetry,
dedicated to verifying the Villagrasa/Baiocco 2026 paper), but with a
completely different scope and physics — see below.

## Why QBBC, not Geant4-DNA

`geant4dna_icsd` uses Geant4-DNA physics (options 2/4/6), valid **only
up to ~1 MeV and only for electrons** in liquid water (see that
project's `PhysicsList.cc`, `SetMaxEnergy(1*MeV)`). Clinical CNAO
beams are 60-250 MeV protons and 120-400 MeV/u carbon-12 ions: tens to
hundreds of MeV per nucleon, an energy domain and transport regime
(including carbon fragmentation nuclear interactions) for which
Geant4-DNA has no validated model. This project uses **QBBC**,
Geant4's reference physics list, explicitly recommended (and used in
the official "Hadrontherapy" example) for protons and ions at clinical
energies.

The two physics models are neither interchangeable nor directly
composable in a single project without serious multi-scale work (see
"Future developments" below) — the same distinction discussed in the
`Analisi_critica_nanodosimetria.pdf` document, where nanodosimetry
(Geant4-DNA) and macroscopic primary-beam transport are explicitly two
separate levels of the same physical chain.

## Geometry

- **World**: tissue-equivalent water block, 4mm x 4mm x 4mm.
- **Ring**: 100 spheres ("Neuron"), 10um radius (typical cortical
  soma, ~20um diameter), arranged on a circle of 500um radius (~1mm
  diameter ring) in the Z=0 plane, at angular position `xi = i/100`,
  `i=0..99` — the same indexing used for message positions in the
  spiking dialogues (e.g. x=0.15, x=0.65).
- **Beam**: pencil beam along +Z, lateral position (x,y) sampled
  uniformly on a disk large enough to cover realistic clinical spot
  spread (5mm radius, after a correction discussed under "Limits"
  below). A real clinical CNAO spot has a sigma of a few mm, much
  wider than the ring itself: at this microscopic scale the beam is
  therefore effectively uniform, and any dose difference observed
  between positions in the results is genuine microdosimetric
  stochasticity (real statistical fluctuation in the number and type
  of interactions), not a beam-shape artifact.

**Stated limit**: the physical scale (1mm ring, 20um soma) is a
literature-based choice for a typical cortical microcircuit, not a
measurement of the network simulated in Brian2 (which is a topological
model, not a metric one — it never had its own physical scale). This
is a conceptual bridge between two separately validated models
(clinical beam physics on one side, ring-attractor network dynamics on
the other), not a claim that real cortical neurons sit on a literal
CNAO beam path.

## What is recorded

For each neuron, over the whole duration of a run (accumulated across
all simulated primary events):
- **total deposited energy** (all tracks, primary and secondary,
  including carbon nuclear interactions)
- **number of events in which the PRIMARY track** (proton or ion,
  `ParentID==0`) physically crosses that neuron
- **dose-averaged LET** (added in report section 5.31, same quantity
  used in the Kanai mixed-beam model for clinical biological-dose
  calculation at NIRS/HIMAC and CNAO)

Output: one CSV per run (`neuron_dose_run0.csv`), renamed per
configuration by `analysis/analyze_ring_dose.py`.

## Results: plateau vs true Bragg peak (`results_picco/`)

With the corrected geometry (5mm sampling disk, a real clinical spot
sigma) and the ring positioned EXACTLY at the Bragg peak for each
energy (found with `braggProfile/`, validated within 1% against
literature values):

| Configuration | Peak depth | Plateau dose (MeV) | Peak dose (MeV) | Ratio |
|---|---|---|---|---|
| 70 MeV proton | 39.5mm | 0.118 | 0.953 | 8.11x |
| 150 MeV proton | 156.5mm | 0.099 | 0.499 | 5.06x |
| 225 MeV proton | 315.5mm | 0.104 | 0.996 | 9.56x |
| 120 MeV/u carbon-12 | 34.5mm | 7.068 | 20.996 | 2.97x |
| 260 MeV/u carbon-12 | 135.5mm | 3.504 | 11.258 | 3.21x |
| 400 MeV/u carbon-12 | 275.5mm | 3.254 | 4.412 | 1.36x |

All six configurations show a physically higher peak than plateau, as
expected — the defining signature of hadrontherapy's clinical
advantage (concentrating dose on the target, sparing healthy tissue
traversed beforehand). Rigor note: the first attempt for 225 MeV
protons (the deepest case, with the narrowest peak) gave, at 50,000
events, an anomalous 0.24x ratio (peak LOWER than plateau) — a pure
Monte Carlo statistics artifact from insufficient events on a very
narrow peak at large depth, not a systematic error: resolved with
500,000 events (ratio 9.56x, consistent with the other five cases).
Lesson kept in mind for any future extension to even greater depths.

Note: an earlier version of this table used the original 550um
sampling disk, later found to be a methodological artifact at depth
(see "Limits" below and report section 5.25) — multiple Coulomb
scattering dilutes the primary beam fluence within too narrow a disk
the deeper the beam travels. The table above already reflects the
corrected 5mm-disk geometry.

## Limits of this first level (and what was found and fixed)

- **Low per-neuron statistics — a real, load-bearing limit**: with a
  ~0.03% geometric direct-crossing probability per neuron per event,
  a few hundred thousand primaries give only ~1-2 expected direct hits
  per neuron. This is not just noise in the point-by-point map — it
  means the "top-10 highest-dose neurons" ranking from a SINGLE Geant4
  run is dominated by Poisson counting noise, not a genuine
  deterministic spatial pattern (report section 5.32): repeating the
  identical physical configuration across 4 independent seeds gave an
  average pairwise top-10 overlap of only 1.0/10. The ring-TOTAL dose
  stays stable (±15% run to run), only its per-neuron distribution
  does not. Partially addressed by averaging the dose map over
  multiple seeds (report section 5.33) rather than trusting a single
  run; a fully clean point-by-point map would still need 2-3 orders of
  magnitude more events per run.
- **No biological damage model**: this project computes only deposited
  energy (physical dose), not a biological endpoint (DNA breaks, cell
  death) — the same limit explicitly stated in the Villagrasa/Baiocco
  paper for the dose→damage→survival step (see
  `Analisi_critica_nanodosimetria.pdf`, section 5). A first, simplified
  biological-weighting step (RBE via dose-averaged LET, anchored to the
  NIRS/Kanai clinical reference point) was added in section 5.31 — see
  below — but it is explicitly NOT the full clinical microdosimetric
  kinetic model (MKM, Kase/Inaniwa), which would require survival-curve
  parameters not available here.
- **Bragg peak depth — RESOLVED**: `braggProfile/` (separate
  sub-project) finds the true peak for a given energy by scanning dose
  vs. depth in a sliced phantom (validated within 1% against known
  reference values: 150 MeV protons → peak at 156.5mm, literature
  ~157.7mm). `cnaoRingImpact` accepts ring depth as a parameter. During
  the plateau-vs-peak comparison a real artifact emerged: the original
  550um lateral sampling disk is too narrow at greater depths, where
  multiple beam scattering dilutes primary fluence in that region —
  giving LOWER dose at the peak than at plateau, the opposite of the
  expected physics. Diagnosed by cross-checking against independent
  `braggProfile` data (same energy, 5x5cm-wide slab, insensitive to the
  artifact) and fixed by widening the disk to 5mm (a real CNAO clinical
  spot sigma, always much wider than multiple scattering even several
  cm deep). After the fix: peak/plateau ratio 2.81x for 70 MeV protons,
  in the physically correct direction.
- **No direct (coupled) link with the Brian2 simulation**: this is a
  pure physics output (per-position dose CSV), analyzed separately from
  network function — not a real-time coupled physics-neural
  co-simulation (which would require translating dose/ionization into
  a neuronal damage/inactivation model with thresholds that, to date,
  have no consolidated value in the literature for a single neuron).

## Link to Brian2: from physical data to functional knockout

`neuroni-progetto/dialogo_spiking_cnao_knockout.py` and its later
variants take the 10 highest-dose neurons from a dose map produced
here and silence them (outgoing synapses removed) in the
multi-message spiking network (report section 5.10). Result, verified
across 8 seeds and against a negative control (a random knockout is
not enough): the message that uses the silenced neurons collapses, the
one that doesn't remains intact — a causal effect, not an artifact,
demonstrated by deliberately reversing it: with the corrected neuron
set found after discovering the sampling-disk artifact above (which by
chance included more neurons from the other message's zone), the
effect flips exactly as predicted. A clean dose-response gradient later
emerged across four independent physical scenarios (different
particle, energy, and depth each time — report sections 5.22-5.28),
and was further refined once section 5.32 revealed that individual
single-run neuron rankings were statistically noisy: rebuilding the
knockout set from a multi-seed-averaged dose map (section 5.33) showed
that an overlap imbalance of at least ~2 neurons is needed for the
effect to rise above the network's own baseline noise floor. Two
follow-up findings extend the picture further: a genuine functional
recovery under extended training after knockout (section 5.29,
computational analogue of post-lesion compensatory plasticity), and a
discovered point of no return when an entire message's encoding zone
is silenced rather than a dose-derived subset (section 5.30). A real,
literature-validated MKM (Microdosimetric Kinetic Model) RBE
calculation — the same formalism used clinically at NIRS/HIMAC — was
added in section 5.34, replacing the simplified linear RBE(LET)
approximation of section 5.31; validated against Kase et al.'s own
published measurements before being trusted, it confirms the earlier
conclusion (the dose-response gradient is robust to RBE weighting)
rather than overturning it. Full details: report, sections 5.22-5.34.

## Recent findings (sections 5.28-5.34)

- **5.28 — Gradient confirmed with protons**: Bragg-peak dose map for
  70 MeV protons, the sharpest overlap imbalance so far (4-vs-0) →
  7.07x spike-count ratio, 8/8 seeds — the strongest and most
  consistent effect in the series.
- **5.29 — Functional recovery after damage**: 5x extended training
  (300 repetitions) after the most severe knockout (5.28): the damaged
  message recovers progressively (+294% between the first and second
  half of extended training, verified across 4 seeds), without ever
  matching the undamaged message — residual STDP plasticity that
  compensates, not repairs.
- **5.30 — A point of no return**: silencing ALL 11 neurons of a
  message's encoding zone (not just a dose-derived subset) makes the
  5.29 recovery disappear completely, confirmed across 5 independent
  runs (1 single run + 4 MC seeds) — compensatory plasticity requires
  at least one surviving functional exit path.
- **5.31 — RBE-LET model (Kanai/NIRS-Chiba)**: added dose-averaged LET
  tracking per neuron (same formalism used clinically at NIRS/HIMAC).
  The RBE gap between particles is sharp (carbon ~3.3, proton ~1.15 at
  their respective peaks, consistent with the clinical RBE=3.0
  reference at 80 keV/um), but within the same depth the neuron ranking
  barely changes (9/10 overlap) — the 5.22-5.28 gradient holds up to
  this refinement.
- **5.32 — A limit discovered: the per-neuron dose map is not
  reproducible**: added an explicit seed to `cnaoRingImpact` (previously
  clock-based, different every run). Repeating the same physical
  configuration across 4 independent seeds gave an average pairwise
  top-10 overlap of only 1.0/10 — dominated by Poisson noise (~1-2
  expected direct hits per neuron out of 500k events). The ring's
  TOTAL dose stays stable (±15%); only its per-neuron distribution is
  noisy.
- **5.33 — The gradient on solid statistical ground**: rebuilt the
  knockout set from the 4-seed-AVERAGED dose map of 5.32 (1-vs-0
  overlap with the encoding zones) → 1.29x ratio, 6/8 seeds,
  statistically indistinguishable from the 0-vs-0 negative control
  (1.39x, 6/8, section 5.24). A threshold emerges: an overlap imbalance
  of AT LEAST 2 neurons is needed for the effect to rise above the
  network's baseline noise — the gradient is refined, not invalidated.
- **5.34 — The real Kase/NIRS MKM model**: replaced the linear RBE(LET)
  approximation of 5.31 with the actual saturation-corrected MKM used
  clinically at NIRS/HIMAC (Kase et al. 2011, J Radiat Res 52:59-68 —
  the same paper describing the TEPC used routinely for biological-dose
  QA at HIMAC), with real HSG cell parameters. Validated before
  publishing: applied to Kase et al.'s own published measurements, it
  reproduces their reported RBE10 within 10-15% (e.g. at the NIRS
  clinical reference point, LET_d=80 keV/um → RBE10 = 2.67 vs. the
  clinical reference 3.0). Applied to our data: RBE10 up to 3.44 for
  carbon at the peak (physically sensible), ~1.0 for the proton (as
  expected) — the 5.31 conclusion holds with the real clinical model
  too, not just the crude linear approximation. Reusable script:
  `analysis/mkm_rbe.py`.

## Future developments

1. **More statistics** (1e6-1e7 primaries) for a clean point-by-point
   map, parallelizable per configuration as already done for
   `nanoICSD` (`scan_energies.sh`) — the most direct way to shrink the
   per-neuron Poisson noise found in section 5.32.
2. **Local nanodosimetric level**: re-simulate, with Geant4-DNA, the
   detailed track ONLY in neurons that receive a direct crossing
   (multi-scale "track re-simulation" approach, the same principle used
   by TOPAS-nBio) — would connect this project to `geant4dna_icsd`
   instead of leaving them separate.
3. **Absolute neuronal damage threshold**: the knockout already done
   uses a relative criterion (top dose in a scenario), not an absolute
   threshold — which remains unconsolidated in the literature for a
   single neuron, the same limit stated in the Villagrasa/Baiocco
   paper. A sensitivity test across different criteria (mean dose
   instead of peak, knockout size) is the natural next step.
4. ~~Full clinical MKM-based RBE model~~ **DONE** (section 5.34) — see
   above; the remaining gap to the full clinical implementation is a
   measured y-spectrum (TEPC or track-structure simulation) instead of
   the LET-based single-value proxy used here.
5. ~~Knockout at the true Bragg peak~~ **DONE** (section 5.27) — see
   the dose-response gradient discussion above.

## Build and run

Same prerequisites as `nanoICSD` (Geant4 11.x), but here the standard
datasets are enough (no need for G4EMLOW/Geant4-DNA):

```bash
source /path/to/geant4-install/bin/geant4.sh
mkdir build && cd build
cmake -DGeant4_DIR=/path/to/geant4-install/lib/cmake/Geant4 ..
make -j$(nproc)

# single configuration: particle (0=proton,1=carbon-12), total kinetic energy MeV,
# ring depth mm, [explicit seed -- default: system clock, different on every run]
./cnaoRingImpact macros/run.mac 0 150            # 150 MeV proton, plateau (2mm)
./cnaoRingImpact macros/run.mac 1 3120           # 260 MeV/u carbon-12 (260*12=3120 MeV total)
./cnaoRingImpact macros/run.mac 0 70 39.5 12345  # 70 MeV proton at the peak (39.5mm), fixed seed
                                                  # -- for multi-seed MC verification, see section 5.32

python3 ../analysis/analyze_ring_dose.py results/
```
