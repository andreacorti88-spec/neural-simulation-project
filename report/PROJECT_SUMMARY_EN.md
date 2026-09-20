# neuroni-progetto — English project summary

**Author:** Andrea Corti — University of Pavia (TLB / Physics)
**Repository:** [neural-simulation-project](https://github.com/andreacorti88-spec/neural-simulation-project)
**Period:** covered through September 2026

*This is a synthesis, not a translation of the full report. The complete Italian report (130+ pages, every experiment explained in detail, full source code, equations) is [`resoconto_progetto_ESTESO.pdf`](resoconto_progetto_ESTESO.pdf) / [`.docx`](resoconto_progetto_ESTESO.docx). This document exists so the project can be read and evaluated by an English-speaking audience — in particular, colleagues at CNAO and at NIRS/HIMAC (Chiba) — without requiring Italian.*

## What this project is

A self-directed, incrementally-built computational neuroscience project spanning single-neuron communication up to a 40-million-neuron network, and — in its most recent and most relevant thread for this audience — a direct integration between clinical hadrontherapy beam physics (simulated in Geant4, using CNAO's actual clinical beam parameters) and a spiking neural network model.

The unifying discipline across the whole project, stated explicitly and followed consistently: **every result is verified by actually running the code**, surprising findings are re-tested across multiple independent random seeds before being trusted, and negative or null results are documented with the same care as positive ones. Several sections of the report exist specifically to record a mistake, an artifact, or a limitation that was found and then corrected — this is treated as equally valuable as a clean positive result, not as something to hide.

## Part 1 — The neuroscience arc (brief)

Progressive build-up, each stage verified before moving to the next:

1. **Basic communication** between two, then many, spiking neurons (Izhikevich model, Brian2): synaptic threshold, short-term plasticity, signal propagation, temporal summation.
2. **Learning**: spike-timing-dependent plasticity (STDP), pattern classification.
3. **Scale**: a balanced excitatory/inhibitory network pushed up to 40 million neurons (57% of mouse-brain scale), including a documented thermal-throttling instability found and worked around at full scale.
4. **Spatial structure**: small-world connectivity, a 3D brain-like shape, interactive visualizations.
5. **Attractor dynamics**: a verified ring attractor (persistence of activity, multistability, robustness to noise), chains of linked attractors.
6. **Decision-making**: competition between alternatives, reinforcement learning, a contextual decision policy.
7. **Final synthesis**: spatial structure combined with biophysically realistic synapses and a synaptic density matching the most advanced published model to date (Kuriyama, Akira et al. 2025).
8. **Learned bidirectional communication in the spiking domain**: an initially completely failed attempt at two ring-attractor populations learning to communicate, root-caused and fixed with three specific mechanisms (a transmission switch, homeostatic normalization, an emergency suppression circuit), verified across 8 independent seeds, then extended to multiple simultaneous messages and to a three-population relay network.

This arc is the substrate on which the CNAO/Geant4 work (Part 2) is built: the "ring attractor" architecture used throughout Part 2 is the same two-population spiking network from step 8, with the same 100 excitatory neurons per ring, indexed at the same angular positions `xi = i/100`.

## Part 2 — Physical radiation impact on the neural architecture (Geant4 + CNAO)

This is the thread most directly relevant to a CNAO/HIMAC audience, and the part of the project built specifically with the goal of connecting computational neuroscience to actual clinical hadrontherapy physics.

### Setup

A Geant4 (v11.4.0) C++ project, `cnaoRingImpact`, models a ring of 100 tissue-equivalent spheres (10 µm radius, typical cortical soma size) arranged on a 500 µm-radius ring — the physical counterpart of the 100-excitatory-neuron ring attractor in the Brian2 model, at matching angular positions. A clinical-scale proton or carbon-12 pencil beam, using CNAO's actual published beam parameters (protons 70/150/225 MeV; carbon-12 120/260/400 MeV/u), is fired through this geometry using **QBBC**, Geant4's reference physics list for clinical hadron energies (explicitly *not* Geant4-DNA, which is only valid up to ~1 MeV and has no validated ion-transport model at clinical energies — this distinction, and why it matters, is documented in detail in the repository).

Per-neuron physical dose (deposited energy) and, in later sections, dose-averaged LET are recorded for every run.

### Key findings, in the order they were established

1. **Bragg peak validated against literature.** A companion tool, `braggProfile`, scans dose vs. depth in a sliced water phantom and locates the true Bragg peak for each clinical energy, validated within 1% against published reference values (150 MeV protons: peak found at 156.5 mm, literature ~157.7 mm).

2. **A real methodological artifact found and corrected.** The first plateau-vs-peak dose comparison gave a physically backwards result (peak dose *lower* than plateau dose). Cross-checked against the independent `braggProfile` depth-scan data, which showed the physics itself was correct, the problem was isolated to the ring geometry: the original lateral beam-sampling disk (550 µm) was too narrow at depth, where multiple Coulomb scattering broadens the beam beyond it. Fixed by widening the sampling disk to 5 mm (a realistic clinical spot sigma). All six plateau-vs-peak configurations then showed the physically correct result (peak dose higher than plateau, by 1.4x to 9.6x depending on configuration) — the defining physical signature of hadrontherapy's clinical advantage.

3. **A functional knockout experiment linking physics to network function.** The 10 highest-physical-dose neurons from a given Geant4 run are "silenced" in the Brian2 network (all outgoing synapses removed — the most direct functional model of a neuron that can no longer influence anything downstream, regardless of what it receives). Applied to a two-message communication task (messages at ring positions x=0.15 and x=0.65), the message whose encoding neurons overlap the silenced set collapses while the other remains intact — verified across 8 independent seeds, and cross-checked against a random-knockout negative control that produces a much weaker, statistically marginal effect.

4. **A clean dose-response gradient across four independent physical scenarios.** Repeating the knockout experiment with dose maps from different particles, energies, and depths (each producing a different physical set of "highest-dose" neurons) revealed a consistent relationship: effect strength scales with how unbalanced the silenced-neuron overlap is between the two messages' own encoding zones — from a 4-vs-0 neuron overlap (strongest effect, 7.07x spike-count ratio, protons at the Bragg peak, consistent across all 8/8 seeds) down to a near-1-vs-1 overlap (statistically null effect, carbon-12 at the peak).

5. **Functional recovery under extended training — a compensatory-plasticity finding.** After the most severe knockout tested, training the network for 5x longer than the original protocol shows the damaged message's response growing progressively and substantially (roughly +290% between the first and second half of the extended training, verified across 4 seeds) — without ever fully catching up to the undamaged message. Interpreted cautiously as a simplified computational analogue of post-lesion compensatory plasticity: the silenced neurons never recover (structurally permanent, as in real radiation damage), but the surviving plasticity reorganizes to partially route around the damage.

6. **A discovered point of no return.** Pushing the knockout to its logical extreme — silencing the *entire* stimulus-encoding zone of a message (11 neurons, rather than a dose-derived subset of 10) — makes the recovery from finding 5 disappear completely, confirmed across 5 independent runs (1 single run + 4 Monte Carlo seeds). Mechanistically sensible: with the whole encoding population's outputs gone, there is no longer any surviving neuron that both receives that message's stimulus *and* can still influence anyone else — compensatory plasticity requires at least one intact exit path, not just more training time.

7. **An RBE-LET biological-dose refinement, tied directly to the methodology pioneered at NIRS/HIMAC (Chiba).** Physical dose (deposited energy) is not the clinically relevant quantity in carbon-ion therapy — biological effectiveness depends strongly on LET (linear energy transfer), which is exactly the problem NIRS addressed with the Kanai mixed-beam model, still the basis of clinical dose specification at HIMAC and, in a similar spirit, at CNAO (clinical RBE = 3.0 at the "neutron-equivalent point," 80 keV/µm dose-averaged LET). Dose-averaged LET tracking was added to the Geant4 code, and a simplified RBE(LET) conversion — explicitly *not* the full clinical microdosimetric kinetic model, which would need cell-survival parameters not available here, but a linear approximation anchored to the same NIRS reference point — was applied to recompute the "highest-dose" neuron ranking in terms of biological rather than physical dose. Result: carbon-12 at its Bragg peak reaches RBE ≈ 3.3, protons at their peak stay at RBE ≈ 1.15 (both consistent with the literature), but *within* a single fixed depth the neuron ranking barely changes (9/10 overlap for both particles) — because LET depends mainly on depth (residual range), not on lateral position across the thin ring. The existing dose-response gradient (finding 4) is therefore confirmed as robust to this refinement, not undermined by it.

8. **A limitation found, documented, and then corrected — arguably the most important methodological result of this thread.** Investigating why repeated Geant4 runs of the identical physical configuration gave different results, it emerged that `cnaoRingImpact` seeded its random generator from the system clock, with no way to reproduce a run. Adding an explicit seed parameter and repeating the same carbon-12-at-peak configuration across 4 independent seeds showed that the "top-10 highest-dose neuron" ranking is essentially uncorrelated between runs (average pairwise overlap 1.0 out of 10) — dominated by Poisson counting noise, since each 10 µm neuron intercepts on average only 1-2 direct primary particles out of 500,000 simulated events. The ring's *total* deposited dose is stable run-to-run (±15%); only its distribution across the 100 individual neurons is noise-dominated. This does not invalidate the causal *mechanism* found in finding 3 (independently verified across dozens of Brian2 Monte Carlo seeds each time), but it does mean the specific neuron-level correspondence claimed in earlier sections to "the CNAO beam's dose map" was overstated. Rebuilding the knockout set from a 4-seed-averaged dose map (rather than a single noisy run) and re-running the full 8-seed verification showed a 1-vs-0 neuron overlap producing an effect statistically indistinguishable from the zero-overlap negative control — revealing that a minimum overlap imbalance of roughly 2 neurons is needed before the physical signal rises above the network's own baseline noise floor. The dose-response gradient from finding 4 is thereby *refined*, not invalidated, by this honestly-reported limitation.

9. **The real Kase/NIRS microdosimetric kinetic model (MKM), validated against published measurements before being trusted.** Finding 7 used a simplified linear RBE(LET) approximation, explicitly flagged as not the actual clinical model. This finding replaces it with the real saturation-corrected MKM (Hawkins; Kase et al. 2006, *Radiat Res* 166:629-638; Kase et al. 2011, *J Radiat Res* 52:59-68 — the same 2011 paper that describes the TEPC-based biological-dose measurement used routinely for quality assurance at HIMAC), implemented with the real published HSG-cell parameters (α₀=0.13 Gy⁻¹, β=0.05 Gy⁻², domain radius rd=0.42 µm, saturation parameter y₀=150 keV/µm, reference dose D10,R=5.0 Gy). Two simplifications remain, both explicitly stated: dose-averaged LET is used as a proxy for the dose-mean lineal energy y_D that would normally come from a measured TEPC spectrum or a full track-structure microdosimetric simulation, and the saturation correction is evaluated in its single-value (delta-spectrum) closed form rather than integrated over a full lineal-energy spectrum. Critically, **the implementation was validated against Kase et al.'s own published data before any result was trusted**: applying the model to their reported y* values for a 290 MeV/u carbon-ion SOBP beam reproduces their published RBE10 within 10-15% (e.g. y*=58.7 keV/µm → RBE10 calculated 2.23, vs. ~2.0-2.3 read from their figure; at the NIRS clinical reference point, LET_d=80 keV/µm → RBE10 calculated 2.67, vs. the clinical reference value of 3.0). Applied to this project's own carbon-12 and proton dose maps: RBE10 reaches 3.44 for carbon at the Bragg peak (physically sensible, consistent with published carbon RBE curves near the overkill region) and stays close to 1.0 for protons (as expected for a low-LET beam) — and the within-depth neuron-ranking conclusion from finding 7 holds up under the real clinical model, not just the earlier linear approximation. A reusable script (`geant4_cnao_ring/analysis/mkm_rbe.py`) implements the model with its validation check built in.

### Why this thread matters for a CNAO/HIMAC audience specifically

- It uses CNAO's actual published clinical beam parameters, not arbitrary values.
- It engages directly with the same LET/RBE methodology (the Kanai mixed-beam model) that NIRS pioneered at HIMAC and that underlies carbon-ion clinical dose specification worldwide.
- The most scientifically mature part of this thread is not a positive result — it is finding and fixing a Monte Carlo statistics problem in the physics simulation itself, the kind of methodological literacy that matters in a dosimetry-adjacent research environment.
- The whole thread is explicit about what it is *not*: not a validated biological damage model, not the full clinical microdosimetric kinetic model, not a claim about real cortical tissue — it is a computational bridge built to reason about the *principle* (physical dose → functional damage → possible recovery), with every simplification stated rather than hidden.

## Repository map

| Location | Content |
|---|---|
| `scripts/` | All Python/Brian2 experiment scripts |
| `figures/` | Plots from every experiment |
| `visualizations/` | Interactive 3D HTML visualizations |
| `report/` | Full Italian report (PDF/DOCX) and this English summary |
| `geant4_cnao_ring/` | The Geant4 C++ project described in Part 2 above, with its own [README](../geant4_cnao_ring/README.md) ([English](../geant4_cnao_ring/README_EN.md)) |

## Honest limitations of the project as a whole

- The spiking network is a deliberately simplified ring-attractor model (Izhikevich neurons, Brian2), not a data-fitted model of real cortical tissue.
- The physical scale of the Geant4 geometry (1mm ring, 20µm soma) is a literature-based choice for a typical cortical microcircuit, not a measurement of the (scale-free, topological) Brian2 model.
- No biological damage model connects physical/biological dose to actual neuronal dysfunction — the "knockout" is a deliberately simplified functional proxy (complete removal of outgoing synapses above a dose threshold), not a validated radiobiological damage model.
- The MKM implementation used in finding 9 is the real Kase/NIRS clinical model with real HSG parameters, validated against published data — but it still substitutes dose-averaged LET for a genuinely measured or simulated lineal-energy spectrum, and uses the single-value rather than full-spectrum saturation correction; both approximations are stated explicitly and both cause the small (10-15%) validation gap against Kase et al.'s own published numbers.
- As documented in finding 8, the per-neuron physical dose map itself carries a real, quantified Monte Carlo statistics limitation that future work should address with either substantially more simulated events per run or systematic multi-seed averaging.
