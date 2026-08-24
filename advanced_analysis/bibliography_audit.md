# Bibliography audit for the PRL draft

Audit date: 2026-08-24.  Scope: all 20 entries in `arxiv_theory_paper/ref.bib` and all citation keys used by `arxiv_theory_paper/main.tex`.  No manuscript or BibTeX file was edited.

## Sources and method

The 19 DOI-backed records were queried through the Crossref DOI registry (`https://api.crossref.org/works/{DOI}`), with the DOI landing page retained as the primary-publication link below.  The arXiv record was checked on the official abstract page and HTML full text for arXiv:2607.07228v1.  Metadata were compared field by field: author identity/order, title, journal, year, volume, issue, page or article number, DOI, and arXiv identifier/version.  Case, punctuation, journal abbreviations, and initials were not treated as errors when the underlying identity was unambiguous.

Verdicts use the following meaning: **Verified** = no substantive metadata correction is required; **Check suggested** = metadata are usable but a canonical spelling, article-number convention, or citation-scope note is advisable; **Needs fix** = a substantive mismatch requiring correction; **Unverifiable** = no primary metadata could be established.

## Per-entry verification

| BibTeX key | DOI/arXiv primary record | Metadata comparison | Verdict |
|---|---|---|---|
| `Hsu2016BIC` | [10.1038/natrevmats.2016.48](https://doi.org/10.1038/natrevmats.2016.48) | Nature Reviews Materials **1**, 16048 (2016); authors and title agree. Crossref spells the last author `Soljačić`; the entry's `Soljacic, Marin` is an unaccented spelling. | Check suggested |
| `Koshelev2023Review` | [10.3367/UFNe.2021.12.039120](https://doi.org/10.3367/UFNe.2021.12.039120) | Physics-Uspekhi **66**, 494–517 (2023); title, author order, pages, and DOI agree. Crossref includes additional middle initials and issue 05, which are optional for this style. | Verified |
| `Azzam2021Review` | [10.1002/adom.202001469](https://doi.org/10.1002/adom.202001469) | Advanced Optical Materials **9**, article 2001469 (2021); authors, title, article number, and DOI agree. | Verified |
| `Yang2014Analytical` | [10.1103/PhysRevLett.113.037401](https://doi.org/10.1103/PhysRevLett.113.037401) | Physical Review Letters **113**, 037401 (2014); authors, title, article number, and DOI agree. | Verified |
| `Zhen2014Topological` | [10.1103/PhysRevLett.113.257401](https://doi.org/10.1103/PhysRevLett.113.257401) | Physical Review Letters **113**, 257401 (2014); authors, title, article number, and DOI agree. Crossref gives `Soljačić` rather than `Soljacic`. | Check suggested |
| `Hsu2013Observation` | [10.1038/nature12289](https://doi.org/10.1038/nature12289) | Nature **499**, 188–191 (2013); authors, title, pages, year, and DOI agree. | Verified |
| `Sadrieva2017LeakyResonances` | [10.1021/acsphotonics.6b00860](https://doi.org/10.1021/acsphotonics.6b00860) | ACS Photonics **4**, 723–727 (2017); all nine authors, title, pages, year, and DOI agree. | Verified |
| `Abdrabou2022FrequencyPerturbation` | [10.1103/PhysRevA.106.013523](https://doi.org/10.1103/PhysRevA.106.013523) | Physical Review A **106**, 013523 (2022); authors, title, article number, year, and DOI agree. | Verified |
| `Rayleigh1907` | [10.1098/rspa.1907.0051](https://doi.org/10.1098/rspa.1907.0051) | Proceedings of the Royal Society A **79**, 399–416 (1907); title, volume, pages, year, and DOI agree. Crossref stores the civil name `John William Strutt`, whereas the established publication name is Lord Rayleigh; the current `Rayleigh, Lord` is a conventional citation form, not a wrong paper. | Check suggested |
| `Wojcik2021ThresholdScattering` | [10.1103/PhysRevLett.127.277401](https://doi.org/10.1103/PhysRevLett.127.277401) | Physical Review Letters **127**, 277401 (2021); authors, title, article number, year, and DOI agree. | Verified |
| `Assouar2018AcousticMetasurfaces` | [10.1038/s41578-018-0061-4](https://doi.org/10.1038/s41578-018-0061-4) | Nature Reviews Materials **3**, 460–472 (2018); six authors, title, pages, year, and DOI agree. | Verified |
| `Ni2019Splitter` | [10.1103/PhysRevB.100.104104](https://doi.org/10.1103/PhysRevB.100.104104) | Physical Review B **100**, 104104 (2019); authors, title, article number, year, and DOI agree. | Verified |
| `Torrent2018Anomalous` | [10.1103/PhysRevB.98.060101](https://doi.org/10.1103/PhysRevB.98.060101) | Physical Review B **98**, 060101 (2018); title, author, volume, year, and DOI agree. Crossref's article-number field is `060101`; APS commonly displays the Rapid Communication designation `060101(R)`. The current `pages={060101(R)}` is therefore a style choice, not a DOI mismatch. | Check suggested |
| `Jin2019Grating` | [10.1103/PhysRevApplied.11.011004](https://doi.org/10.1103/PhysRevApplied.11.011004) | Physical Review Applied **11**, 011004 (2019); authors, title, article number, year, and DOI agree. | Verified |
| `Bernard2022NegativeReflection` | [10.1103/PhysRevApplied.17.024059](https://doi.org/10.1103/PhysRevApplied.17.024059) | Physical Review Applied **17**, 024059 (2022); five authors, title, article number, year, and DOI agree. Crossref gives the accented surname `Léon`; the entry uses `Leon`. | Check suggested |
| `Fan2021Multifunctional` | [10.1103/PhysRevApplied.16.044029](https://doi.org/10.1103/PhysRevApplied.16.044029) | Physical Review Applied **16**, 044029 (2021); authors, title, article number, year, and DOI agree. | Verified |
| `Cao2024UnderwaterAbnormal` | [10.1103/PhysRevApplied.21.034015](https://doi.org/10.1103/PhysRevApplied.21.034015) | Physical Review Applied **21**, 034015 (2024); six authors, title, article number, year, and DOI agree. | Verified |
| `Farhat2024UltrahighQ` | [10.1002/advs.202402917](https://doi.org/10.1002/advs.202402917) | Advanced Science **11**, article 2402917 (2024); six authors, title, article number, year, and DOI agree. Crossref uses `Martínez`; the entry's initials are unambiguous but omit the accent. | Check suggested |
| `Yang2024MergedFP` | [10.1103/PhysRevB.110.184108](https://doi.org/10.1103/PhysRevB.110.184108) | Physical Review B **110**, 184108 (2024); authors, title, article number, year, and DOI agree. | Verified |
| `Karavaev2026Rayleigh` | [arXiv:2607.07228v1](https://arxiv.org/abs/2607.07228) | Official arXiv record: *Rayleigh Bound States in the Continuum*, Ilya Karavaev, Mingzhao Song, Andrey Bogdanov; submitted 2026-07-08, v1, category physics.optics. The current eprint, year, archive prefix, class, and v1 note agree. The abstract and full text explicitly discuss Rayleigh anomalies, poles/eigenfrequencies, and the collision of a leaky pole with a Rayleigh branch point. | Verified |

Summary: 14 entries are fully verified as written; 6 have only canonical-spelling or publisher-style checks suggested; no substantive title/DOI/year/volume/page mismatch and no unverifiable entry was found. The count treats `Rayleigh1907`, `Torrent2018Anomalous`, and the three diacritic cases as usable but worth normalizing.

## Exact optional BibTeX normalizations

These are not applied here because the requested writable scope is audit-only. They improve canonical spelling or make the source convention explicit.

1. Preserve the established `Lord Rayleigh` display while recording the registry name only if a fully canonical Crossref form is desired:

```bibtex
author = {Strutt, John William},
```

The current `author = {Rayleigh, Lord}` is acceptable for this historical paper and should not be changed merely to satisfy Crossref.

2. Restore diacritics in author fields (requires the manuscript's UTF-8 toolchain, or equivalent TeX accent commands):

```bibtex
author = {Hsu, Chia Wei and Zhen, B. and Stone, A. D. and Joannopoulos, J. D. and Soljačić, Marin},
author = {Zhen, B. and Hsu, C. W. and Lu, L. and Stone, A. D. and Soljačić, Marin},
author = {Bernard, Simon and Chikh-Bled, Feriel and Kourchi, Hasna and Chati, Farid and Léon, Fernand},
author = {Farhat, M. and Achaoui, M. and Martínez, J. A. I. and Addouche, M. and Wu, Y. and Khelif, A.},
```

If accent-safe source is preferred, use TeX forms such as `Solja\v{c}i\'c`, `L\'eon`, and `Mart\'{\i}nez`; the exact encoding should be chosen consistently with the compiler. The current initials are still bibliographically identifiable.

3. For a strict Crossref article-number representation of the Rapid Communication, use:

```bibtex
pages = {060101},
note  = {Rapid Communication},
```

The current `060101(R)` follows common APS display practice and is not a required correction.

4. For a self-contained arXiv record, the following optional fields make the source/version explicit:

```bibtex
url  = {https://arxiv.org/abs/2607.07228},
note = {arXiv:2607.07228v1 [physics.optics]},
```

## Citation-support audit

The following checks concern whether the cited source supports the sentence in which it is used, independent of metadata correctness.

| `main.tex` location | Citation group | Support judgment |
|---|---|---|
| Introduction, BIC definition and photonic context | `Hsu2016BIC`, `Koshelev2023Review`, `Azzam2021Review`, `Yang2014Analytical`, `Zhen2014Topological` | Collectively appropriate: reviews support the definition/context; Yang supports analytical BIC treatment; Zhen supports the topological-BIC context. |
| Introduction, threshold bookkeeping | `Rayleigh1907`, `Wojcik2021ThresholdScattering` | Appropriate as a two-source combination: Rayleigh supplies the grating/diffraction-threshold foundation, and Wojcik supplies the scattering-matrix singularity near channel openings. The acoustic extension is the manuscript's own application, not a result directly established by Wojcik. |
| Introduction, acoustic metagratings | `Assouar2018AcousticMetasurfaces`, `Ni2019Splitter`, `Torrent2018Anomalous`, `Jin2019Grating` | Appropriate. The review supports the general metasurface statement; the three papers support compact acoustic-grating functions and diffraction engineering. |
| Introduction, waterborne anomalous/negative reflection | `Bernard2022NegativeReflection`, `Fan2021Multifunctional`, `Cao2024UnderwaterAbnormal` | Appropriate. These papers directly concern waterborne/underwater acoustic metagratings and anomalous or negative diffraction/reflection. |
| Introduction, underwater high-Q/merged BIC behavior | `Farhat2024UltrahighQ`, `Yang2024MergedFP` | Appropriate: Farhat is ultrasound BIC/high-Q evidence; Yang is an underwater phononic-crystal merged-BIC result. |
| Introduction, conventional BIC pole behavior | `Hsu2013Observation`, `Sadrieva2017LeakyResonances`, `Abdrabou2022FrequencyPerturbation` | Sadrieva and Abdrabou support leaky-pole/perturbative behavior. Hsu2013 is a direct observation of trapped light and is historically relevant, but it is not by itself a source for the general phrase “pole trajectories”; consider moving it to a trapped-light sentence or retaining it only as historical context. |
| Introduction, electromagnetic Rayleigh-BIC precedent | `Karavaev2026Rayleigh` | Directly supported. The official full text states that the BIC is identified by a real-axis scattering pole with simultaneous numerator zero and explicitly describes a leaky-pole/Rayleigh-branch-point collision. |

No citation key used by `main.tex` is absent from this report. `suppl.tex` contains no `\cite{...}` key at audit time.

## Recommended action for the primary thread

No urgent BibTeX repair is required for DOI resolvability or bibliographic identity. Before submission, normalize the four author diacritics if the chosen LaTeX engine supports them, decide whether to use the Crossref article number `060101` or the APS Rapid Communication form `060101(R)`, and tighten the Hsu2013 sentence so that its citation is not presented as the sole support for pole trajectories. The arXiv citation is valid as a 2026-v1 preprint; it should be clearly labeled as a preprint rather than a peer-reviewed journal article.
