# TMEM87A expression in colorectal cancer

Analysis code and derived data for *Epithelial RNA detection and tissue protein abundance of TMEM87A in colorectal cancer*, by Yin Chen, Jinglin Wang, Yashi Peng, Fei Li, Zaili Yang and Qian Yang.

Release **v1.0.0** corresponds to the first submission materials. Cite the [version-specific Zenodo record](https://doi.org/10.5281/zenodo.24535179). The companion GitHub tag is [v1.0.0](https://github.com/wangjinglin0905-gif/tmem87a-crc-expression/tree/v1.0.0).

| Supplement | Contents |
| --- | --- |
| 1 | Discovery single-cell comparisons, bulk RNA, survival models, paired CPTAC protein and matched RNA |
| 2 | GSE166555 patient pairs, epithelial compartment shares, RNA-inferred clone contrasts and derived cell measurements |
| 3 | Spatial measurements, library-size overlap, image provenance and core reproduction script |
| 4 | Laser-microdissected (PXD002137) and pathology-guided regional (PXD041412) proteomic measurements, contrasts, scripts and checks |
| 5 | Patient-level decomposition of epithelial RNA into detection and mean expression among detected cells for TMEM87A and TMEM63A |

## Reproduce selected results

Use Python 3.12 in a dedicated environment. The dependency file lists the four library versions tested for the three commands. Working from the root of this release, run:

```console
python -m pip install -r requirements.txt
python verify_release.py results
```

The entry point checks the release hashes, validates all five supplement archives and extracts them into your output directory. It checks 18 patient-paired summaries, two CPTAC targets and matched RNA associations, four spatial donors and 50 depth comparisons. The external proteomics script verifies the laser-microdissection estimates and regional patient summaries. The RNA decomposition script regenerates both tables from 110 patient–gene pairs and checks 44 external patient–gene–tissue combinations against cell-level values. The script compares regenerated tables with their archived copies. Some native text writers use platform-specific line endings; the comparison accepts either LF or CRLF.

The numbered ZIPs are the manuscript supplements. Their READMEs document inputs, analysis units, data accessions, methods and limits. For convenience, the same five Python files contained in the supplements are also available under `analysis/`. The source-extraction scripts in supplement 4 accept original public proteomics search outputs. The three default checks begin at the supplied derived measurements. They do not reprocess raw sequencing reads or mass spectra.

## Data sources and licences

The project used GEO, TCGA, CPTAC/LinkedOmics, Human Protein Atlas, ProteomeXchange and author-deposited annotations. [Source and rights information](SOURCES.md) describes where to obtain the inputs. Source identifiers and original terms continue to apply to third-party values and annotations. Original analysis code is released under [MIT](LICENSE-CODE.txt). Authors' original results and documentation are under [CC BY 4.0](LICENSE-DATA.txt). The package contains numerical derivatives and provenance. Full-resolution H&E images remain with the original providers.
