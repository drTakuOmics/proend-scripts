# PAN2 structural follow-up: negative-stain TEM and AlphaFold 3 model

> **Status:** Repository extension in progress. This page summarizes exploratory structural follow-up to the published ProEnd PAN2 results. The TEM observations and AF3 model should not be interpreted as an experimentally resolved atomic structure.

## Background

ProEnd identified the archaeal PAN2 protein (MJ1494; [UniProt Q58889](https://www.uniprot.org/uniprotkb/Q58889/entry)) as an HbYX-containing proteasome regulator. The published ProEnd study reported PAN2 interaction with and activation of the archaeal 20S proteasome: [Salcedo et al., BMC Genomics (2024)](https://link.springer.com/article/10.1186/s12864-024-10864-4).

This extension adds two complementary structural observations:

1. negative-stain TEM followed by CryoSPARC particle analysis; and
2. an AlphaFold 3 model of a PAN2 hexamer positioned on the T20S alpha ring.

The aim is to document the follow-up workflow and provide a structural hypothesis consistent with the previously published biochemical results.

---

## 1. Negative-stain TEM dataset

A purified preparation containing PAN2 and T20S was examined by negative-stain transmission electron microscopy. **TODO:** confirm the exact sample composition, nucleotide condition, stain, incubation conditions, and whether the preparation was purified as a preassembled complex or mixed before grid preparation.

Eight micrographs were used for this exploratory CryoSPARC workflow.

![CryoSPARC workflow for the PAN2-T20S negative-stain TEM dataset](figures/pan2_cryosparc_workflow.png)

**Figure 1. CryoSPARC workflow used for exploratory PAN2-T20S particle analysis.** Eight micrographs were imported, particles were detected by blob picking, extracted, classified into 100 two-dimensional classes, and manually curated. Four classes containing 42 particles were retained as candidate PAN2-T20S views for visualization. The low retained particle number means that these classes are treated as qualitative observations rather than a high-resolution structural reconstruction.

### Workflow summary

| Step | CryoSPARC output | Reported value |
|---|---|---:|
| Import micrographs | Micrographs | 8 |
| Acquisition | Accelerating voltage | 80 kV |
| Acquisition | Exposure metadata | 1.00 e/Å |
| Acquisition | Magnification | 29,070x |
| Acquisition | Image dimensions | 5056 x 5056 px |
| Blob picker | Candidate particles | 21,534 |
| Blob picker | Mean candidates per micrograph | 2,692 |
| Blob picker | Particle diameter range | 100-150 |
| Extraction | Extracted particles | 16,616 |
| Extraction | Mean extracted particles per micrograph | 2,077 |
| Extraction | Extraction box | 440 px |
| 2D classification | Particles classified | 13,605 |
| 2D classification | Number of classes | 100 |
| 2D selection | Retained classes | 4 |
| 2D selection | Retained particles | 42 |

**TODO:** confirm the units used for the blob-picking diameter and exposure metadata exactly as reported by the microscope/CryoSPARC project.

---

## 2. Representative particles and 2D classes

The particle population was heterogeneous. Most classes were excluded during manual curation, while a small subset showed shapes compatible with PAN2-like, T20S-like, or putative PAN2-T20S assemblies.

![Representative PAN2-T20S particle views or selected 2D classes](figures/pan2_selected_2d_classes.png)

**Figure 2. Representative selected particle views from the exploratory negative-stain dataset.** The selected classes contain particles with morphology consistent with a PAN2-T20S assembly in different orientations. Because only 42 particles were retained across four classes, these observations support the presence of candidate assemblies but do not establish stoichiometry, symmetry, resolution, or a unique three-dimensional conformation.

### Interpretation of the TEM result

The selected particles provide visual evidence compatible with PAN2 association with the T20S proteasome. The result is presented conservatively because the dataset contains only eight micrographs and a small number of retained candidate particles. No three-dimensional reconstruction is claimed in this repository update.

---

## 3. AlphaFold 3 PAN2 hexamer-T20S alpha-ring model

To explore a possible molecular arrangement, AlphaFold 3 was used to model the PAN2 hexamer together with the seven alpha subunits of the archaeal T20S proteasome.

### Model composition

| Component | Copy number |
|---|---:|
| PAN2 | 6 chains |
| T20S alpha subunit | 7 chains |
| ATP | 5 |
| ADP | 1 |
| Mg2+ | 6 |

The model therefore represents a **PAN2 hexamer-T20S alpha-ring complex**, not the complete four-ring 20S proteasome.

![Rotation of the AF3 PAN2 hexamer-T20S alpha-ring model](figures/pan2_t20s_af3_rotation.gif)

**Figure 3. Rotation of the static AlphaFold 3 PAN2 hexamer-T20S alpha-ring prediction.** The animation is intended to show the overall predicted arrangement and does not represent molecular dynamics, ATP hydrolysis, or an experimentally observed conformational trajectory.

### Predicted interface

![Predicted PAN2 C-terminal interface with the T20S alpha ring](figures/pan2_t20s_af3_interface.png)

**Figure 4. Predicted PAN2-T20S alpha-ring interface.** The AF3 model positions the PAN2 C-terminal region near the intersubunit pockets of the T20S alpha ring, providing a structural hypothesis for HbYX-dependent engagement. This predicted interface is compatible with proteasome activation but has not been validated at atomic resolution by the TEM dataset.

### AF3 confidence information

Add the values from the selected AF3 output before finalizing this section.

| Metric | Value |
|---|---:|
| AF3 model/sample identifier | TODO |
| Ranking score | TODO |
| Overall pTM | TODO |
| Overall ipTM | TODO |
| PAN2-alpha-ring pairwise ipTM | TODO |
| Interface PAE summary | TODO |
| Local pLDDT for PAN2 C termini | TODO |

A PAE heatmap or confidence panel may be added only if it remains readable in the scrolling documentation layout.

---

## 4. Integrated interpretation

The published biochemical data established PAN2 as a candidate archaeal 20S proteasome regulator. The negative-stain TEM dataset now provides a small set of particle views compatible with PAN2-T20S assemblies, while the AF3 prediction proposes a possible arrangement of the PAN2 hexamer on the T20S alpha ring and a candidate C-terminal interface.

Together, these results extend the ProEnd PAN2 case study from sequence-based discovery and biochemical activity toward a structural working model. The TEM data support the presence of candidate assemblies at the particle level, whereas the AF3 model supplies an atomic hypothesis that remains to be experimentally validated.

---

## 5. Reproducibility files

The following files should be deposited with this extension when approved for public release:

```text
docs/pan2-extension/
├── README.md
├── CODEX_TASK.md
├── figures/
│   ├── pan2_cryosparc_workflow.png
│   ├── pan2_selected_2d_classes.png
│   ├── pan2_t20s_af3_rotation.gif
│   └── pan2_t20s_af3_interface.png
├── cryosparc/
│   ├── workflow.json
│   └── processing_notes.md
└── af3/
    ├── input_summary.md
    ├── confidence_summary.csv
    └── interface_notes.md
```

Raw microscope data, unpublished source images, and AF3 output archives should only be added if their public release has been approved.

---

## Remaining author checks

- [ ] Confirm the exact TEM sample composition and nucleotide condition.
- [ ] Confirm stain, grid type, microscope model, and acquisition metadata.
- [ ] Confirm the units associated with `1.00 e/Å` and the picking diameter `100-150`.
- [ ] Upload the final workflow image.
- [ ] Upload the final selected-particle or 2D-class panel.
- [ ] Upload the AF3 rotation GIF and static interface image.
- [ ] Add AF3 pTM, ipTM, pairwise-interface, PAE, and pLDDT values that are available.
- [ ] Review the wording before linking this page from the main README.
