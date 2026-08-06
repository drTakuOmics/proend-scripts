# Codex task: finalize the PAN2 TEM and AF3 repository extension

Work only on branch `docs/pan2-tem-af3-update`.

## Goal

Finalize the new PAN2 structural follow-up documentation using the files that the author uploads under `docs/pan2-extension/`.

The central evidence chain is:

1. ProEnd previously identified PAN2/MJ1494 as an HbYX-containing archaeal proteasome regulator and reported biochemical interaction/activity.
2. Negative-stain TEM followed by CryoSPARC analysis yielded a small number of particle classes consistent with PAN2-T20S assemblies.
3. AlphaFold 3 produced a PAN2 hexamer-T20S alpha-ring model containing 5 ATP, 1 ADP, and 6 Mg2+, with a predicted C-terminal interface near the alpha-ring pockets.

## Files to edit

- `docs/pan2-extension/README.md`
- `README.Rmd`
- `README.md`

Do not modify the analysis scripts.

## Expected author-supplied assets

Look for these files, or the closest clearly named equivalents:

```text
docs/pan2-extension/figures/pan2_cryosparc_workflow.png
docs/pan2-extension/figures/pan2_selected_2d_classes.png
docs/pan2-extension/figures/pan2_t20s_af3_rotation.gif
docs/pan2-extension/figures/pan2_t20s_af3_interface.png
```

Optional supporting files may be present under:

```text
docs/pan2-extension/cryosparc/
docs/pan2-extension/af3/
```

## Required edits

1. Replace only the relevant `TODO` entries in `docs/pan2-extension/README.md` when the corresponding information is available in uploaded files or author notes.
2. Verify every image path and ensure each figure renders on GitHub.
3. Preserve the evidence distinction:
   - published biochemical result;
   - exploratory TEM/CryoSPARC observation;
   - AF3 structural prediction.
4. Do not describe the AF3 model as a solved structure or the predicted interface as experimentally confirmed.
5. Do not claim a 3D reconstruction from the TEM dataset.
6. Add a compact preview section to both `README.Rmd` and `README.md`, preferably before the general conclusion. Keep it to one paragraph, one preview image or GIF, and one link to `docs/pan2-extension/README.md`.
7. Keep `README.Rmd` and `README.md` scientifically synchronized.

## Fixed workflow values

Use these values unless the author supplies a correction:

- 8 micrographs
- 80 kV
- 29,070x
- image dimensions: 5056 x 5056 px
- 21,534 blob-picked candidate particles
- 2,692 candidate particles per micrograph
- particle diameter range: 100-150, unit still requiring author confirmation
- 16,616 extracted particles
- 2,077 extracted particles per micrograph
- extraction box: 440 px
- 13,605 particles used for 2D classification
- 100 classes
- 4 selected classes
- 42 selected particles

Do not silently reinterpret the ambiguous `1.00 e/Å` metadata. Keep it exactly as supplied or mark it for confirmation.

## Fixed AF3 composition

- PAN2: 6 chains
- T20S alpha subunit: 7 chains
- ATP: 5
- ADP: 1
- Mg2+: 6

Call this a `PAN2 hexamer-T20S alpha-ring model`, not a full PAN2-20S atomic structure.

## Style

The dedicated extension page should read as a scrolling repository case study, not as a publication figure legend dump. Use short sections, one visual per narrative step, and concise captions.

## Validation

Before finishing:

- confirm all linked files exist;
- ensure no unsupported scientific claims were introduced;
- confirm `README.Rmd` and `README.md` agree;
- list unresolved metadata as explicit TODOs;
- provide a concise summary of changes and remaining author decisions in the commit or PR description.
