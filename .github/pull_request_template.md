## Change and scientific contract

- What changes, and why?
- Which inputs, units, coordinate frames, replicate populations or saved-output paths are affected?
- Which existing outputs/checkpoints must remain compatible?

## Validation

- [ ] Targeted tests and affected CLI path exercised (list commands/results).
- [ ] Python 3.11 / 3.13 checks reviewed.
- [ ] Real-trajectory checks reviewed when relevant; limitations stated.

## Plotting changes (if applicable)

- [ ] Opened **MDscope figure review** Actions artifact: `mdscope-figure-review/index.html`.
- [ ] Reviewed base/head plots for axes, labels, units, legends, missing data, clipping and color mappings.
- [ ] Included a link to the visual artifact and described the intended differences (or explained why there should be none).
- [ ] Added a deterministic saved-data example if a changed renderer was not covered by the gallery.

Synthetic previews verify rendering only; they do not demonstrate statistical validity or physical convergence.
