# Changelog

## splicelogic 0.0.67

- Initial release as a Bioconductor package.
- Translate differential transcript usage results into discrete splice
  events.
- Event detection functions:
  [`calc_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md),
  [`calc_included_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md),
  [`calc_mx_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md),
  [`calc_retained_introns()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md),
  [`calc_a5ss()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md),
  [`calc_a3ss()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md),
  and
  [`calc_all_events()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md).
- Input preprocessing with
  [`preprocess_input()`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md)
  to standardize exon metadata.
- Input preparation with
  [`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
  to combine transcript annotations and DTU results.
