# Package index

## Package overview

- [`splicelogic`](https://thelovelab.github.io/splicelogic/reference/splicelogic-package.md)
  [`splicelogic-package`](https://thelovelab.github.io/splicelogic/reference/splicelogic-package.md)
  : splicelogic: differential transcripts to splice events

## Input preprocessing

- [`preprocess_input()`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md)
  : Preprocess input GRanges object for splicing event calculation
- [`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
  : Prepare exon ranges from a TxDb and DTU results table

## Event detection

- [`calc_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
  [`calc_included_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
  [`calc_mx_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
  [`calc_retained_introns()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
  [`calc_alt_ss()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
  [`calc_a5ss()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
  [`calc_a3ss()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
  [`calc_all_events()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
  : Calculate splice events from a GRanges object

## Mock data generation

- [`create_mock_data()`](https://thelovelab.github.io/splicelogic/reference/create_mock_data.md)
  : Create mock GRanges data for splicing event testing
- [`generate_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_mx()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_retained_introns()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_a5ss()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_a3ss()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  : Generate mock splice events in a GRanges object
