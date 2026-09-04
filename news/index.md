# Changelog

## splicelogic 1.1.4

- Added
  [`prepare_exons_by_partition()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons_by_partition.md)
  for building input from two sets of transcripts, without needing a
  estimate column from a DTU analysis.
- Added
  [`find_atss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  and
  [`find_ates()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  for detecting alternative transcription start and end sites. Both are
  included in
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md).
- [`find_a5ss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  and
  [`find_a3ss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  now label events per moved boundary, so an exon that moves both of its
  boundaries is reported as both an a5ss and an a3ss.
- Added
  [`get_seq()`](https://thelovelab.github.io/splicelogic/reference/get_seq.md)
  to extract exon sequences from a `GRanges` against a genome.
- Results now carry `event_type`, `event_tx_id` and `event_estimate`,
  and any columns registered via `additional_columns` are propagated to
  the partner transcript with an `event_` prefix.
  - Added long-form aliases for every finder
    ([`find_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
    [`find_included_exons()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
    [`find_mutually_exclusive_exons()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
    [`find_retained_introns()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
    [`find_alternative_5_prime_splice_sites()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
    [`find_alternative_3_prime_splice_sites()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
    [`find_alternative_transcription_start_sites()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
    [`find_alternative_transcription_end_sites()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)).

## splicelogic 0.99.0

Submit to Bioconductor 3.23
