# Changelog

## splicelogic 1.2.0

- New event types: alternative transcription start and end sites
  ([`find_atss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_ates()`](https://thelovelab.github.io/splicelogic/reference/find_events.md))
  and spliced introns
  ([`find_si()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)).
- New input option:
  [`prepare_exons_by_partition()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons_by_partition.md)
  builds the input from two user-defined sets of transcripts (up and
  down). Unlike
  [`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md),
  it does not need DTU results.
- [`get_seq()`](https://thelovelab.github.io/splicelogic/reference/get_seq.md)
  extracts exon sequences against a genome.
- Plotting functions for event results:
  [`plot_transcripts()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts.md),
  [`plot_transcripts_by_gene()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts_by_gene.md),
  [`plot_transcripts_by_pair()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts_by_pair.md),
  [`plot_transcript_lfc()`](https://thelovelab.github.io/splicelogic/reference/plot_transcript_lfc.md)
  and
  [`event_palette()`](https://thelovelab.github.io/splicelogic/reference/event_palette.md).
- Long-form aliases for every finder,
  e.g. [`find_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/find_events.md).
- [`find_se()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_ie()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_mxe()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  and
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  pass extra arguments on to
  [`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).
- Results carry `event_type`, `event_tx_id`, `event_estimate` and
  `event_<col>` for any `additional_columns`.
- The default `type` is now `"over"`; `"boundary"` matches only the
  splice site facing the candidate; and `"in"` is renamed `"identical"`.

## splicelogic 1.1.5

- Added plotting functions for event results:
  [`plot_transcripts()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts.md),
  [`plot_transcripts_by_gene()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts_by_gene.md),
  [`plot_transcripts_by_pair()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts_by_pair.md),
  [`plot_transcript_lfc()`](https://thelovelab.github.io/splicelogic/reference/plot_transcript_lfc.md)
  and
  [`event_palette()`](https://thelovelab.github.io/splicelogic/reference/event_palette.md)
  (require `ggplot2`).
- Added
  [`find_si()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  for spliced introns, the reverse of
  [`find_ri()`](https://thelovelab.github.io/splicelogic/reference/find_events.md).
- The default `type` of
  [`find_se()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_ie()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_mxe()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  and
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  is now `"over"` (was `"boundary"`).
- `type = "boundary"` now matches only the splice site of each flanking
  exon that faces the candidate exon, on either strand.
- The `"in"` match type is renamed `"identical"`.
- [`find_se()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_ie()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_mxe()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  and
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  pass extra arguments on to
  [`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html),
  e.g. `minoverlap = 5L`.
- `plyranges` moved to Depends, so it is attached with `splicelogic`.

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
- Added
  [`find_si()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  (alias
  [`find_spliced_introns()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)),
  the reverse of
  [`find_ri()`](https://thelovelab.github.io/splicelogic/reference/find_events.md):
  an intron that an up transcript splices out but a down transcript
  retains.
  [`find_ri()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  gains an `inverse` argument, as
  [`find_se()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  has for
  [`find_ie()`](https://thelovelab.github.io/splicelogic/reference/find_events.md).
  Included in
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md).
- [`find_a5ss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  and
  [`find_a3ss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  now label events per moved boundary, so an exon that moves both of its
  boundaries is reported as both an a5ss and an a3ss.
- The default `type` of
  [`find_se()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_ie()`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
  [`find_mxe()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  and
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  is now `"over"`. `"boundary"` now requires the flank to share its
  inner boundary — the splice site facing the candidate — rather than
  either its start or its end, and is strand-aware. The `"in"` match
  type is renamed `"identical"`. Extra arguments (`...`) are passed on
  to
  [`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  when matching flanking exons.
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
