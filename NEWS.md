
# splicelogic 1.2.0

* New event types: alternative transcription start and end sites
  (`find_atss()`, `find_ates()`) and spliced introns (`find_si()`).
* New input option: `prepare_exons_by_partition()` builds the input from two sets of
  transcripts, without a DTU analysis.
* `get_seq()` extracts exon sequences against a genome.
* Plotting functions for event results: `plot_transcripts()`,
  `plot_transcripts_by_gene()`, `plot_transcripts_by_pair()`,
  `plot_transcript_lfc()` and `event_palette()`.
* Long-form aliases for every finder, e.g. `find_skipped_exons()`.
* `find_se()`, `find_ie()`, `find_mxe()` and `find_all_events()` pass
  extra arguments on to `GenomicRanges::findOverlaps()`.
* Results carry `event_type`, `event_tx_id`, `event_estimate` and
  `event_<col>` for any `additional_columns`.
* The default `type` is now `"over"`; `"boundary"` matches only the
  splice site facing the candidate; and `"in"` is renamed `"identical"`.

# splicelogic 1.1.5

* Added plotting functions for event results: `plot_transcripts()`,
  `plot_transcripts_by_gene()`, `plot_transcripts_by_pair()`,
  `plot_transcript_lfc()` and `event_palette()` (require `ggplot2`).
* Added `find_si()` for spliced introns, the reverse of `find_ri()`.
* The default `type` of `find_se()`, `find_ie()`, `find_mxe()` and
  `find_all_events()` is now `"over"` (was `"boundary"`).
* `type = "boundary"` now matches only the splice site of each flanking
  exon that faces the candidate exon, on either strand.
* The `"in"` match type is renamed `"identical"`.
* `find_se()`, `find_ie()`, `find_mxe()` and `find_all_events()` pass
  extra arguments on to `GenomicRanges::findOverlaps()`, e.g.
  `minoverlap = 5L`.
* `plyranges` moved to Depends, so it is attached with `splicelogic`.

# splicelogic 1.1.4

* Added `prepare_exons_by_partition()` for building input from two sets
  of transcripts, without needing a estimate column from a DTU analysis.
* Added `find_atss()` and `find_ates()`
  for detecting alternative transcription start and end sites. Both
  are included in `find_all_events()`.
* Added `find_si()` (alias `find_spliced_introns()`), the reverse of
  `find_ri()`: an intron that an up transcript splices out but a down
  transcript retains. `find_ri()` gains an `inverse` argument, as
  `find_se()` has for `find_ie()`. Included in `find_all_events()`.
* `find_a5ss()` and `find_a3ss()` now label events per moved boundary,
  so an exon that moves both of its boundaries is reported as both an
  a5ss and an a3ss.
* The default `type` of `find_se()`, `find_ie()`, `find_mxe()` and
  `find_all_events()` is now `"over"`. `"boundary"` now requires the
  flank to share its inner boundary — the splice site facing the
  candidate — rather than either its start or its end, and is
  strand-aware. The `"in"` match type is renamed `"identical"`. Extra
  arguments (`...`) are passed on to `GenomicRanges::findOverlaps()`
  when matching flanking exons.
* Added `get_seq()` to extract exon sequences from a `GRanges` against
  a genome.
* Results now carry `event_type`,
  `event_tx_id` and `event_estimate`, and any columns registered via
  `additional_columns` are propagated to the partner transcript with an
  `event_` prefix.
  * Added long-form aliases for every finder (`find_skipped_exons()`,
  `find_included_exons()`, `find_mutually_exclusive_exons()`,
  `find_retained_introns()`, `find_alternative_5_prime_splice_sites()`,
  `find_alternative_3_prime_splice_sites()`,
  `find_alternative_transcription_start_sites()`,
  `find_alternative_transcription_end_sites()`).

# splicelogic 0.99.0

Submit to Bioconductor 3.23 
