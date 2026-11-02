# Package index

## Package overview

- [`splicelogic`](https://thelovelab.github.io/splicelogic/reference/splicelogic-package.md)
  [`splicelogic-package`](https://thelovelab.github.io/splicelogic/reference/splicelogic-package.md)
  : splicelogic: differential transcripts to splice events

## Input preprocessing

- [`preprocess()`](https://thelovelab.github.io/splicelogic/reference/preprocess.md)
  : Preprocess input GRanges object for splicing event calculation
- [`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
  : Prepare exon ranges from a TxDb and DTU results table
- [`prepare_exons_by_partition()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons_by_partition.md)
  : Prepare exons from two transcript partitions

## Event detection

- [`find_se()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_ie()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_included_exons()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_mxe()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_mutually_exclusive_exons()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_ri()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_retained_introns()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_si()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_spliced_introns()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_a5ss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_alternative_5_prime_splice_sites()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_a3ss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_alternative_3_prime_splice_sites()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_atss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_alternative_transcription_start_sites()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_ates()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_alternative_transcription_end_sites()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  : Find splice events from annotated exons

## Plotting

- [`plot_transcript_lfc()`](https://thelovelab.github.io/splicelogic/reference/plot_transcript_lfc.md)
  : Plot a per-transcript fold change for the transcripts of splicelogic
  events
- [`plot_transcripts()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts.md)
  : Plot the transcripts involved in splicelogic events
- [`plot_transcripts_by_gene()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts_by_gene.md)
  : Plot the transcripts involved in splicelogic events, one plot per
  gene
- [`plot_transcripts_by_pair()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts_by_pair.md)
  : Plot the transcripts involved in splicelogic events, one plot per
  pair
- [`event_palette()`](https://thelovelab.github.io/splicelogic/reference/event_palette.md)
  : Default fill colours for splicelogic event types

## Sequence extraction

- [`get_seq()`](https://thelovelab.github.io/splicelogic/reference/get_seq.md)
  : Extract sequences for GRanges

## Mock data generation

- [`create_mock_data()`](https://thelovelab.github.io/splicelogic/reference/create_mock_data.md)
  : Create mock GRanges data for splicing event testing
- [`generate_se()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_mxe()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_ri()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_a5ss()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_a3ss()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_atss()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  [`generate_ates()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
  : Generate mock splice events in a GRanges object
