# splicelogic 0.0.67

* Initial release as a Bioconductor package.
* Translate differential transcript usage results into
  discrete splice events.
* Event detection functions: `calc_skipped_exons()`,
  `calc_included_exons()`, `calc_mx_exons()`,
  `calc_retained_introns()`, `calc_a5ss()`, `calc_a3ss()`,
  and `calc_all_events()`.
* Input preprocessing with `preprocess_input()` to
  standardize exon metadata.
* Input preparation with `prepare_exons()` to combine transcript
  annotations and DTU results.