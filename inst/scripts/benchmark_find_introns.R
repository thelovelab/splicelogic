library(GenomicRanges)
library(microbenchmark)
library(bench)
library(ggplot2)
library(patchwork)

basic_introns_lapply <- function(gr) {
  grl <- split(gr, gr$tx_id)
  lapply(grl, gaps) 
}

# function to find introns given a GRanges object of exons using lapply and gaps
# keeps metadata columns gene_id and tx_id
introns_lapply <- function(gr) {
  # make GRangesList of exons split by tx_id 
  grl <- split(gr, gr$tx_id)
  # for each tx_id, find the introns using gaps() and keep metadata columns gene_id and tx_id
  intr_list <- lapply(names(grl), function(tx_id) {
    gr_tx <- grl[[tx_id]]
    intr <- gaps(gr_tx)
    mcols(intr)$gene_id <- unique(gr_tx$gene_id)
    mcols(intr)$tx_id   <- tx_id
    intr
  }) #output is a list of GRanges objects, one per tx_id
  # make a GRangesList of introns and unlist it to get a GRanges object
  GRangesList(intr_list) |> unlist(use.names = FALSE)
}

find_introns_v0 <- function(gr) {
  introns <- gr |>
    plyranges::group_by(tx_id) |>
     # introns are between exons, so we can use the start of the next exon and the end of the current exon to define them
    plyranges::mutate(
      intron_start = end + 1,
      intron_end = dplyr::lead(start) - 1,
    ) |>
    plyranges::filter(!is.na(intron_start) & !is.na(intron_end)) |>
    plyranges::mutate(
      start = intron_start,
      end = intron_end,
      intron = TRUE
    ) |>
    plyranges::ungroup() |>
    plyranges::select(-intron_start, -intron_end, -exon_rank)
  return(introns)
}

benchmark_and_plot <- function(sizes, x_var) {
  my_palette <- c(
  "#1f77b4",  # (31,119,180)
  "#d62728",  # (214,39,40)
  "#17becf",   # (23,190,207)
  "#ffbf00",  # (255,191,0)
  "#b5e188"  # (181,225,136)
)

  scale_results <- lapply(names(sizes), function(k) {
    s <- sizes[[k]]
    grk <- create_mock_data(s$genes, s$tx_per_gene, s$exons)

    mbk <- bench::mark(
      basic_introns_lapply  = basic_introns_lapply(grk),
      introns_lapply = introns_lapply(grk),
      find_introns_v0 = find_introns_v0(grk),
      find_introns = find_introns(grk),
      iterations = 1L,
      check = FALSE
    )
    sm <- summary(mbk, unit = "ms")  # optional: force a unit
    sm[[x_var]] <- s[[x_var]]
    sm
  })

  df <- do.call(rbind, scale_results) %>%
    as_tibble() %>%
    mutate(
      expression = as.character(expression),
      median_ms = as.numeric(median) * 1000, # convert bench_time → numeric ms
      mem_kb = as.numeric(mem_alloc) / 1024 # KB
    )

  p_time <- ggplot(df,  aes(x = .data[[x_var]], 
                    y = median_ms, 
                    color = expression, 
                    group = expression))  +
    geom_point(size = 3, alpha = 0.9) +
    geom_line(linewidth = 1, alpha = 0.7) +
    theme_classic(base_size = 14) +
    labs(
      x = x_var,
      y = "Median time (ms)",
      color = "Function",
      title = paste("Benchmark intron-finding functions"),
      subtitle = paste("Median run time across increasing number of ", x_var)
    )+
    scale_color_manual(values = my_palette)

  # memory plot
  p_mem <- ggplot(df, aes(x = .data[[x_var]], 
                    y = mem_kb, 
                    color = expression, 
                    group = expression))  +
    geom_point(size = 3, alpha = 0.9) +
    geom_line(linewidth = 1, alpha = 0.7) +
    theme_classic(base_size = 14) +
    labs(
      x = x_var,
      y = "Memory allocated (KB)",
      color = "Function",
      subtitle = paste("Memory allocated across increasing number of ", x_var)
    )+
    scale_color_manual(values = my_palette)

  plot <- wrap_plots(list(p_time, p_mem), ncol = 2)
  # save the plot to a file
  out_dir <- file.path("inst", "scripts", "plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  filename <- file.path(out_dir, sprintf("plots_benchmark_find_introns_by_%s.png", x_var))

  ggsave(filename = filename, plot = plot, width = 16, height = 6, bg = "white")

}


sizes_genes <- list(
  minimal = list(genes=1,  tx_per_gene=1, exons=2),
  verysmall  = list(genes=3,  tx_per_gene=5, exons=5),        
  small  = list(genes=10,  tx_per_gene=5, exons=5),
  medium = list(genes=50,  tx_per_gene=5, exons=5),
  large = list(genes=100, tx_per_gene=5, exons=5)
)
sizes_exons <- list(
  minimal = list(genes=1,  tx_per_gene=1, exons=2),
  large = list(genes=10, tx_per_gene=20, exons=100),
  medium = list(genes=10,  tx_per_gene=20, exons=50),
  small  = list(genes=10,  tx_per_gene=20, exons=10),
  verysmall  = list(genes=10,  tx_per_gene=20, exons=5)
)

benchmark_and_plot(sizes_genes, "genes")
benchmark_and_plot(sizes_exons, "exons")

