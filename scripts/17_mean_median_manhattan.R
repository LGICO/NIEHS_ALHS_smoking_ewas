library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

project_dir <- "output"

manifest_dir <- "references"

epic_cpgs <- data.table::fread(
  file.path(manifest_dir, "EPIC.hg38.manifest.gencode.v36.tsv.gz"),
  select="probeID"
)$probeID

epicv2_cpgs <- data.table::fread(
  file.path(manifest_dir, "EPICv2.hg38.manifest.gencode.v41.tsv.gz"),
  select="probeID"
)$probeID

hm450_cpgs <- data.table::fread(
  file.path(manifest_dir, "HM450.hg38.manifest.tsv.gz"),
  select="Probe_ID"
)$Probe_ID

previous_array_cpgs <- unique(sub(
  "_.*$", "", c(epic_cpgs, epicv2_cpgs, hm450_cpgs)
))

plot_stacked_manhattan <- function(current_file, former_file, output_file) {
  current <- data.table::fread(current_file) %>%
    dplyr::transmute(
      CpG=CpG_ID,
      chr=as.integer(sub("^chr", "", CpG_chrm)),
      pos=as.numeric(CpG_beg),
      pvalue=as.numeric(P),
      gene=genesUniq
    ) %>%
    dplyr::filter(chr %in% 1:22, !is.na(pos), !is.na(pvalue))
  
  former <- data.table::fread(former_file) %>%
    dplyr::transmute(
      CpG=CpG_ID,
      chr=as.integer(sub("^chr", "", CpG_chrm)),
      pos=as.numeric(CpG_beg),
      pvalue=as.numeric(P),
      gene=genesUniq
    ) %>%
    dplyr::filter(chr %in% 1:22, !is.na(pos), !is.na(pvalue))
  
  bonf <- 0.05 / nrow(current)
  
  chromosome_offsets <- dplyr::bind_rows(current, former) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(chr_length=max(pos), .groups="drop") %>%
    dplyr::arrange(chr) %>%
    dplyr::mutate(offset=cumsum(chr_length)-chr_length)
  
  plot_df <- dplyr::bind_rows(
    current=current,
    former=former,
    .id="status"
  ) %>%
    dplyr::left_join(chromosome_offsets, by="chr") %>%
    dplyr::mutate(
      rel_pos=pos+offset,
      logged_p=-log10(pmax(pvalue, .Machine$double.xmin)),
      status=factor(status, levels=c("current", "former"))
    ) %>%
    dplyr::arrange(chr, pos)
  
  axis_df <- plot_df %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(chr_center=(min(rel_pos)+max(rel_pos))/2, .groups="drop")
  
  select_labels <- function(data, group_name) {
    data %>%
      dplyr::filter(
        status == group_name,
        pvalue <= bonf,
        !CpG %in% previous_array_cpgs
      ) %>%
      dplyr::group_by(chr) %>%
      dplyr::slice_max(logged_p, n=1, with_ties=FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::slice_max(logged_p, n=12, with_ties=FALSE) %>%
      dplyr::mutate(
        gene=dplyr::if_else(is.na(gene), "", gene),
        label=dplyr::if_else(gene == "", CpG, paste0(CpG, "\n", gene))
      )
  }
  
  current_labels <- select_labels(plot_df, "current")
  former_labels <- select_labels(plot_df, "former")
  
  current_max <- ceiling(max(plot_df$logged_p[plot_df$status == "current"], na.rm=TRUE))
  former_max <- ceiling(max(plot_df$logged_p[plot_df$status == "former"], na.rm=TRUE))
  
  current_plot <- ggplot2::ggplot(
    plot_df %>% dplyr::filter(status == "current"),
    ggplot2::aes(x=rel_pos, y=logged_p, color=factor(chr))
  ) +
    ggplot2::geom_point(size=1) +
    ggplot2::scale_color_manual(values=rep(c("#3261cf", "#88a1db"), 11)) +
    ggplot2::scale_x_continuous(
      labels=axis_df$chr,
      breaks=axis_df$chr_center,
      expand=ggplot2::expansion(mult=0.01)
    ) +
    ggplot2::scale_y_continuous(
      limits=c(0, current_max+5),
      expand=ggplot2::expansion(mult=c(0.02, 0))
    ) +
    ggplot2::geom_hline(
      yintercept=-log10(1e-5), color="grey", linetype="dashed", linewidth=0.5
    ) +
    ggplot2::geom_hline(
      yintercept=-log10(bonf), color="#150c6b", linewidth=0.5
    ) +
    ggrepel::geom_text_repel(
      data=current_labels,
      ggplot2::aes(label=label),
      color='black',
      size=5,
      segment.size=0.2,
      point.padding=0.3,
      box.padding=1.5,
      force=4,
      force_pull=0.1,
      max.overlaps=Inf,
      show.legend=FALSE
    ) +
    ggplot2::annotate(
      "text", x=Inf, y=current_max+4,
      label=sprintf("Bonferroni cutoff: P = %.3e", bonf),
      hjust=1.05, vjust=1, size=4, color="#150c6b"
    ) +
    ggplot2::labs(
      x=NULL,
      y=expression(atop("Current smoker vs. never smoker", -log[10](P)))
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position="none",
      axis.title.x=ggplot2::element_blank(),
      plot.margin=ggplot2::margin(b=0, l=10),
      axis.text=ggplot2::element_text(size=16),
      axis.title=ggplot2::element_text(size=24, face="bold")
    )
  
  former_plot <- ggplot2::ggplot(
    plot_df %>% dplyr::filter(status == "former"),
    ggplot2::aes(x=rel_pos, y=logged_p, color=factor(chr))
  ) +
    ggplot2::geom_point(size=1) +
    ggplot2::scale_color_manual(values=rep(c("#a33636", "#d99898"), 11)) +
    ggplot2::scale_x_continuous(
      labels=axis_df$chr,
      breaks=axis_df$chr_center,
      expand=ggplot2::expansion(mult=0.01)
    ) +
    ggplot2::scale_y_continuous(
      limits=c(0, former_max+5),
      expand=ggplot2::expansion(mult=c(0.02, 0))
    ) +
    ggplot2::geom_hline(
      yintercept=-log10(1e-5), color="grey", linetype="dashed", linewidth=0.5
    ) +
    ggplot2::geom_hline(
      yintercept=-log10(bonf), color="#150c6b", linewidth=0.5
    ) +
    ggrepel::geom_text_repel(
      data=former_labels,
      ggplot2::aes(label=label),
      color="black",
      size=5,
      segment.size=0.2,
      point.padding=0.3,
      box.padding=1.5,
      force=4,
      force_pull=0.1,
      max.overlaps=Inf,
      show.legend=FALSE
    ) +
    ggplot2::labs(
      x="Chromosome",
      y=expression(atop("Former smoker vs. never smoker", -log[10](P)))
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position="none",
      plot.margin=ggplot2::margin(t=0, l=10),
      axis.text=ggplot2::element_text(size=16),
      axis.title=ggplot2::element_text(size=24, face="bold")
    )
  
  grDevices::pdf(output_file, height=18, width=20)
  gridExtra::grid.arrange(current_plot, former_plot, nrow=2)
  grDevices::dev.off()
}

run_collapsed_plots <- function(collapse_method) {
  output_dir <- file.path(project_dir, paste0("cpg_collapsed_", collapse_method))
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  
  plot_stacked_manhattan(
    current_file=file.path(output_dir, "weighted_euro_current_vs_never.tsv"),
    former_file=file.path(output_dir, "weighted_euro_former_vs_never.tsv"),
    output_file=file.path(
      output_dir,
      paste0("weighted_euro_stacked_manhattan_", collapse_method, ".pdf")
    )
  )
  
  plot_stacked_manhattan(
    current_file=file.path(output_dir, "weighted_trans_current_vs_never.tsv"),
    former_file=file.path(output_dir, "weighted_trans_former_vs_never.tsv"),
    output_file=file.path(
      output_dir,
      paste0("weighted_trans_stacked_manhattan_", collapse_method, ".pdf")
    )
  )
}

run_collapsed_plots("mean")
run_collapsed_plots("median")

