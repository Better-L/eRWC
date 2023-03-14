# function to check/load/install required packages, derived from http://tinyurl.com/jjwyzph
checkpackages <- function(x){
  for (i in x) {
    if (!require(i, character.only = TRUE)) install.packages(i, dependencies = TRUE); library(i, character.only = TRUE )
  }
}

checkBiocpackages <- function(x){
  for (i in x) {
    if (!require(i, character.only = TRUE)) BiocManager::install(i); library(i, character.only = TRUE )
  }
}

# a function derived from "R in action" book
mystats <- function(x, na.omit=TRUE){
  if (na.omit)
    x <- x[!is.na(x)]
  n <- round(length(x), digits = 1)
  mx <- round(max(x), digits = 1)
  mn <- round(min(x), digits = 1)
  m <- round(mean(x), digits = 1)
  s <- round(sd(x), digits = 1)
  return(c(n = n, mean = m, stdev = s, max = mx, min = mn))
}

test_two_grps <- function(i, df, sample_names){
  desc_str <- paste(i, sample_names[1], sample_names[2])
  message(strrep("-", 80))
  message(i, appendLF = TRUE)
  eq_var = FALSE
  data <- filter(df, variable == i)
  f_var <-  with(data, var.test(CFU_mL[sample == sample_names[1]], CFU_mL[sample == sample_names[2]]))
  f_var$data.name <- desc_str
  print(f_var)
  message(strrep(" ", 5), strrep("-",66))
  if (f_var$p.value >= 0.05) {
    eq_var = TRUE
  }
  t_var <-  with(data, t.test(CFU_mL[sample == sample_names[1]], CFU_mL[sample == sample_names[2]], var.equal = eq_var))
  t_var$data.name <- desc_str
  print(t_var)
}

test_more_grps <- function(i, df){
  message(strrep("-", 80))
  message(i, appendLF = TRUE)
  data <- filter(df, variable == i)
  b_var <-  with(data, bartlett.test(CFU_mL ~ sample))
  b_var$data.name <- i
  print(b_var)
  message(strrep(" ", 5), strrep("-",66))
  if (b_var$p.value >= 0.05) {
    # ANOVA
    aov_var <- with(data, aov(CFU_mL ~ sample))
    aov_var$data.name <- i
    print(summary(aov_var))
    pval <- unlist(summary(aov_var))['Pr(>F)1']
    if (pval <= 0.05) {
      message(strrep(" ", 5), strrep("-",66))
      print(TukeyHSD(aov_var, "sample"))
      DFerr <- df.residual(aov_var)
      MSerr <- deviance(aov_var)/DFerr
      plot(HSD.test(y = aov_var, trt = "sample", DFerror = DFerr, MSerror = MSerr, alpha = 0.05, main = i, console = FALSE), main = i, las = 2)
    }
  } 
  else {
    kru_var <- with(data, kruskal.test(CFU_mL ~ sample))
    kru_var$data.name <- i
    print(kru_var)
    if (kru_var$p.value <= 0.05) {
      message(strrep(" ", 5), strrep("-",66))
      pht <- dunnTest(CFU_mL ~ sample, data = data, method = "bonferroni")
      pht$data.name <- paste(i, "(CFU_mL ~ sample)")
      print(pht, dunn.test.results = TRUE)
    }
  }
}

load_ngs_data <- function(fp) {
  tax_table <- read.xlsx(fp, sheet = "tax_table", rowNames = TRUE)
  tax_table$ASV <- row.names(tax_table) # we keep the ASV vector
  row.names(tax_table) <- tax_table$Sequence # necessary to construct the phylogenetic tree afterwards
  
  ASV_table <- read.xlsx(fp, sheet = "ASV_table", rowNames = TRUE)
  row.names(ASV_table) <- tax_table$Sequence # necessary to construct the phylogenetic tree afterwards
  
  sample_data <- read.xlsx(fp, sheet = "sample_data", rowNames = TRUE)
  
  # Creating the phyloseq object---------------------------------------------------------------------------
  # Follow https://joey711.github.io/phyloseq/index.html \Tutorials
  # We have ASV, but we will refer to OTU to align with Phyloseq dictionary
  # OTU and taxa tables must be matrix
  otumat <- as.matrix(ASV_table, stringsAsFactors = FALSE)
  taxmat <- as.matrix(tax_table, stringsAsFactors = FALSE)
  
  OTU = otu_table(otumat, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  SAM = sample_data(sample_data)
  
  physeq = phyloseq(OTU, TAX, SAM)
  
  return(physeq)
}

import_and_select_PCA_data <- function(fp_piv, fp_meta, redundant_f=FALSE, rm_samples=FALSE, select_samples=FALSE) {
  meta <- read.csv(fp_meta)
  data <- read.csv(fp_piv)
  if (!isFALSE(select_samples)) {
    s_filter <- meta[[select_samples[[1]]]] == select_samples[[2]]
    meta <- dplyr::filter(meta, s_filter)
    data <- dplyr::filter(data, s_filter)
  }
  if (!isFALSE(rm_samples)) {
    s_filter <- !(meta[, 1] %in% rm_samples)
    meta <- dplyr::filter(meta, s_filter)
    data <- dplyr::filter(data, s_filter)
  }
  if (!isFALSE(redundant_f)) {
    s_filter <- !(names(data) %in% redundant_f)
    data <- dplyr::select(data, names(data)[s_filter])
  }
  
  # Removal of constant features(i.e. no differences among samples)
  data <- data[vapply(data, function(x) length(unique(x)) > 1, logical(1L))]
  all_data <- merge(meta, data, by = names(meta)[1])
  return(all_data)
}

get_pivot_data <- function(dataset, fp_meta) {
  meta_len <- length(names(read.csv(fp_meta))) + 1
  data <- dataset[,c(1, meta_len:length(dataset))]
  data <- data %>% remove_rownames %>% column_to_rownames(var = "sampleID")
  return(data)
}

draw_PCA <- function(data, fp_meta, title = "PCA", color_col="AdjunctCulture", xaxis_tick_label_size = 12, cpal = FALSE) {
  if (isFALSE(cpal)) {
    # colorblind palette
    cpal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
              "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  }
  
  PCA_plt_data <- PCA(get_pivot_data(data, fp_meta), graph = FALSE)
  PCA_plt <- fviz_pca_ind(PCA_plt_data, 
                          repel = TRUE, 
                          label = "none", 
                          geom = "point",
                          palette = cpal,
                          title = title
  ) +
    geom_point(size = 5, aes(colour = data[, color_col])) +
    theme(text = element_text(size = xaxis_tick_label_size),
          plot.title = element_text(size = xaxis_tick_label_size + 2, hjust = 0.5),
          legend.text = element_text(size = xaxis_tick_label_size),
          legend.title = element_blank())
  return(PCA_plt)
}

draw_PCA_contrib <- function(data, fp_meta, title, xaxis_tick_label_size = 12, cpal = FALSE) {
  if (isFALSE(cpal)) {
    # colorblind palette
    cpal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
              "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  }
  PCA_plt_data <- PCA(get_pivot_data(data, fp_meta), graph = FALSE)
  PCA_plt <- fviz_pca_var(PCA_plt_data, 
                          col.var = "contrib",
                          gradient.cols = cpal,
                          repel = TRUE, # Avoid text overlapping
                          title = title
  ) + 
    theme(text = element_text(size = xaxis_tick_label_size),
          plot.title = element_text(size = xaxis_tick_label_size + 2, hjust = 0.5),
          legend.text = element_text(size = xaxis_tick_label_size),
          legend.title = element_blank())
  return(PCA_plt)
}

supplementary_plots <- function(data, fp_meta, title, color_columns=c("SubExperiment", "culture_type")) {
  # colorblind palette
  cpal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  data <- data %>%
    mutate(SubExperiment = paste("Day", substr(SubExperiment, nchar(SubExperiment), nchar(SubExperiment))))
  culture_dict <- list("None" = "control", "# remove e in eRWCeRM.y " = "y", "eRM.H.y" = "y", "eRM.HS.y" = "y", "eRM.o" = "o", "eRM.H.o" = "o", "eRM.HS.o" = "o")
  data <- data %>%
    mutate(AdjunctCulture = recode_factor(AdjunctCulture, !!!culture_dict))
  p1 <- draw_PCA(data, fp_meta, title = title, color_col = color_columns[1], cpal = cpal) #c('#b2182b','#fee391','#67a9cf','#2166ac'))
  p2 <- draw_PCA(data, fp_meta, title = title, color_col = color_columns[2], cpal = cpal) #c('#b2182b','#fee391','#67a9cf','#2166ac'))
  p3 <- draw_PCA_contrib(data, fp_meta, title = paste("Variables PCA - ", title), cpal = c("#00AFBB", "#E7B800", "#FC4E07"))
  plt <- ((p1 + p2 ) / (p3 + guide_area()) ) +  plot_layout(guides = 'auto') + 
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 18, face = "bold"))
  return(plt)
}

draw_ha_stacked_plot <- function(df, sampleid_vector, 
                                 abund_thresh = 1.0, sp_res_th = 0.6,
                                 x_labels, yaxis_tick_label_size = 10,
                                 yaxis_label_size = 10, xaxis_tick_label_size = 12,
                                 Xaxis_rot = 0, ordered_sp_labels,
                                 df_sorter, species_colors, wrap_group, x_factor,
                                 ordered_species) {
  
  plt_title <- paste(c("≥", abund_thresh, "% average abundance"), collapse = " ")
  ha_Plot <- ggplot(arrange(df, !!!df_sorter) %>%
                      mutate(Species = factor(Species, levels = rev(ordered_species))),
                    aes(x = factor(df[[x_factor]], levels = sampleid_vector),
                        y = Abundance, fill = Species)) + 
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = c("black", species_colors[1:length(unique(df$Species)) - 1]), 
                      labels = rev(ordered_sp_labels)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10.0), expand = expansion(mult = c(0, .01))) +
    scale_x_discrete(labels = x_labels) + 
    theme(axis.title.x = element_blank()) + 
    theme(axis.text.x = element_text(angle = Xaxis_rot, hjust = 1, size = xaxis_tick_label_size)) +
    theme(axis.text.y = element_text(size = yaxis_tick_label_size)) +   
    theme(axis.title.y = element_text(size = yaxis_label_size)) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, )) + #face = "bold")) +
    theme(plot.margin = margin(0, 0, 0.3, 0, "cm")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.border = element_rect(fill = NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    
    theme(legend.text = element_markdown(size = 12),
          legend.title = element_blank(),
          legend.position = "right",
          legend.justification = "left",
          legend.direction = "vertical") +
    
    ylab("Relative abundance (%)") +
    facet_wrap(wrap_group, scales = "free_x", nrow = 1, strip.position = "bottom") +
    theme(panel.spacing = unit(0.25, units = "cm"),
          strip.placement = "outside",                      # Place facet labels outside x axis labels.
          strip.background = element_rect(fill = "white"),  # Make facet label background white.
          axis.title = element_blank(),               # Remove x and y axis titles.
          strip.text.x = element_text(size = xaxis_tick_label_size),
    ) +
    ggtitle(plt_title)
  
  return(ha_Plot)
}

ordered_species_guides <- function(specieslist, thresh_str) {
  # vector to order the data frame, to make "other species" the last element of the legend
  unique_sp <- unique(specieslist)
  unique_sp <- unique_sp[unique_sp != thresh_str]
  ordered_species = append(unique_sp, thresh_str)
  # use labels for element_markdown to use italics for species names (but not for "Other species")
  ordered_sp_labels <- as.vector(
    sapply(
      unique_sp,
      FUN = function(X) if (grepl(" sp.", X)) paste0("*", strsplit(X, " sp."), "* sp.") else paste0("*", X, "*")
    )
  )
  # handle subspecies names
  ordered_sp_labels <- as.vector(
    sapply(
      ordered_sp_labels,
      FUN = function(X) if (grepl(" subsp.", X)) paste0(unlist(strsplit(X, " subsp. "))[1] , "* subsp. *", unlist(strsplit(X, " subsp. "))[-1]) else X
    )
  )
  ordered_sp_labels = append(ordered_sp_labels, thresh_str)
  
  return(list(ordered_species, ordered_sp_labels))
}

draw_la_stacked_plot <- function(df, sampleid_vector, 
                                 abund_thresh = 1.0, sp_res_th = 0.6,
                                 x_labels, yaxis_tick_label_size = 10,
                                 yaxis_label_size = 10, xaxis_tick_label_size = 12,
                                 Xaxis_rot = 0, ordered_sp_labels, species_colors,
                                 wrap_group, x_factor,
                                 ordered_species) {
  
  plt_title <- paste(c("<", abund_thresh, "% average abundance"), collapse = " ")
  la_Plot <- ggplot(df %>%
                      mutate(Species = factor(Species, levels = rev(ordered_species))),
                    aes(x = factor(df[[x_factor]], levels = sampleid_vector),
                        y = Abundance, fill = Species)) + 
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = species_colors, labels = rev(ordered_sp_labels)) +
    scale_y_continuous(breaks = seq(0, 100, by = 2.0), expand = expansion(mult = c(0, .01))) +
    scale_x_discrete(labels = x_labels) + 
    theme(axis.title.x = element_blank()) + 
    theme(axis.text.x = element_text(angle = Xaxis_rot, hjust = 1, size = xaxis_tick_label_size)) +
    theme(axis.text.y = element_text(size = yaxis_tick_label_size)) +   
    theme(axis.title.y = element_text(size = yaxis_label_size)) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, )) + #face = "bold")) +
    theme(plot.margin = margin(0, 0, 0.3, 0, "cm")) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.border = element_rect(fill = NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    theme(legend.text = element_markdown(size = 12),
          legend.title = element_blank(),
          legend.position = "right",
          legend.justification = "left",
          legend.direction = "vertical") +
    ylab("Relative abundance (%)") +
    facet_wrap(wrap_group, scales = "free_x", nrow = 1, strip.position = "bottom") +
    theme(panel.spacing = unit(0.25, units = "cm"),
          strip.placement = "outside",                      # Place facet labels outside x axis labels.
          strip.background = element_rect(fill = "white"),  # Make facet label background white.
          axis.title = element_blank(),               # Remove x and y axis titles.
          strip.text.x = element_text(size = xaxis_tick_label_size)) + 
    ggtitle(plt_title)
  
  return(la_Plot)
}