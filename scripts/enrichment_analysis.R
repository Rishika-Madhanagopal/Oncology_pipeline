# Load required libraries
library(gprofiler2)
library(dplyr)
library(ggplot2)

# 1. Load top predictive genes 
top_genes <- read.csv("Oncology_pipeline/data/top_predictive_features.csv")$gene

# 2. Run g:Profiler enrichment (human genes)
gostres <- gost(top_genes, organism = "hsapiens")

# 3. Flatten list columns before saving
gost_df <- gostres$result %>%
  mutate(across(where(is.list), ~sapply(., toString)))

# 4. Save results to CSV 
dir.create("result", showWarnings = FALSE)
write.csv(gost_df, "data/gprofiler_enrichment.csv", row.names = FALSE)

# 5. Plot top 10 enriched terms
top_terms <- gost_df %>%
  arrange(p_value) %>%
  slice(1:10) %>%
  mutate(term_name = factor(term_name, levels = rev(term_name)))  # For plot ordering

ggplot(top_terms, aes(x = term_name, y = -log10(p_value), fill = source)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 10 Enriched Pathways (g:Profiler)",
    x = "Pathway",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12)

# 6. Save plot
ggsave("result/gprofiler_enrichment_plot.png", width = 10, height = 6, dpi = 300)
