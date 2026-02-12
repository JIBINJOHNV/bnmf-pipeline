
library(navmix)
library(GenomicSEM)

set.seed(42)

cluster_out = navmix(z_mat)



# ✅ 1️⃣ Variant Assignment Table
z_assign <- cluster_out$fit$z
g_matrix <- cluster_out$fit$g

variant_ids <- rownames(g_matrix)

variant_cluster_table <- data.frame(
  variant = variant_ids,
  cluster = z_assign
)

# Label noise properly
variant_cluster_table$cluster_label <- ifelse(
  variant_cluster_table$cluster == ncol(g_matrix),
  "Noise",
  paste0("Cluster_", variant_cluster_table$cluster)
)

head(variant_cluster_table)


# ✅ 2️⃣ Cluster Size Summary Table
cluster_size_table <- as.data.frame(table(variant_cluster_table$cluster_label))
colnames(cluster_size_table) <- c("Cluster", "Number_of_variants")

cluster_size_table


# ✅ 3️⃣ Trait Profile Table (Cluster Means)
mu_table <- as.data.frame(cluster_out$fit$mu)

mu_table$Trait <- rownames(mu_table)
mu_table <- mu_table[, c(ncol(mu_table), 1:(ncol(mu_table)-1))]

head(mu_table)



# ✅ 4️⃣ Posterior Probability Table (Soft Assignments)
g_table <- as.data.frame(cluster_out$fit$g)
g_table$variant <- rownames(g_table)

g_table <- g_table[, c(ncol(g_table), 1:(ncol(g_table)-1))]

head(g_table)
# High confidence variants only
high_confidence <- g_table[g_table$Cluster.1 > 0.9 |
                           g_table$Cluster.2 > 0.9 |
                           g_table$Cluster.3 > 0.9 |
                           g_table$Cluster.4 > 0.9, ]

# ✅ 5️⃣ Combine Hard + Soft Assignment

combined_table <- merge(
  variant_cluster_table,
  g_table,
  by = "variant"
)

head(combined_table)


# ✅ 6️⃣ Optional: Separate Tables Per Cluster
cluster1_variants <- subset(combined_table, cluster_label == "Cluster_1")
cluster2_variants <- subset(combined_table, cluster_label == "Cluster_2")
cluster3_variants <- subset(combined_table, cluster_label == "Cluster_3")
cluster4_variants <- subset(combined_table, cluster_label == "Cluster_4")
noise_variants    <- subset(combined_table, cluster_label == "Noise")




R <- cor(z_mat, use = "pairwise.complete.obs")



import polars as pl



import polars as pl

def impute_missing_with_ld_proxies_polars(
    target_z_df: pl.DataFrame,
    reference_z_df: pl.DataFrame,
    ld_proxy_df: pl.DataFrame
):
    """
    Independently imputes each trait column by finding the highest R2 proxy 
    present in the reference data for that specific trait.
    """
    # Create a copy to avoid modifying the original
    updated_target_df = target_z_df.clone()
    trait_cols = [c for c in updated_target_df.columns if c != "uniq_id"]
    # 1. Pre-sort proxies by R2 descending.
    # This ensures that when we drop duplicates later, we keep the strongest proxy.
    sorted_proxies = ld_proxy_df.sort("r2", descending=True)
    imputation_log_list = []
    for trait in trait_cols:
        # A. Get only the variants that are actually missing for THIS trait
        missing_variants = updated_target_df.filter(pl.col(trait).is_null()).select("uniq_id")
        if missing_variants.is_empty():
            continue
        # B. Find ALL potential partners for these missing variants
        # Join: [Missing Leads] -> [All Possible Proxies] -> [Reference Z-scores for this trait]
        candidates = (
            missing_variants
            .join(sorted_proxies, left_on="uniq_id", right_on="query_id", how="inner")
            .join(
                reference_z_df.select(["uniq_id", trait]), 
                left_on="partner_id", 
                right_on="uniq_id", 
                how="inner"
            )
            # Crucial: Only keep proxies that actually have data for this trait
            .filter(pl.col(trait).is_not_null())
            # For each missing Lead, keep ONLY the first row (the one with highest R2)
            .unique(subset=["uniq_id"], keep="first")
        )
        if not candidates.is_empty():
            # C. Map the imputed values back into the target column
            # We rename the column to avoid collision during the join
            updated_target_df = updated_target_df.join(
                candidates.select(["uniq_id", trait]).rename({trait: f"{trait}_imp"}),
                on="uniq_id",
                how="left"
            ).with_columns(
                pl.coalesce([pl.col(trait), pl.col(f"{trait}_imp")]).alias(trait)
            ).drop(f"{trait}_imp")
            # D. Log the successful imputations
            log_entry = candidates.select([
                pl.lit(trait).alias("trait"),
                pl.col("uniq_id").alias("original_snp"),
                pl.col("partner_id").alias("proxy_snp"),
                pl.col("r2"),
                pl.col(trait).alias("imputed_value")
            ])
            imputation_log_list.append(log_entry)
    # Combine all logs into a single summary DataFrame
    if imputation_log_list:
        imputation_log_df = pl.concat(imputation_log_list)
    else:
        imputation_log_df = pl.DataFrame(schema={
            "trait": pl.Utf8, "original_snp": pl.Utf8, "proxy_snp": pl.Utf8, 
            "r2": pl.Float64, "imputed_value": pl.Float64
        })
    return updated_target_df, imputation_log_df