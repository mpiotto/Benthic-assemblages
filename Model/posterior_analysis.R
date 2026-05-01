#=============================================================================
# Posterior analysis for the JSDM (Antarctic benthic assemblages)
#-----------------------------------------------------------------------------

# Outline:
#   0.  Setup, theme, helpers, species-level metadata (group)
#   1.  Coefficient extraction (betas)
#   2.  species-specific responses to depth 
#   3.  Predictions over bathymetry (continuous and discrete)
#   4.  Occupancy and conditional cover plots (species level)
#   5.  Dominant species summaries
#   6.  Taxa-level (morphological group) occupancy and cover
#   7.  Sediment-relevant trait groups plot
#   8.  Trait effects (z / gamma)
#   9.  Trait-response scatter plots (feeding, habitat, reproductive)
#  10.  Gamma slope effects (bathymetry)
#  11.  Rho parameter
#  
#
# Figure references (paper / supplementary):
#   Fig.   2A   -> Section 6.1   (dominant taxa absolute occupancy)
#   Fig.   2B   -> Section 6.2   (taxa relative occupancy)
#   Fig.   2C   -> Section 6.3   (taxa conditional cover)
#   Fig.   3A   -> Section 9.1   (feeding strategy boxplot, occupancy)
#   Fig.   3B   -> Section 10    (feeding strategy gamma, occupancy)
#   Fig.   3C   -> Section 9.2   (habitat scatter, occupancy)
#   Fig.   3D   -> Section 10    (habitat gamma, occupancy)
#   Fig.   3E   -> Section 9.2   (reproductive scatter, occupancy)
#   Fig.   3F   -> Section 10    (reproductive gamma, occupancy)
#   Fig.   4    -> Section 7     (sediment-relevant trait groups)
#   Fig.   5A   -> Section 4.1   (dominant species occupancy curves)
#   Fig.   SA5A -> Section 9.1   (feeding strategy boxplot, cover|present)
#   Fig.   SA5B -> Section 10    (feeding strategy gamma, cover|present)
#   Fig.   SA5C -> Section 9.2   (habitat scatter, cover|present)
#   Fig.   SA5D -> Section 10    (habitat gamma, cover|present)
#   Fig.   SA5E -> Section 9.2   (reproductive scatter, cover|present)
#   Fig.   SA5F -> Section 10    (reproductive gamma, cover|present)
#   Fig.   SA6  -> Section 11    (rho posterior)
#   Fig.   SA7  -> Section 2     (species-specific responses to depth)
#   Fig.   SA8  -> Section 4.2   (dominant species conditional cover)
#   Fig.   SA9  -> Section 4.1   (subdominant ascidians, occupancy)
#   Fig.   SA10 -> Section 4.2   (subdominant ascidians, cover|present)
#=============================================================================

# ---- 0. Setup ------------------------------------------------------------

rm(list = ls())

library(rstan)
library(matrixStats)
library(tidyverse)   # dplyr, tidyr, ggplot2, tibble, purrr
library(viridis)

WORK_DIR <- ""
setwd(WORK_DIR)

load("inputdata.RData")
load("model_fit.RData")

draws <- rstan::extract(fit)


# ---- 0.1  Reusable theme and palettes ------------------------------------

theme_pub <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      axis.text        = element_text(colour = "black"),
      axis.line        = element_line(colour = "black"),
      legend.title     = element_blank(),
      panel.grid.minor = element_blank()
    )
}

YEAR_COLS  <- c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF")
GROUP_COLS <- c("lightblue", "red", "grey80", "orange", "darkorange3",
                "skyblue", "brown", "black", "hotpink", "purple",
                "blue", "yellow")


# ---- 0.2  Generic helpers ------------------------------------------------

col_summary <- function(M, probs = c(0.025, 0.975)) {
  data.frame(
    mean = colMeans(M, na.rm = TRUE),
    L    = colQuantiles(M, probs = probs[1], na.rm = TRUE),
    U    = colQuantiles(M, probs = probs[2], na.rm = TRUE)
  )
}
row_summary <- function(M, probs = c(0.025, 0.975)) {
  data.frame(
    mean = rowMeans(M, na.rm = TRUE),
    L    = rowQuantiles(M, probs = probs[1], na.rm = TRUE),
    U    = rowQuantiles(M, probs = probs[2], na.rm = TRUE)
  )
}


# ---- 0.3  Species-level metadata -----------------------------------------
# Cross feeding strategy with habitat (binarized to two levels).
#
# Final assignment used downstream:
#   G1 = infaunal filter-feeder
#   G2 = short epifaunal filter-feeder
#   G3 = tall epifaunal filter-feeder
#   G4 = predator
#   G5 = opportunistic

trait_df <- data.frame(species = species,
                       fed   = tr$fstrategy,
                       hab   = tr$habitat,
                       stringsAsFactors = FALSE)
trait_df$hab[trait_df$hab > 2] <- 2

combos        <- unique(trait_df[, c("fed", "hab")])
combos$group  <- paste0("G", seq_len(nrow(combos)))

trait_df <- merge(trait_df, combos, by = c("fed", "hab"))
trait_df <- trait_df[order(trait_df$species), ]

trait_df <- trait_df %>%
  mutate(group = case_when(
    group == "G3" ~ "G6",
    group == "G2" ~ "G3",
    group == "G1" ~ "G2",
    group == "G5" ~ "G1",
    group == "G4" ~ "G5",
    TRUE          ~ group)) %>%
  mutate(group = case_when(
    group == "G6" ~ "G4",
    TRUE          ~ group))

# Restore original species order (the rest of the script assumes alignment with `species`)
trait_df <- trait_df[match(species, trait_df$species), ]
group    <- trait_df$group

GROUP_NICE <- c(G1 = "infaunal filter-feeder",
                G2 = "short epifaunal filter-feeder",
                G3 = "tall epifaunal filter-feeder",
                G4 = "predator",
                G5 = "opportunistic")


# ---- 1. Coefficient extraction -------------------------------------------
# Layout: 16 coefficients per species, 1-8 presence (psi), 9-16 cover (mu).
#   
# pres    cover        parameter
#   1 and  9 : intercept (~15 m, year 1994)
#   2 and 10 : intercept differential, 1998
#   3 and 11 : intercept differential, 2010
#   4 and 12 : intercept differential, 2022
#   5 and 13 : slope w.r.t. bathymetry
#   6 and 14 : slope differential, 1998
#   7 and 15 : slope differential, 2010
#   8 and 16 : slope differential, 2022

COEF_NAMES <- c(
  "94_psi", "98_psi", "10_psi", "22_psi",
  "deep_psi", "d98_psi", "d10_psi", "d22_psi",
  "94_mu",  "98_mu",  "10_mu",  "22_mu",
  "deep_mu", "d98_mu", "d10_mu", "d22_mu"
)
K_PRES  <- 8
K_COV   <- 8
K_TOTAL <- K_PRES + K_COV
stopifnot(length(COEF_NAMES) == K_TOTAL)

# draws$betas is [n_draws, nsp * K_TOTAL] in column-major order: cols 1..nsp
# correspond to coef 1 across species, nsp+1..2*nsp to coef 2, etc.
# Reshape into a 3-D array [draw, species, coef].
n_draws   <- nrow(draws$betas)
betas_arr <- array(draws$betas, dim = c(n_draws, nsp, K_TOTAL))
dimnames(betas_arr) <- list(NULL, species, COEF_NAMES)

# Convenience views:
#   post_betas[[k]] : n_draws x nsp  (one matrix per coefficient)
#   species_B[[j]]    : n_draws x K     (one matrix per species)
post_betas <- lapply(seq_len(K_TOTAL), function(k) betas_arr[, , k])
names(post_betas) <- COEF_NAMES

species_B <- lapply(seq_len(nsp), function(j) betas_arr[, j, ])
names(species_B) <- species


# ---- 2. Species-specific responses to depth (b_deep) ---------------------
# Year-specific depth slopes per species: deep_psi (+ d{yy}_psi for yy != 94).

depth_slopes_by_year <- list(
  "1994" = post_betas[["deep_psi"]],
  "1998" = post_betas[["deep_psi"]] + post_betas[["d98_psi"]],
  "2010" = post_betas[["deep_psi"]] + post_betas[["d10_psi"]],
  "2022" = post_betas[["deep_psi"]] + post_betas[["d22_psi"]]
)

b_deep <- bind_rows(
  lapply(names(depth_slopes_by_year), function(yr) {
    s <- col_summary(depth_slopes_by_year[[yr]])
    data.frame(year = yr, species = species, group = group, s)
  })
)

# Fig. SA7 
p_b_deep <- b_deep %>%
  filter(group %in% c("G2", "G3")) %>%
  mutate(group = factor(group, levels = c("G2", "G3"))) %>%
  ggplot(aes(x = species, y = mean, color = year)) +
  geom_linerange(aes(ymin = L, ymax = U),
                 position = position_dodge(width = 0.6),
                 alpha = 0.65, size = 3) +
  geom_point(size = 4.5, alpha = 0.7,
             position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "darkgrey", size = 1) +
  labs(x = "Species", y = "Species responses to depth") +
  scale_colour_manual(values = viridis(4)) +
  theme_pub() +
  theme(axis.text.x  = element_text(angle = 90, size = 18),
        axis.text.y  = element_text(size = 18),
        axis.title   = element_text(size = 20),
        legend.text  = element_text(size = 16),
        strip.text   = element_text(size = 16),
        legend.position = "bottom",
        panel.spacing.y = unit(0.5, "cm"))
print(p_b_deep)


# ---- 3. Predictions over bathymetry --------------------------------------
# For each species and each bathymetry value, simulate marginal occupancy
# (averaging over the transect random effect via a quadrature of `draws_z`)
# and conditional cover (Beta mean), then unconditional cover.
# One simulator, called twice: continuous gradient + discrete sampled depths.
# Predictions are kept in memory only -- not written to disk.

draws_z     <- qnorm(p = ppoints(500))
sigma_ranef <- sqrt(draws$sigma_tr^2 + draws$sigma_tr_sp^2)
n_iter      <- nrow(species_B[[1]])

# Years and the coefficient indices that build them up.
# Each entry: c(int, slope, d_int_or_NA, d_slope_or_NA)
year_coef_idx <- list(
  "94" = c(int = 1, slope = 5, d_int = NA, d_slope = NA),
  "98" = c(int = 1, slope = 5, d_int = 2,  d_slope = 6),
  "10" = c(int = 1, slope = 5, d_int = 3,  d_slope = 7),
  "22" = c(int = 1, slope = 5, d_int = 4,  d_slope = 8)
)

simulate_species <- function(b_post, bathy, sigma_ranef, draws_z) {
  # b_post : n_iter x K_TOTAL (posterior of one species, see species_B[[j]])
  # Returns a named list of 12 matrices [length(bathy) x n_iter]:
  # Psi94..21, mu94..21, cov94..21
  n_iter <- nrow(b_post)
  nb     <- length(bathy)

  out_names <- c(paste0("Psi", names(year_coef_idx)),
                 paste0("mu",  names(year_coef_idx)),
                 paste0("cov", names(year_coef_idx)))
  out <- replicate(length(out_names),
                   matrix(NA_real_, nb, n_iter), simplify = FALSE)
  names(out) <- out_names

  for (j in seq_len(n_iter)) {
    eps <- draws_z * sigma_ranef[j]              # transect ranef draws
    for (k in seq_along(year_coef_idx)) {
      idx_p <- year_coef_idx[[k]]
      idx_m <- idx_p + 8                          # cover block offset

      # Presence: marginalize random effect by averaging plogis(eta + eps)
      eta_p <- b_post[j, idx_p["int"]] +
               (if (is.na(idx_p["d_int"]))   0 else b_post[j, idx_p["d_int"]]) +
               (b_post[j, idx_p["slope"]] +
                  (if (is.na(idx_p["d_slope"])) 0 else b_post[j, idx_p["d_slope"]])) * bathy
      out[[k]][, j] <- vapply(eta_p,
                              function(e) mean(plogis(e + eps)),
                              numeric(1))

      # Conditional cover: deterministic (no random effect in cover)
      eta_m <- b_post[j, idx_m["int"]] +
               (if (is.na(idx_m["d_int"]))   0 else b_post[j, idx_m["d_int"]]) +
               (b_post[j, idx_m["slope"]] +
                  (if (is.na(idx_m["d_slope"])) 0 else b_post[j, idx_m["d_slope"]])) * bathy
      out[[4 + k]][, j] <- plogis(eta_m)
    }
  }
  # Unconditional cover
  for (k in 1:4) out[[8 + k]] <- out[[k]] * out[[4 + k]]
  out
}

simulate_all <- function(bathy) {
  setNames(
    lapply(seq_along(species), function(i) {
      cat(sprintf("Species %d / %d\n", i, nsp))
      simulate_species(species_B[[i]], bathy, sigma_ranef, draws_z)
    }),
    species
  )
}

# Continuous bathymetric gradient
bathy_cont <- seq(sdepth_min, sdepth_max, length.out = 10)
Preds      <- simulate_all(bathy_cont)
#load("Preds.Rdata")
# run simulations or continue with our Preds

# Discrete sampled depths (15, 20, 25, 30 m, standardized like X)
bathy_disc <- c(15, 20, 25, 30) - depth_avg / depth_sd
Preds_hist <- simulate_all(bathy_disc)
#load("Preds_disc.Rdata")
# run simulations or continue with our Preds_hist

# Stack predictions into a tidy data frame
preds_to_df <- function(Preds, slots, param_label, bathy_grid) {
  # slots: integer vector of length 4 selecting four matrices per species
  # (1:4 occupancy, 5:8 cover|present, 9:12 unconditional cover).
  years <- c("1994", "1998", "2010", "2022")
  bind_rows(lapply(seq_along(species), function(i) {
    bind_rows(lapply(seq_along(slots), function(k) {
      M <- Preds[[i]][[ slots[k] ]]
      data.frame(
        species  = species[i],
        morpho = morpho[morph_id[i]],
        type   = type[[i]],
        param  = param_label,
        year   = years[k],
        row_summary(M),
        bathy  = bathy_grid
      )
    }))
  }))
}


# ---- 4. Dominant species occupancy and conditional cover plots -----------

OUT1 <- preds_to_df(Preds, 1:4, "occupancy", bathy_cont)   # occupancy
OUT2 <- preds_to_df(Preds, 5:8, "coverp",    bathy_cont)   # cover | present

# Back-transform standardized bathymetry to depth (m)
OUT1$sbathy <- OUT1$bathy * depth_sd + depth_avg
OUT2$sbathy <- OUT2$bathy * depth_sd + depth_avg


# Species labels
SP_LABELS <- c(MOLGULA      = "Molgula pedunculata",
               MALACOBELEMNON = "Malacobelemnon daytoni",
               LATERNULA    = "Laternula elliptica",
               SEROLIS      = "Paraserolis spp.",
               CNEMIDOCARPA = "Cnemidocarpa verrucosa",
               ASCIDIA      = "Ascidia challengeri",
               CORELLA      = "Corella antarctica")

relabel_species <- function(df, ord) {
  df %>% mutate(species = recode(species, !!!SP_LABELS),
                species = factor(species, levels = ord))
}

# 4.1 Occupancy curves per species
plot_occ_curves <- function(df, ord) {
  df %>% relabel_species(ord) %>%
    ggplot(aes(x = sbathy, y = mean, col = year)) +
    geom_ribbon(aes(ymin = L, ymax = U, fill = year), alpha = 0.3) +
    geom_line(size = 1) +
    scale_colour_manual(values = YEAR_COLS) +
    scale_fill_manual(values   = YEAR_COLS) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    facet_grid(rows = vars(species), cols = vars(year)) +
    labs(x = "Depth (m)", y = "Occupancy") +
    theme_pub() +
    theme(axis.text  = element_text(size = 18),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          strip.text  = element_text(size = 16),
          legend.position = "bottom",
          panel.spacing.y = unit(0.5, "cm"))
}

# Fig. 5A 
ord_com_occ <- c("Molgula pedunculata", "Malacobelemnon daytoni",
                 "Laternula elliptica", "Paraserolis spp.")
print(plot_occ_curves(OUT1 %>% filter(type == "common"), ord_com_occ))

# Fig. SA9 
ord_asc_occ <- c("Cnemidocarpa verrucosa", "Ascidia challengeri", "Corella antarctica")
print(plot_occ_curves(OUT1 %>% filter(species %in% c("CNEMIDOCARPA", "ASCIDIA", "CORELLA")),
                      ord_asc_occ))


# 4.2 Conditional-cover plots per species (at discrete depths)
coverp_df <- preds_to_df(Preds_hist, 5:8, "coverp", c(15, 20, 25, 30)) %>%
  filter(!(year == "1998" & bathy %in% c(15, 25))) %>%
  rename(depth = bathy)

plot_cov_bars <- function(df, ord) {
  df %>% relabel_species(ord) %>%
    ggplot(aes(x = factor(depth), y = mean, col = year)) +
    geom_linerange(aes(ymin = L, ymax = U),
                   position = position_dodge(width = 0.5),
                   alpha = 0.8, size = 2) +
    geom_point(size = 4.5, alpha = 0.9,
               position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = YEAR_COLS) +
    facet_grid(rows = vars(species)) +
    labs(x = "Depth (m)", y = "Percentage conditional cover") +
    theme_pub(base_size = 13) +
    theme(legend.position = "bottom")
}

# Fig. SA8 
ord_com_cov <- c("Molgula pedunculata", "Malacobelemnon daytoni",
                 "Laternula elliptica", "Paraserolis spp.")
print(plot_cov_bars(coverp_df %>% filter(type == "common"), ord_com_cov))

# Fig. SA10 
ord_asc_cov <- c("Ascidia challengeri", "Cnemidocarpa verrucosa", "Corella antarctica")
print(plot_cov_bars(coverp_df %>%
                      filter(species %in% c("CNEMIDOCARPA", "ASCIDIA", "CORELLA")),
                    ord_asc_cov))


# ---- 5. Dominant species summaries ---------------------------------------
# Slot indices in Preds_hist:
#   1..4  : occupancy (Psi) for 1994, 1998, 2010, 2022
#   5..8  : cover|present (mu) for the same four years
# Depth indices: 1=15 m, 2=20 m, 3=25 m, 4=30 m

summarise_at_depth <- function(sp, slot, depth_idx) {
  v <- Preds_hist[[sp]][[slot]][depth_idx, ]
  data.frame(mean = mean(v),
             L    = quantile(v, 0.025),
             U    = quantile(v, 0.975))
}

ratio_summary <- function(sp_num, slot_num, sp_den, slot_den, depth_ids) {
  M_num <- Preds_hist[[sp_num]][[slot_num]][depth_ids, , drop = FALSE]
  M_den <- Preds_hist[[sp_den]][[slot_den]][depth_ids, , drop = FALSE]
  r     <- colMeans(M_num / M_den)
  data.frame(mean = mean(r),
             L    = quantile(r, 0.025),
             U    = quantile(r, 0.975))
}

# Examples:
#   summarise_at_depth("MOLGULA", 1, 3)               # occ Molgula 1994 at 25 m
#   ratio_summary("MOLGULA", 1, "MOLGULA", 2, 3:4)    # 1994/1998 occ avg 25-30 m

# Dominant ascidians: mean occupancy at 30 m, 1994 vs 2021
asc_top4 <- c("MOLGULA", "CORELLA", "CNEMIDOCARPA", "ASCIDIA")
ma94 <- colMeans(do.call(rbind, lapply(asc_top4,
                                       function(s) Preds_hist[[s]][[1]][4, ])))
ma21 <- colMeans(do.call(rbind, lapply(asc_top4,
                                       function(s) Preds_hist[[s]][[4]][4, ])))
asc_occ_ratio <- list(mean = mean(ma94 / ma21),
                      ci   = quantile(ma94 / ma21, c(0.025, 0.975)))


# ---- 6. Taxa-level (morphological group) occupancy and cover -------------

# Probability that AT LEAST ONE species in a group is present:
#   Occ = 1 - prod(1 - psi_j)
process_group <- function(group_name, preds, slots = 1:4) {
  members <- preds[names(preds) == group_name]
  if (length(members) == 0) stop("No species in group ", group_name)
  combined <- lapply(slots, function(s) {
    mats <- lapply(members, function(sp) 1 - sp[[s]])
    Reduce("*", mats)
  })
  1 - do.call(rbind, combined)
}

# 16 rows per taxon (4 years x 4 depths)
build_focc <- function(M) {
  data.frame(row_summary(M),
             year  = rep(c(1994, 1998, 2010, 2022), each = 4),
             depth = rep(c(15, 20, 25, 30), times = 4))
}

# 6.1 Higher-taxon occupancy
preds_morpho <- setNames(Preds_hist, morpho[morph_id])
taxa_levels  <- c("ascidiacea",   "porifera",    "asteroidea",
                  "isopoda",      "gastropoda",  "bivalvia",
                  "ophiuroidea",  "ctenophora",  "pennatulaceae",
                  "nematoda",     "echinoidea",  "actiniaria")
Occ_taxa <- setNames(lapply(taxa_levels, process_group, preds = preds_morpho),
                     taxa_levels)

focc_g <- bind_rows(lapply(names(Occ_taxa), function(nm) {
  cbind(build_focc(Occ_taxa[[nm]]), morfo = nm)
}))

# Fig. 2A 
p_dom <- focc_g %>%
  filter(morfo %in% c("ascidiacea", "pennatulaceae", "isopoda", "bivalvia"),
         !(year == 1998 & depth %in% c(15, 25))) %>%
  ggplot(aes(x = factor(depth), y = mean, color = morfo, fill = morfo)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.85),
           alpha = 0.8) +
  geom_errorbar(aes(ymin = L, ymax = U),
                position = position_dodge(width = 0.85),
                width = 0.25, size = 0.8) +
  scale_fill_manual(values   = c("red", "orange", "black", "blue")) +
  scale_colour_manual(values = c("red", "orange", "black", "blue")) +
  facet_grid(rows = vars(year)) +
  labs(x = "Depth (m)", y = "Occupancy") +
  theme_pub(base_size = 13) +
  theme(legend.position = "none", strip.text = element_text(size = 14))
print(p_dom)


# 6.2 Relative occupancy (each taxon / sum of all taxa)
sum_Occ   <- Reduce("+", Occ_taxa)
re_focc_g <- bind_rows(lapply(names(Occ_taxa), function(nm) {
  cbind(build_focc(Occ_taxa[[nm]] / sum_Occ), morfo = nm)
}))

# Fig. 2B 
p_rel <- re_focc_g %>%
  filter(!(year == 1998 & depth %in% c(15, 25))) %>%
  ggplot(aes(x = factor(depth), y = mean, fill = morfo)) +
  facet_grid(cols = vars(year)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = GROUP_COLS) +
  labs(x = "Depth (m)", y = "Relative occupancy") +
  theme_pub() +
  theme(legend.position = "none", strip.text = element_text(size = 16))
print(p_rel)


# 6.3 Taxa-level conditional cover (weighted by within-group relative occ)
occ_yearly  <- lapply(Preds_hist, function(s) do.call(rbind, s[1:4]))
covp_yearly <- lapply(Preds_hist, function(s) do.call(rbind, s[5:8]))
names(occ_yearly)  <- morpho[morph_id]
names(covp_yearly) <- morpho[morph_id]

# Sum of within-group occupancies (one matrix per group)
group_occ_sum <- tapply(seq_along(occ_yearly), names(occ_yearly),
                        function(idx) Reduce("+", occ_yearly[idx]))

# Relative within-group occupancy per species
rel_occ <- mapply(function(M, gname) M / group_occ_sum[[gname]],
                  occ_yearly, names(occ_yearly), SIMPLIFY = FALSE)

# Group-level conditional cover = sum_j (rel_occ_j * cov|p_j)
weighted_cov <- mapply(`*`, covp_yearly, rel_occ, SIMPLIFY = FALSE)
group_cov    <- tapply(seq_along(weighted_cov), names(weighted_cov),
                       function(idx) Reduce("+", weighted_cov[idx]))

cov_mat <- do.call(rbind, group_cov)
covp_taxa <- data.frame(
  mean   = rowMeans(cov_mat),
  morpho = rep(names(group_cov), each = nrow(group_cov[[1]])),
  year   = rep(c("1994", "1998", "2010", "2022"), each = 4),
  depth  = rep(c("15", "20", "25", "30"), length.out = nrow(cov_mat))
) %>% filter(!(year == "1998" & depth %in% c("15", "25")))

# Fig. 2C 
p_cov_taxa <- ggplot(covp_taxa, aes(x = depth, y = mean * 100, fill = morpho)) +
  facet_grid(cols = vars(year)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = GROUP_COLS) +
  labs(x = "Depth (m)", y = "Percentage conditional cover") +
  theme_pub() +
  theme(legend.position = "none", strip.text = element_text(size = 16))
print(p_cov_taxa)


# ---- 7. Sediment-relevant trait groups  ---------------------------

preds_grp     <- setNames(Preds_hist, group)
groups_to_run <- intersect(c("G1", "G2", "G3", "G4", "G5"), unique(group))
Occ_groups    <- setNames(lapply(groups_to_run, process_group, preds = preds_grp),
                          groups_to_run)

focc_g_traits <- bind_rows(lapply(names(Occ_groups), function(g) {
  cbind(build_focc(Occ_groups[[g]]), group = g)
})) %>%
  filter(!(year == 1998 & depth %in% c(15, 25))) %>%
  mutate(group = factor(GROUP_NICE[group], levels = GROUP_NICE))

# Fig. 4 
p_grp <- ggplot(focc_g_traits,
                aes(x = factor(depth), y = mean,
                    color = group, fill = group)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8), alpha = 0.75) +
  geom_errorbar(aes(ymin = L, ymax = U),
                position = position_dodge(width = 0.8),
                width = 0.25, size = 0.8) +
  scale_fill_manual(values   = viridis(5)) +
  scale_colour_manual(values = viridis(5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  facet_grid(rows = vars(factor(year))) +
  labs(x = "Depth (m)", y = "Occupancy") +
  guides(color = guide_legend(ncol = 2),
         fill  = guide_legend(ncol = 2)) +
  theme_pub() +
  theme(legend.position = "bottom",
        panel.spacing.y = unit(0.5, "cm"))
print(p_grp)


# ---- 8. Trait effects (z / gamma) ----------------------------------------
# draws$z is [n_iter, K_TOTAL * L_J]. Within each covariate block, the L_J
# columns are the trait coefficients. First K_PRES*L_J cols correspond to
# psi covariates, the next K_COV*L_J to mu.
#
# Trait columns (order TT was built in):
#   1 = intercept, 2 = predator, 3 = opportunistic,
#   4 = rstrategy (continuous), 5 = habitat (continuous)

L_J <- 5
n_z <- nrow(draws$z)

z_psi <- draws$z[, 1:(K_PRES * L_J)]
z_mu  <- draws$z[, (K_PRES * L_J + 1):(K_TOTAL * L_J)]

# Reshape into [iter, trait, covariate]
arr_psi <- array(z_psi, dim = c(n_z, L_J, K_PRES))
arr_mu  <- array(z_mu,  dim = c(n_z, L_J, K_PRES))

COV_NAMES <- c("y94", "dy98", "dy09", "dy21", "d94", "s98", "s09", "s21")

trait_slice <- function(arr, trait_id) {
  out <- arr[, trait_id, ]
  colnames(out) <- COV_NAMES
  out
}

# Effects on occupancy
int_psi    <- trait_slice(arr_psi, 1)
pred_psi   <- int_psi + trait_slice(arr_psi, 2)
oppor_psi  <- int_psi + trait_slice(arr_psi, 3)
rep_psi    <- trait_slice(arr_psi, 4)
hab_psi    <- trait_slice(arr_psi, 5)
list_gamma_occ <- list(filter = int_psi, predator = pred_psi,
                       opportunistic = oppor_psi)

# Effects on cover|present
int_mu    <- trait_slice(arr_mu, 1)
pred_mu   <- int_mu + trait_slice(arr_mu, 2)
oppor_mu  <- int_mu + trait_slice(arr_mu, 3)
rep_mu    <- trait_slice(arr_mu, 4)
hab_mu    <- trait_slice(arr_mu, 5)
list_gamma_cov <- list(filter = int_mu, predator = pred_mu,
                       opportunistic = oppor_mu)


# ---- 9. Trait scatter plots (depth response by trait) --------------------
# Year-specific depth slopes per species (psi reused from Section 2; mu new).

depth_slopes_mu_by_year <- list(
  "1994" = post_betas[["deep_mu"]],
  "1998" = post_betas[["deep_mu"]] + post_betas[["d98_mu"]],
  "2010" = post_betas[["deep_mu"]] + post_betas[["d10_mu"]],
  "2022" = post_betas[["deep_mu"]] + post_betas[["d22_mu"]]
)

build_response_df <- function(slope_list, trait_col) {
  bind_rows(lapply(names(slope_list), function(yr) {
    s <- col_summary(slope_list[[yr]])
    data.frame(year = yr, species = species, morfo = tr$morfo,
               trait = trait_col, s)
  }))
}

# 9.1 Feeding strategy (categorical) -- boxplot per strategy
plot_fed_box <- function(slope_list, ylab) {
  df <- build_response_df(slope_list, tr$fstrategy) %>%
    rename(fed = trait) %>%
    mutate(fed = recode(fed, "filter" = "filter feeder"),
           fed = factor(fed, levels = c("filter feeder", "predator", "opportunistic")))
  ggplot(df, aes(x = fed, y = mean, color = year, fill = year)) +
    geom_boxplot(alpha = 0.5) +
    geom_hline(yintercept = 0, size = 1.3, linetype = "dashed", color = "darkgrey") +
    scale_colour_manual(values = viridis(4)) +
    scale_fill_manual(values   = viridis(4)) +
    labs(x = "Feeding strategy", y = ylab) +
    theme_pub(base_size = 13) +
    theme(legend.position = "none")
}

# Fig. 3A 
print(plot_fed_box(depth_slopes_by_year,    "Response to depth (occupancy)"))

# Fig. SA5A 
print(plot_fed_box(depth_slopes_mu_by_year, "Response to depth (cover|present)"))


# 9.2 Continuous traits (habitat, reproductive) -- scatter with linear smooth
plot_trait_scatter <- function(slope_list, trait_vec, xlab, ylab) {
  df <- build_response_df(slope_list, trait_vec)
  df$strait <- df$trait * sd(trait_vec) + mean(trait_vec)
  ggplot(df, aes(x = strait, y = mean, color = year)) +
    geom_linerange(aes(ymin = L, ymax = U),
                   position = position_dodge(width = 0.5),
                   alpha = 0.55, size = 0.5) +
    geom_point(size = 4, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, size = 1.3, linetype = "dashed",
               color = "darkgrey") +
    geom_smooth(method = "lm", se = FALSE) +
    scale_colour_manual(values = viridis(4)) +
    labs(x = xlab, y = ylab) +
    theme_pub(base_size = 13) +
    theme(legend.position = "none")
}

# Fig. 3C 
print(plot_trait_scatter(depth_slopes_by_year,    tr$habitat,
                         "Habitat use", "Response to depth (occupancy)"))
# Fig. SA5C
print(plot_trait_scatter(depth_slopes_mu_by_year, tr$habitat,
                         "Habitat use", "Response to depth (cover|present)"))
# Fig. 3E 
print(plot_trait_scatter(depth_slopes_by_year,    tr$rstrategy,
                         "Reproductive strategy", "Response to depth (occupancy)"))
# Fig. SA5E 
print(plot_trait_scatter(depth_slopes_mu_by_year, tr$rstrategy,
                         "Reproductive strategy", "Response to depth (cover|present)"))


# ---- 10. Gamma slope effects (year-by-year depth response) ---------------
# For each strategy / continuous trait, build year-specific intercepts and
# slopes from the 8 gamma covariate columns:
#   cols 1..4 -> year intercepts: y94, y94+dy98, y94+dy10, y94+dy22
#   cols 5..8 -> year bathy slopes: d94, d94+s98, d94+s10, d94+s22

build_gamma_year_slopes <- function(gamma_mat) {
  # gamma_mat: n_iter x 8, columns y94, dy98, dy09, dy21, d94, s98, s10, s22
  out <- matrix(NA_real_, nrow(gamma_mat), 8)
  out[, 1] <- gamma_mat[, 1]
  out[, 2] <- gamma_mat[, 1] + gamma_mat[, 2]
  out[, 3] <- gamma_mat[, 1] + gamma_mat[, 3]
  out[, 4] <- gamma_mat[, 1] + gamma_mat[, 4]
  out[, 5] <- gamma_mat[, 5]
  out[, 6] <- gamma_mat[, 5] + gamma_mat[, 6]
  out[, 7] <- gamma_mat[, 5] + gamma_mat[, 7]
  out[, 8] <- gamma_mat[, 5] + gamma_mat[, 8]
  colnames(out) <- c("y94", "y98", "y10", "y22",
                     "d94", "d98", "d10", "d22")
  out
}

# 10.1 Multi-strategy gamma (categorical traits): one row per strategy
build_gamma_df_multi <- function(gamma_list) {
  bind_rows(lapply(names(gamma_list), function(nm) {
    M <- build_gamma_year_slopes(gamma_list[[nm]])
    data.frame(strategy = nm,
               col_summary(M),
               responseto = colnames(M),
               year   = rep(c("1994", "1998", "2010", "2022"), 2),
               factor = rep(c("year", "depth"), each = 4))
  })) %>%
    mutate(strategy = recode(strategy, filter = "Filter",
                             predator = "Predator",
                             opportunistic = "Opportunistic"),
           strategy = factor(strategy,
                             levels = c("Filter", "Predator", "Opportunistic")))
}

# 10.2 Single-trait gamma (continuous traits): one row per year
build_gamma_df_single <- function(gamma_mat) {
  M <- build_gamma_year_slopes(gamma_mat)
  data.frame(col_summary(M),
             responseto = colnames(M),
             year   = rep(c("1994", "1998", "2010", "2022"), 2),
             factor = rep(c("year", "depth"), each = 4))
}

# Plot for the multi-strategy case (feeding strategy)
plot_gamma_strategy <- function(df) {
  df %>% filter(factor == "depth") %>%
    ggplot(aes(x = factor(year), y = mean, fill = strategy)) +
    geom_hline(yintercept = 0, size = 1.5, linetype = "dashed",
               color = "darkgrey") +
    geom_point(aes(shape = strategy),
               position = position_dodge(width = 0.5),
               size = 5, color = "black", fill = "black") +
    geom_errorbar(aes(ymin = L, ymax = U),
                  position = position_dodge(width = 0.5),
                  width = 0.2, size = 0.75, color = "black") +
    scale_shape_manual(values = c(16, 17, 15)) +
    scale_fill_manual(values  = viridis(3)) +
    labs(x = "Year",
         y = "Effect of trait on species responses to bathy") +
    theme_pub() +
    theme(axis.text  = element_text(size = 17),
          axis.title = element_text(size = 15),
          legend.position    = c(0.01, 0.99),
          legend.justification = c(0, 1),
          legend.key.size    = unit(0.5, "lines"),
          legend.background  = element_rect(fill = "transparent", color = NA),
          legend.direction   = "horizontal")
}

# Plot for the single-trait case (reproductive, habitat)
plot_gamma_single <- function(df) {
  df %>% filter(factor == "depth") %>%
    ggplot(aes(x = factor(year), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey", size = 1) +
    geom_linerange(aes(ymin = L, ymax = U), size = 0.8) +
    geom_point(size = 4) +
    labs(x = "Year",
         y = "Effect of the trait on species responses to depth") +
    theme_pub(base_size = 13)
}

#  Fig. 3B 
print(plot_gamma_strategy(build_gamma_df_multi(list_gamma_occ)))
#  Fig. SA5B 
print(plot_gamma_strategy(build_gamma_df_multi(list_gamma_cov)))
#  Fig. 3D 
print(plot_gamma_single(build_gamma_df_single(hab_psi)))
#  Fig. SA5D 
print(plot_gamma_single(build_gamma_df_single(hab_mu)))
#  Fig. 3F 
print(plot_gamma_single(build_gamma_df_single(rep_psi)))
#  Fig. SA5F 
print(plot_gamma_single(build_gamma_df_single(rep_mu)))


# ---- 11. Rho parameter ---------------------------------------------------

# Fig. SA6 
p_rho <- ggplot(data.frame(rho = draws$rho), aes(x = rho)) +
  geom_density(color = "black", fill = "grey", alpha = 0.65, size = 0.4) +
  geom_vline(xintercept = mean(draws$rho),
             color = "#414487FF", size = 1.2, linetype = "dashed") +
  labs(x = expression(rho), y = "Frequency") +
  theme_pub(base_size = 12)
print(p_rho)

cat(sprintf("rho: mean = %.3f, 95%% CrI = [%.3f, %.3f]\n",
            mean(draws$rho),
            quantile(draws$rho, 0.025),
            quantile(draws$rho, 0.975)))
