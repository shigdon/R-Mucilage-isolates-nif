# Kp NIF Regulon Complex Heatmap
library(circlize)
library(ComplexHeatmap)

# 35% NIF Regulon gene model coverage
Kp_0.3_chm_1 <- 
Heatmap(Kp_hmm_0.3_mat,
        name = "HMM Hits",
        col = colormap,
        cluster_columns = F,
        heatmap_legend_param = list(at = seq(0, 50, 2)),
        top_annotation = Kp_Function_cha,
        bottom_annotation = Kp_operon_cha,
        row_names_gp = gpar(fontsize = 1.5))

Kp_0.3_chm_2 <-
  Heatmap(Kp_hmm_0.3_df$Genus, name = "Genus", col = genus_colors_0.3, 
          width = unit(5, "mm"))

Kp_0.3_chm_3 <-
  Heatmap(Kp_hmm_0.3_df$N_ratio, name = "15N/14N", col = inferno(direction = -1, 6), width = unit(5, "mm"))

Kp_0.3_list <- Kp_0.3_chm_1 + Kp_0.3_chm_2 + Kp_0.3_chm_3

draw(Kp_0.3_list, heatmap_legend_side = "right", annotation_legend_side = "right")

## 50 % NIF Regulon gene model coverage
Kp_0.5_chm_1 <- 
  Heatmap(Kp_hmm_0.5_mat,
          name = "HMM Hits",
          col = colormap_Kp_0.5,
          cluster_columns = F,
          heatmap_legend_param = list(at = seq(0, 20, 1)),
          top_annotation = Kp_Function_cha,
          bottom_annotation = Kp_operon_cha,
          row_names_gp = gpar(fontsize = 1.5))

Kp_0.5_chm_2 <-
  Heatmap(Kp_hmm_0.5_df$Genus, name = "Genus", col = genus_colors, 
          width = unit(5, "mm"))

Kp_0.5_chm_3 <-
  Heatmap(Kp_hmm_0.5_df$N_ratio, name = "15N/14N", col = inferno(direction = -1, 6), width = unit(5, "mm"))

Kp_0.5_list <- Kp_0.5_chm_1 + Kp_0.5_chm_2 + Kp_0.5_chm_3

draw(Kp_0.5_list, heatmap_legend_side = "right", annotation_legend_side = "right")

## 0.85% NIF Regulon gene model coverage

draw(Kp_0.9_list, heatmap_legend_side = "right", annotation_legend_side = "right")

# Filtered

## Kp NIF Regulon - 85% coverage - 15N/14N >= 1.2
draw(nonKp_0.9cov_hiN_list, heatmap_legend_side = "right", annotation_legend_side = "left")

## Lactococcus

### 35% coverage

draw(Kp_hmm_0.3_lactococcus_list, heatmap_legend_side = "right", annotation_legend_side = "left")

### 50% coverage
draw(Kp_hmm_0.5_lactococcus_list, heatmap_legend_side = "right", annotation_legend_side = "left")

## 85% coverage
#draw(Kp_hmm_0.9_lactococcus_list, heatmap_legend_side = "right", annotation_legend_side = "left")

## 85% coverage - nonSantos

### nonSantos -> Kp 85% hmm model coverage
Kp_hmm_0.9_nonSantos_chm <- 
draw(Kp_hmm_0.9_nonSantos_list, heatmap_legend_side = "right", annotation_legend_side = "right")

### Santos Positive -> Kp 85% hmm model coverage
Kp_hmm_0.9_santos_pos_chm <- 
draw(Kp_hmm_0.9_santos_pos_list, heatmap_legend_side = "right")#, annotation_legend_side = "left")

### nonSantos -> Kp 85% hmm coverage, ranked by 15N/14N ratio - descending
Kp_hmm_0.9_nonSantos_rnk_chm <- 
  draw(Kp_hmm_0.9_nonSantos_rnk_list, heatmap_legend_side = "right", annotation_legend_side = "right")

### Santos Positive -> 85% hmm coverage, ranked by 15N/14N ratio - descending
Kp_hmm_0.9_santos_pos_rnk_chm <- 
draw(Kp_hmm_0.9_santos_pos_rnk_list, heatmap_legend_side = "right", annotation_legend_side = "right")

## FIG 2A -  Plot Complexheatmap panel of NIF Gene Mining
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(Kp_hmm_0.9_santos_pos_chm, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(Kp_hmm_0.9_santos_pos_rnk_chm, newpage = FALSE)
upViewport()

## FIG 2B -  Plot Complexheatmap panel of NIF Gene Mining
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(Kp_hmm_0.9_nonSantos_chm, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(Kp_hmm_0.9_nonSantos_rnk_chm, newpage = FALSE)
upViewport()

## FIG 2C -  Plot Complexheatmap panel of NIF Gene Mining
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(Kp_hmm_0.9_hemiSantos_chm, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(Kp_hmm_0.9_hemiSantos_rnk_chm, newpage = FALSE)
upViewport()

# alt. nif plots

## Santos Positive Isolates

alt_nif_santos_0.9cov_rnk_chm <- 
draw(alt_nif_santos_0.9cov_rnk_list, heatmap_legend_side = "right", annotation_legend_side = "right")



test_plot
anf_VD_plot <- plot(V_alt_nif, doWeights = FALSE)

draw(Kp_operon_cha, 1:14)
draw(Kp_genus_colors, 1:29)

# alt. nif sourmash compare complex heatmap
draw(chm_smash_alt_nif_list, heatmap_legend_side = "right", annotation_legend_side = "right")


# Mono Prok Plots - Fig 2 + Fig 3

## FIG 2A -  Plot Complexheatmap panel of NIF Gene Mining
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(Kp_hmm_0.9_sp_mono_prok_chm, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(Kp_hmm_0.9_sp_mono_prok_rnk_chm, newpage = FALSE)
upViewport()

## FIG 2B -  Plot Complexheatmap panel of NIF Gene Mining
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(Kp_hmm_0.9_ns_mono_prok_chm, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(Kp_hmm_0.9_ns_mono_prok_rnk_chm, newpage = FALSE)
upViewport()

## FIG 2C -  Plot Complexheatmap panel of NIF Gene Mining
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(Kp_hmm_0.9_hs_mono_prok_chm, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(Kp_hmm_0.9_hs_mono_prok_rnk_chm, newpage = FALSE)
upViewport()

### 15N ranked

## Fig3A
anf_VD_mono_prok_plot <- plot(V_alt_nif_mono_prok, doWeights = FALSE)

#Fig3B
draw(chm_smash_alt_nif_mono_prok_list, heatmap_legend_side = "right", annotation_legend_side = "right")



# PGP_Assay CHM

pgp_master_chm

pgp_master_single_ids_chm
