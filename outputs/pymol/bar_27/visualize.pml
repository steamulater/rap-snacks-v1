# bar_27/visualize.pml
# native_ala 0 hits across all DBs — Ganja Burn lyric AAs encode novel fold
# FoldSeek hits — concordance: 1344  native_ala: 0  free_design: 1440

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_27

load concordance_model_0.pdb, bar_27_concordance
color cyan, bar_27_concordance
show cartoon, bar_27_concordance
hide lines, bar_27_concordance
load native_ala_model_0.pdb, bar_27_native_ala
color tv_orange, bar_27_native_ala
show cartoon, bar_27_native_ala
hide lines, bar_27_native_ala
load free_design_model_0.pdb, bar_27_free_design
color purple, bar_27_free_design
show cartoon, bar_27_free_design
hide lines, bar_27_free_design

align bar_27_native_ala, bar_27_concordance
align bar_27_free_design, bar_27_concordance

# Representation options — uncomment to switch
# show surface, all    # surface view (slow)
# show sticks, all     # all-atom view
# set transparency, 0.5, all

# pLDDT colouring by b-factor (Boltz stores pLDDT per-residue in B column)
# spectrum b, blue_white_red, all, minimum=0, maximum=1

# Show as rainbow from N→C terminus
# spectrum count, rainbow, all

zoom all
orient all

# === What to look for in bar_27 ===
# native_ala 0 hits across all DBs — Ganja Burn lyric AAs encode novel fold
# bar_27 HEADLINE: native_ala has 0 FoldSeek hits across ALL databases
# Ganja Burn lyric amino acids (BJOZXU→Ala) encode a structurally novel fold
# concordance 1344 hits, free_design 1440 hits — but the lyric text itself is novel
# Compare native_ala shape to concordance — distinct topology expected
# Try: show surface, bar_27_native_ala  → inspect the novel fold

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
