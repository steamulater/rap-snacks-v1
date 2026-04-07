# bar_6/visualize.pml
# TIM barrel — all three buckets hit known folds; concordance top_prob=0.997
# FoldSeek hits — concordance: 2596  native_ala: 1565  free_design: 1922

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_6

load concordance_model_0.pdb, bar_6_concordance
color cyan, bar_6_concordance
show cartoon, bar_6_concordance
hide lines, bar_6_concordance
load native_ala_model_0.pdb, bar_6_native_ala
color tv_orange, bar_6_native_ala
show cartoon, bar_6_native_ala
hide lines, bar_6_native_ala
load free_design_model_0.pdb, bar_6_free_design
color purple, bar_6_free_design
show cartoon, bar_6_free_design
hide lines, bar_6_free_design

align bar_6_native_ala, bar_6_concordance
align bar_6_free_design, bar_6_concordance

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

# === What to look for in bar_6 ===
# TIM barrel — all three buckets hit known folds; concordance top_prob=0.997
# bar_6 concordance is a TIM barrel-like fold (2596 hits, top_prob=0.997)
# bar_6 free_design converges to PII signalling protein (prob=1.000)
# All three buckets are structurally similar — lyric encoded a real fold
# Try: align bar_6_free_design, bar_6_concordance  → should show low RMSD

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
