# bar_3/visualize.pml
# native_ala structurally novel (6 hits); lyric AAs encode unusual fold
# FoldSeek hits — concordance: 1071  native_ala: 6  free_design: 1164

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_3

load concordance_model_0.pdb, bar_3_concordance
color cyan, bar_3_concordance
show cartoon, bar_3_concordance
hide lines, bar_3_concordance
load native_ala_model_0.pdb, bar_3_native_ala
color tv_orange, bar_3_native_ala
show cartoon, bar_3_native_ala
hide lines, bar_3_native_ala
load free_design_model_0.pdb, bar_3_free_design
color purple, bar_3_free_design
show cartoon, bar_3_free_design
hide lines, bar_3_free_design

align bar_3_native_ala, bar_3_concordance
align bar_3_free_design, bar_3_concordance

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

# === What to look for in bar_3 ===
# native_ala structurally novel (6 hits); lyric AAs encode unusual fold
# bar_3 native_ala near-novel (6 hits)
# Ala substitution at BJOZXU positions creates an unusual fold
# Inspect native_ala shape vs concordance and free_design

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
