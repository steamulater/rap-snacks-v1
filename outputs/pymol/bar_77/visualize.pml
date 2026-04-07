# bar_77/visualize.pml
# all three known; consistent pattern
# FoldSeek hits — concordance: 1016  native_ala: 1199  free_design: 1210

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_77

load concordance_model_0.pdb, bar_77_concordance
color cyan, bar_77_concordance
show cartoon, bar_77_concordance
hide lines, bar_77_concordance
load native_ala_model_0.pdb, bar_77_native_ala
color tv_orange, bar_77_native_ala
show cartoon, bar_77_native_ala
hide lines, bar_77_native_ala
load free_design_model_0.pdb, bar_77_free_design
color purple, bar_77_free_design
show cartoon, bar_77_free_design
hide lines, bar_77_free_design

align bar_77_native_ala, bar_77_concordance
align bar_77_free_design, bar_77_concordance

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

# === What to look for in bar_77 ===
# all three known; consistent pattern

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
