# bar_0/visualize.pml
# concordance unusual (22 hits); free_design pulls into common fold (1170 hits)
# FoldSeek hits — concordance: 22  native_ala: 273  free_design: 1170

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_0

load concordance_model_0.pdb, bar_0_concordance
color cyan, bar_0_concordance
show cartoon, bar_0_concordance
hide lines, bar_0_concordance
load native_ala_model_0.pdb, bar_0_native_ala
color tv_orange, bar_0_native_ala
show cartoon, bar_0_native_ala
hide lines, bar_0_native_ala
load free_design_model_0.pdb, bar_0_free_design
color purple, bar_0_free_design
show cartoon, bar_0_free_design
hide lines, bar_0_free_design

align bar_0_native_ala, bar_0_concordance
align bar_0_free_design, bar_0_concordance

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

# === What to look for in bar_0 ===
# concordance unusual (22 hits); free_design pulls into common fold (1170 hits)

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
