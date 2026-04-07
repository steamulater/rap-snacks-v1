# bar_11/visualize.pml
# all three known; no strong structural novelty signal
# FoldSeek hits — concordance: 1649  native_ala: 1039  free_design: 1330

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_11

load concordance_model_0.pdb, bar_11_concordance
color cyan, bar_11_concordance
show cartoon, bar_11_concordance
hide lines, bar_11_concordance
load native_ala_model_0.pdb, bar_11_native_ala
color tv_orange, bar_11_native_ala
show cartoon, bar_11_native_ala
hide lines, bar_11_native_ala
load free_design_model_0.pdb, bar_11_free_design
color purple, bar_11_free_design
show cartoon, bar_11_free_design
hide lines, bar_11_free_design

align bar_11_native_ala, bar_11_concordance
align bar_11_free_design, bar_11_concordance

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

# === What to look for in bar_11 ===
# all three known; no strong structural novelty signal

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
