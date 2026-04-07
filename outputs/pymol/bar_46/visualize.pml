# bar_46/visualize.pml
# backbone failure — free_design 0 hits is failure not novelty; dropout
# FoldSeek hits — concordance: 636  native_ala: 418  free_design: 0

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_46

load concordance_model_0.pdb, bar_46_concordance
color cyan, bar_46_concordance
show cartoon, bar_46_concordance
hide lines, bar_46_concordance
load native_ala_model_0.pdb, bar_46_native_ala
color tv_orange, bar_46_native_ala
show cartoon, bar_46_native_ala
hide lines, bar_46_native_ala
load free_design_model_0.pdb, bar_46_free_design
color purple, bar_46_free_design
show cartoon, bar_46_free_design
hide lines, bar_46_free_design

align bar_46_native_ala, bar_46_concordance
align bar_46_free_design, bar_46_concordance

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

# === What to look for in bar_46 ===
# backbone failure — free_design 0 hits is failure not novelty; dropout
# bar_46 is a BACKBONE FAILURE — all three structures are disordered
# free_design pLDDT 0.309 even unconstrained — no sequence can fold this backbone
# Included for comparison — this is what failure looks like
# Try: spectrum b, blue_white_red, all  → expect mostly red (low pLDDT)

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
