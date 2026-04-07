# bar_17/visualize.pml
# MPNN improved pLDDT AND moved away from known folds (free_design 56 hits); best outcome
# FoldSeek hits — concordance: 201  native_ala: 1337  free_design: 56

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_17

load concordance_model_0.pdb, bar_17_concordance
color cyan, bar_17_concordance
show cartoon, bar_17_concordance
hide lines, bar_17_concordance
load native_ala_model_0.pdb, bar_17_native_ala
color tv_orange, bar_17_native_ala
show cartoon, bar_17_native_ala
hide lines, bar_17_native_ala
load free_design_model_0.pdb, bar_17_free_design
color purple, bar_17_free_design
show cartoon, bar_17_free_design
hide lines, bar_17_free_design

align bar_17_native_ala, bar_17_concordance
align bar_17_free_design, bar_17_concordance

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

# === What to look for in bar_17 ===
# MPNN improved pLDDT AND moved away from known folds (free_design 56 hits); best outcome
# bar_17 HEADLINE: free_design improved pLDDT AND moved away from known folds
# native_ala 1337 hits vs free_design only 56 hits
# Alignment between concordance and free_design will show large RMSD
# This is the best outcome: foldable + structurally novel
# Try: align bar_17_free_design, bar_17_concordance  → large RMSD expected

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
