# bar_9/visualize.pml
# native_ala near-novel (3 hits); free_design → Human Saposin A (prob=0.537)
# FoldSeek hits — concordance: 536  native_ala: 3  free_design: 287

bg_color black
set ray_shadows, 0
set antialias, 2
cd /Users/tamukamartin/software/rap-protein/rap-snacks-v1/outputs/pymol/bar_9

load concordance_model_0.pdb, bar_9_concordance
color cyan, bar_9_concordance
show cartoon, bar_9_concordance
hide lines, bar_9_concordance
load native_ala_model_0.pdb, bar_9_native_ala
color tv_orange, bar_9_native_ala
show cartoon, bar_9_native_ala
hide lines, bar_9_native_ala
load free_design_model_0.pdb, bar_9_free_design
color purple, bar_9_free_design
show cartoon, bar_9_free_design
hide lines, bar_9_free_design

align bar_9_native_ala, bar_9_concordance
align bar_9_free_design, bar_9_concordance

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

# === What to look for in bar_9 ===
# native_ala near-novel (3 hits); free_design → Human Saposin A (prob=0.537)
# bar_9 is the top performer (free_design mean pLDDT 0.839)
# native_ala has only 3 FoldSeek hits — structurally unusual
# free_design converges to Human Saposin A fold (helix bundle, prob=0.537)
# Try: show surface, bar_9_free_design  → should see compact helix bundle
# Try: spectrum b, blue_white_red, bar_9_concordance  → see pLDDT gradient

set_view [\
   1,0,0,\
   0,1,0,\
   0,0,1,\
   0,0,-100,\
   0,0,0,1 ]
