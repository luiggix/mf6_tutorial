./data -- files needed to recreate the McDonald Valley model
./make_mv.ipynb -- jupyter notebook to create and run the McDonald Valley model with MODFLOW 6 using a triangular mesh.
./optimize_mv.ipynb -- jupyter notebook to optimize pumping rates subject to head constraints.  This notebook requires that the ./make_mv.ipynb be run first to create the mv model.
./optimize_mv_bmi.ipynb -- jupyter notebook to optimize pumping rates using the linked-library version of MODFLOW 6 and the xmipy/bmi interface.
