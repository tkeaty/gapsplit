=================================================
GapSplit: universal sampling for metabolic models
=================================================

* GapSplit is a sampling algorithm designed to generate uniform, high-coverage sample points on any metabolic model
* regardless of convexity (i.e. logical/integer constraints).

Functions
---------

sample(fname, n_points, lower_bounds=None, upper_bounds=None, n_update=100, n_secondary=0)
   * Generate samples from a given input model.

   INPUT:
      * fname - str
            * String representing path to model file (see gurobipy.read() for acceptable file types).

      * n_points - int
            * Number of desired sample points.

      * lower_bounds - list/ndarray, optional
            * FVA minimums for model. Generated if not provided.

      * upper_bounds - list/ndarray, optional
            * FVA maximums for model. Generated if not provided.

      * n_update - int, optional
            * Refresh rate (in points) for console output of current model coverage and sample count.

      * n_secondary - int, optional
            * Number of additional gaps targeted for splitting.

   OUTPUT:
      * samples - ndarray
            * n_points by n_reactions array of sample points.


Dependencies
------------
 * gurobipy: 7.0 and up (requires download and license from gurobi.com - license provided free for academic users)
 * numpy: 1.14.5
