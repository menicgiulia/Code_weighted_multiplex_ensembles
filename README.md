# Code_weighted_multiplex_ensembles
This repository contains:

1)'code_entropy_duplex.m' a Matlab function for the computation of the entropy value and the related Lagrangian multipliers for a canonical duplex ensemble with given expected multidegree and multistrength.

2)'main_single_instance.m' a Matlab script for the generation of a new duplex. It requires as INPUT, Z,T01, T10, T11, D01b, D10a, D11a, D11b, matrices hat are the results of 'code_entropy_duplex.m'. This script calls the function 'duplexsingleinstance.m'

3)'duplexsingleinstance.m',  see the script "main_single_instance.m" for explanations

This program is distributed by the authors in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. If you use this code you should cite:

[1] G. Menichetti, D. Remondini, G. Bianconi, "Correlations between weights and overlap in ensembles of weighted multiplex networks", Phys. Rev. E 90, 062817

http://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.062817

[2] G. Menichetti, D. Remondini, P. Panzarasa, R. J. Mondragón, G. Bianconi, "Weighted Multiplex Networks", PLoS ONE 9(6): e97857. doi: 10.1371/journal.pone.0097857

http://dx.doi.org/10.1371/journal.pone.0097857

Further details in the supplementary materials.
