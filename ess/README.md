## Enhanced Scatter Search

This repository contains a C implementation of derivation of a Scatter Search algorithm developed by Egea, J. A. et al. 2009, 2010. The first implementation was in matlab in a so called eSS package which is not supported anymore; in fact, the improved and parallelism matlab version is now included in the package called AMIGO. The source code for improved eSS implementation is available in form of M-packages. The C version presented here is based on the latest changes to the algorithm but not provide the parallelize algorithm yet. 

### Enhanced Scatter Search Algorithm 

Enhanced Scatter Search is a meta heuristic population based algorithm with similarities to evolutionary strategies. The algorithms follows the main procedure of original [Scatter Search](ss/README.md) algorithm; however, it deviates in refinement process in which algorithm tries to improve on good solutions over and over, in a process introduced as `goBeyond`.

### Prerequisites

The local search algorithms, Nelder-Mead and Levenberg-Marquardt, are provided using the GNU Scientific Library (GSL). 

### References

- Egea, J. a, Balsa-canto, E., Garcı, G., & Banga, J. R. (2009). Dynamic Optimization of Nonlinear Processes with an Enhanced Scatter Search Method, 4388–4401.
- Egea, J. a., Martí, R., & Banga, J. R. (2010). An evolutionary method for complex-process optimization. Computers and Operations Research, 37(2), 315–324. http://doi.org/10.1016/j.cor.2009.05.003
- Villaverde, A. F., Henriques, D., Smallbone, K., Bongard, S., Schmid, J., Cicin-Sain, D., … Banga, J. R. (2015). BioPreDyn-bench: a suite of benchmark problems for dynamic modelling in systems biology. BMC Systems Biology, 9(1), 8. http://doi.org/10.1186/s12918-015-0144-4