# LM-RBF-CFD

Hydrodynamics code using radial basis functions and Levenberg-Marquardt method

1. A file named `problem.data` must exist in the run directory, of the format  
   &ProblemData  
   x_dim =                 ! Width in CU  
   y_dim =                 ! Height in CU  
   res =                   ! points/CU  
   tol =                   ! Tolerance for LM solver  
   BC =                    ! Boundary conditions to use. [1] = Channel flow, [2] = Couette flow  
   Re =                    ! Reynolds number of fluid  
   /                       ! This line must be included  

2. A file named `constants.data` must exist in the run directory, of the format  
   &ConstantsData  
   uleft = 1               ! x velocity of fluid at left boundary  
   vleft = 0               ! y velocity of fluid at left boundary  
   pleft = 0               ! pressure of fluid at left boundary  
   uwall = 0               ! velocity of wall  
   /                       ! This line must be included  

3. A folder named `out` must exist in the run directory.

4. A file named `alpha.dat`, containing an initial guess of the RBF decomposition coefficients, *may* exist in the run directory. If it exists, it must be of the format  
   4d-1  
   5d1  
   6d-2  
   3d2  
   ... etc, where ui = alpha(:nodes), vi = alpha(nodes:nodes\*2), and pi = alpha(nodes\*2:nodes\*3)  
   If it does not exist, the default guess of all 1s will be used for the decomposition coefficients.
