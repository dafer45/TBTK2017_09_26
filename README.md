# Compilation  
First install TBTK following instructions on https://github.com/dafer45/TBTK/.
Make sure the installed version is v0.9.4. If TBTK has been compiled with CUDA
support this project can be compiled immediately by typing 'make'. If TBTK has
been compiled without CUDA support, the lines

```CPropertyExtractor propertyExtractor(  
	cSolver,  
	NUM_COEFFICIENTS,  
	true,  
	false,  
	true  
);```

in src/main.cpp first need to be modified to 

```CPropertyExtractor propertyExtractor(  
	cSolver,  
	NUM_COEFFICIENTS,  
	true,  
	false,  
	true  
);```

and the makefile will need to be modified in accordance with instructions
provided there. Once this is done, type 'make' to compile the project.

# Setup calculation
To setup the calculation for a certain set of parameters, open the file
"Parameters" and change the values. The first ten values controls the model,
the eleventh parameter controls whether the calculations are done along a 1D
cut or over the full surface, and the last parameters controls the
ChebyshevSolver.

# Run calculation
Type './build/a.out' to run a calculation. Note that the calculations can take
many hours to complete (in 2017) for lattice sizes comparable to those in the
corresponding article (201x201) even if a GPU is used. At the time of writing
the full 2D calculation is likely to be too computationally demanding to be
done on a CPU other than for comparatively small lattice sizes.

# Plot results
Two different plot scripts are provided for the two cases where cut1D is set
to true (1) and false (0). If a calculation was done with cut1D=1, type
'./plot.py TBTKResults.h5 1.57 0 0.005' to plot the results. Here The first
parameter is the file TBTKResults.h5 where the results are stored after the
calculation, the second and third are the polar and azimuthal angles of the
spin-polarization, respectively, and the fourth is a smoothing parameter
(sigma) used to perform Gaussian smoothing in the energy direction. The value
of these parameters can be played with to see the result for different
polarization axes and smoothing. The ploted results are found in the folder
'figures'. If the calculation instead was done with cut1D=0, replace
'./plot.py'with './plot2D.py'.
