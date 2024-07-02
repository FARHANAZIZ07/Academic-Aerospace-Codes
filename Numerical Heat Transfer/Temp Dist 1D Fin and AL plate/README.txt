############
	Code to predict the temperature distribution
	in a 1D cylindrical fin exposed to convection and radiation on its outer surface. The fin has a base
	temperature of 400 K, a thermal conductivity of 401 W / m K, a convection heat transfer coefficient of
	10 W / m2 K, an ambient fluid temperature / surroundings temperature of 273 K, an emissivity of 0.5,
	and is 0.02 meters long and 0.003 meters in diameter. The governing equation is solved using both
	Gauss-Seidel iteration or matrix inversion. Hint: Used over-relaxation to force
	solution to converge. It is very easy for radiative boundary conditions to cause divergence!
	
	a.Provide a plot showing the steady state temperature distribution in the fin for an emissivity of 0
	i.e. no radiation). Note: A convergence study is performed  to discover the correct
	number of nodes to use. Use a convergence criteria of at least 0.0001.
	
	b. Provide a plot showing the steady state temperature distribution in the fin for an emissivity of
	0.5. Note:  A convergence study is performed to discover the correct number of nodes
	to use.
###########
 	Code to predict the temperature distribution as a function of x and y for a 2D slab of aluminum
	(k = 300 W / m K) that is 1 meter by 1 meter. The slab has a uniform temperature of T = 400 K on the
	bottom surface (y = 0) and 300 K on all of the other surfaces. Code  uses Gauss Seidel iteration
	and Matrix inversion.
	
	a. Provide a 2D plot showing the temperature distribution in your wall. Note: 
	performed a convergence study to discover the correct number of nodes to use. Use a convergence
	criteria of at least 0.0001. Note: Used the Seaborn Python package to produce easy 2D plots

#########################