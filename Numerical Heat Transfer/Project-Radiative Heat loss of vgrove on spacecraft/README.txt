################ Project ###############

##### CODE 1#####
 

   Calculations of the view Factors

######## CODE 2##############

	Radiosity calculations using  Gauss Siedel Method

################


	A V-groove has been machined into the side of a spacecraft orbiting the Earth. The V-
	groove must be there for structural support. However, you have concerns that this V-groove will cause
	the spacecraft to lose additional thermal energy to space. Completed the impact of the V-groove
	on the radiative thermal losses from the spacecraft.
	The V-groove is considered isothermal everywhere with a constant wall temperature of 300 K.
	The walls have a known emissivity of e (which varies in this problem) and the V-groove has an angle
	of phi. The surroundings of the V-groove are at a temperature of 0 K. The two sides of the V-groove are
	each 10 cm long. The V-groove is symmetric.

	Discretized the V-groove into equally-sized control surfaces to 10 control surfaces per side.
	Using your Monte Carlo code,  the view factor between a control surface on the left side of
	the V-groove with every control surface on the opposite side of the V-groove is detrermined. 
	Performed a ray convergence study to ensure  the correct number of rays.
	
	Code predicts the radiosity of each control surface along one wall (the V-groove is
	symmetrical, so both walls have the same radiosity profile). The radiosity balance equation given
	below, solved using Gauss-Seidel iteration. In this equation, the subscript i
	indicates the control surface for which you are determining the radiosity, and the subscript j indicates a
	separate control surface on the opposite surface.

		JiAi=e*Ai*sigma*Ti^4 + (1-e)(Σj=1)^N *Jj*Fj-1Aj( some reason my superscripts and subscript didnt work)
	
	The radiosity as a function of wall position (starting from the center of the V-groove and moving out towards the top of the V-groove) for emissivity values of 0.1, 	0.5 and 0.9 and for angles of 30°, 90° and 120° is determined.

	