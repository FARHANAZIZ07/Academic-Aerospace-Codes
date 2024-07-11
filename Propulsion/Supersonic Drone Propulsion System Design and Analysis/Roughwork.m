clear all
clc
syms G R P A_Ac u_inlet A_As M Ac
A_As=((G+1)/2)^(-(G+1)/(2*(G-1)))*((1+((G-1)/2)*M^2)/M)^((G+1)/(2*(G-1)))
Inlet_Massflow_I1= ((sqrt(G/R))*(2/(G+1))^((G+1)/(2*(G-1))))*(1/P)*(1/(A_As))*A_Ac*Ac
A_AL= {[(G+1)/2]^-[(G+1)/(G-1)/2]} / M * [1 + M^2 * (G-1)/2]^[(G+1)/(G-1)/2]