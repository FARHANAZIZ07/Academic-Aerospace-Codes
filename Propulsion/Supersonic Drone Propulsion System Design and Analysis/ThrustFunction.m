function [Po2,A_Ac,Norm_Massflow,To2] = PressureFunction(u_Mach,u_alpha,U_N,inletdata,Engine_data,ALTT,u_inlet,Altitudes,Ac,G,R)
%  PRESSUREFUNCTION Brief summary of this function.
% 
% Detailed explanation of this function.
[ Po2_P0a,A_Ac,Norm_Massflow,To2] = EngineInletmatching(u_Mach,u_alpha,U_N,inletdata,Engine_data,ALTT,u_inlet,Altitudes,Ac,G,R)
Pa=interp1(ALTT.Alt,ALTT.Pa_values,Altitudes);
Poa=Pa*(1+(G-1)/2*u_Mach^2)^(G/(G-1));
Po2= Po2_P0a*Poa;
end