function [Po2,A_Ac,Norm_Massflow,To2,To4,f] = PressureFunction(u_Mach,u_alpha,U_N,inletdata,Engine_data,ALTT,u_inlet,Altitudes,Ac,G,R)
%  PRESSUREFUNCTION Brief summary of this function.
% 
% Detailed explanation of this function.
[ Po2_P0a,A_Ac,Norm_Massflow,To2,OPR,f] = EngineInletmatching(u_Mach,u_alpha,U_N,inletdata,Engine_data,ALTT,u_inlet,Altitudes,Ac,G,R);
Pa=interp1(ALTT.Alt,ALTT.Pa_values,Altitudes);
Poa=Pa*(1+(G-1)/2*u_Mach^2)^(G/(G-1));
Po2= Po2_P0a*Poa;
Gh=1.33;
eta_c=0.95;
eta_b=0.99;
QR=4.83*10^8;
Cp_c=(G*R)/(G-1);
Cp_h=(Gh*R)/(Gh-1);
To3=To2.*(1+1./eta_c.*(OPR^((G-1)./G)-1));
To4=(Cp_c.*To3+eta_b.*QR.*f*10^-2)./(Cp_h.*(1+f.*10^-2));
end