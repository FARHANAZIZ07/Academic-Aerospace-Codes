function[ Po2_P0a,A_Ac,Norm_Massflow,To2,OPR,f] = EngineInletmatching(M,alpha,U_N,inletdata,Engine_data,ALTT,u_inlet,Altitudes,Ac,G,R)
%  ENGINEINLETMATCHING Brief summary of this function.
% 
% Detailed explanation of this function.
Tsls=518.69*5/9;
% M = [1.5:0.2:4];
Inlet_Massflow = zeros(1,length(M));
Norm_Massflow = zeros(1,length(M));
% alpha=[1:0.2:20];
for ii=1:1:length(alpha)
%for k=[1 2]
    for j = 1:length(Altitudes)
    Ta = interp1(ALTT.Alt,ALTT.Ta_values,Altitudes(j));
    for i = 1:length(M)
        [P,A_Ac] = Inletfunction(M(i),alpha(ii),inletdata,u_inlet);
        %% inlet
        A_As = (((G+1)/2)^-((G+1)/(G-1)/2)) / M(i) * (1 + M(i)^2 * (G-1)/2)^((G+1)/(G-1)/2);
        Inlet_Massflow(i) = ((sqrt(G/R))*(2/(G+1))^((G+1)/(2*(G-1))))*(1/P)*(1./(A_As))*A_Ac*Ac(u_inlet);
        %%Engine
        Toa = Ta*(1+(G-1)/2*M(i).^2);
        To2 = Toa;
        thetha = To2./Tsls;
        Norm_RPM = U_N./sqrt(thetha);
        Norm_Massflow(i) = interp1(Engine_data.Norm_RPM_values,Engine_data.Norm_Massflow_Values,Norm_RPM);
        OPR(i)=interp1(Engine_data.Norm_RPM_values,Engine_data.P_R_values,Norm_RPM);
        f(i)=interp1(Engine_data.Norm_RPM_values,Engine_data.f_values,Norm_RPM);
        
    if Inlet_Massflow(i)<Norm_Massflow(i)
               Po2_P0a(i)=P.*(Inlet_Massflow(i)./Norm_Massflow(i));
    else 
            Po2_P0a(i)=P;%critical& sub crititcal
    end
    end...
%end
end
end