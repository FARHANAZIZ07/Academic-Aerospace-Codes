clc; close all; clc;
load('Pressure_Calibrated.txt')
load('pressure_no_wind.txt')
load('pressure_3.txt')
load('wing_1_forces.mat')
load('wing_2_forces.mat')
load('wing_3_forces.mat')
alphas = unique(wing_1.Alpha);
rho = 1.225;
Re=240668.5083;
%% Wing 1 Uncorrected
% Inputs
velocity_1 = 14;
s_1 = 0.129032;
% Calculate dynamic pressure
q_1 = 0.5*rho*velocity_1^2;
% Allocate arrays
lifts_1 = zeros(length(alphas),1);
drags_1 = zeros(length(alphas),1);
cls_1 = zeros(length(alphas),1);
cds_1 = zeros(length(alphas),1);
% Calculate lift and drag from Fx and Fy
wing_1.lift = wing_1.Fy.*cosd(wing_1.Alpha)-wing_1.Fx.*sind(wing_1.Alpha);
wing_1.drag =  wing_1.Fy.*sind(wing_1.Alpha)+wing_1.Fx.*cosd(wing_1.Alpha);
% Calculate time-averaged lift, drag, cl, cd
for i = 1:length(alphas)
   alpha = alphas(i);
   tmp = wing_1(wing_1.Alpha==alpha,:);
   tmp_avg = mean(tmp);
   lifts_1(i) = tmp_avg.lift;
   drags_1(i) = tmp_avg.drag;
   cls_1(i) = lifts_1(i)/(q_1*s_1);
   cds_1(i) = drags_1(i)/(q_1*s_1);
end
%% Wing 2 Uncorrected
% Inputs
velocity_2 = 20;
s_2 = 0.06322568;
% Calculate dynamic pressure
q_2 = 0.5*rho*velocity_2^2;
% Allocate arrays
lifts_2 = zeros(length(alphas),1);
drags_2 = zeros(length(alphas),1);
cls_2 = zeros(length(alphas),1);
cds_2 = zeros(length(alphas),1);
% Calculate lift and drag from Fx and Fy
wing_2.lift = wing_2.Fy.*cosd(wing_2.Alpha)-wing_2.Fx.*sind(wing_2.Alpha);
wing_2.drag =  wing_2.Fy.*sind(wing_2.Alpha)+wing_2.Fx.*cosd(wing_2.Alpha);
% Calculate time-averaged lift, drag, cl, cd
for i = 1:length(alphas)
   alpha = alphas(i);
   tmp = wing_2(wing_2.Alpha==alpha,:);
   tmp_avg = mean(tmp);
   lifts_2(i) = tmp_avg.lift;
   drags_2(i) = tmp_avg.drag;
   cls_2(i) = lifts_2(i)/(q_2*s_2);
   cds_2(i) = drags_2(i)/(q_2*s_2);
end
%% Wing 3 Uncorrected
% Inputs
velocity_3 = 35;
s_3 = 0.02064512;
% Calculate dynamic pressure
q_3 = 0.5*rho*velocity_3^2;
% Allocate arrays
lifts_3 = zeros(length(alphas),1);
drags_3 = zeros(length(alphas),1);
cls_3 = zeros(length(alphas),1);
cds_3 = zeros(length(alphas),1);
% Calculate lift and drag from Fx and Fy
wing_3.lift = wing_3.Fy.*cosd(wing_3.Alpha)-wing_3.Fx.*sind(wing_3.Alpha);
wing_3.drag =  wing_3.Fy.*sind(wing_3.Alpha)+wing_3.Fx.*cosd(wing_3.Alpha);
% Calculate time-averaged lift, drag, cl, cd
for i = 1:length(alphas)
   alpha = alphas(i);
   tmp = wing_3(wing_3.Alpha==alpha,:);
   tmp_avg = mean(tmp);
   lifts_3(i) = tmp_avg.lift;
   drags_3(i) = tmp_avg.drag;
   cls_3(i) = lifts_3(i)/(q_3*s_3);
   cds_3(i) = drags_3(i)/(q_3*s_3);
end
%% Pressure Calcs
dx=1.5*0.0254; % Distance between sensors
P0_volts=mean(pressure_no_wind);
P_volts=mean(pressure_3);
P_volts_c=P_volts-P0_volts;
P_volts_c=P_volts_c(1,2:end-1);
Slopes=[2303.7, 2491.9, 2482.8, 2487.9, 2479.6, 2468.8, 5409.8];
P_Pascals=[];
for i=1:length(P_volts_c)
   P_Pascals(i)=(P_volts_c(i)*Slopes(i)/1000*6894.76)+98651.8;
end
P=P_Pascals;
x=[0 dx 2*dx 3*dx 4*dx 5*dx 6*dx];
p=polyfit(x,P_Pascals,1);
dp=p(1,1);
P_calc=polyval(p,x);
figure
plot(x,P_Pascals,'b',x,P_calc,'b--')
xlabel('X-Position (m)'),ylabel('Pressure (Pa)'),legend 'Data' 'Linear Regression'
%% Plots
figure
plot(alphas,drags_1,alphas,drags_2,alphas,drags_3)
xlabel('alpha'), ylabel('Drag'), title('Force of Drag vs Alpha')
legend('Wing 1','Wing 2','Wing 3')
figure
plot(alphas,lifts_1,alphas,lifts_2,alphas,lifts_3)
xlabel('alpha'), ylabel('Lift'), title('Force of Lift vs Alpha')
legend('Wing 1','Wing 2','Wing 3')
figure
plot(alphas,(lifts_1./drags_1),alphas,(lifts_2./drags_2),alphas,(lifts_3./drags_3))
xlabel('alpha'), ylabel('L/D'), title('L/D vs Alpha')
legend('Wing 1','Wing 2','Wing 3')
%% Factors in Blockage Calculations
c_1=10*0.0254;
c_2=7*0.0254;
c_3=4*0.0254;
t_1_0=0.12*c_1;
t_2_0=0.12*c_2;
t_3_0=0.12*c_3;
t_1=(0.88*c_1*sind(alphas))+t_1_0;
t_2=(0.88*c_2*sind(alphas))+t_2_0;
t_3=(0.88*c_3*sind(alphas))+t_3_0;
a_1=20*0.0254;
a_2=14*0.0254;
a_3=8*0.0254;
h=30*0.0254; % Tunnel Height
A_1=t_1*a_1;
A_2=t_2*a_2;
A_3=t_3*a_3;
A_tunnel=30*30*(0.0254^2);
C=h*h; % Tunnel Cross Sectional Area
BR_1=A_1/A_tunnel;
BR_2=A_2/A_tunnel;
BR_3=A_3/A_tunnel;
tcRatio=0.12;
ctRatio=tcRatio^(-1);
G=0.225; % From graph of Gamma vs t/c
BS_factor2=1.1; % From graph of Lambda vs t/c
Sigma_1=((pi()^2)/48)*((c_1/h)^2);
Sigma_2=((pi()^2)/48)*((c_2/h)^2);
Sigma_3=((pi()^2)/48)*((c_3/h)^2);
BS_factor3=4.25; % From graph of Lambda vs c/t
Volume_1=c_1*t_1_0*a_1;
Volume_2=c_2*t_2_0*a_2;
Volume_3=c_3*t_3_0*a_3;
K1_Thom=0.71;
K=0.9;
K1_Herriot=1;
K3_Herriot=0.915;
T1_Herriot=0.82;
m_1=2.9;
m_2=3.1;
m_3=3.3;
%% Calculating Blockages
% Horizontal Buoyancy 2D
dD_B_Glauert_1=0.5*pi()*BS_factor2*(t_1.^2)*dp*a_1;
dD_B_Glauert_2=0.5*pi()*BS_factor2*(t_2.^2)*dp*a_2;
dD_B_Glauert_3=0.5*pi()*BS_factor2*(t_3.^2)*dp*a_3;
dD_B_Allen_1=((6*(h^2))/(pi()))*G*Sigma_1*dp*a_1;
dD_B_Allen_2=((6*(h^2))/(pi()))*G*Sigma_2*dp*a_2;
dD_B_Allen_3=((6*(h^2))/(pi()))*G*Sigma_3*dp*a_3;
% Horizontal Buoyancy 3D
dD_B_Glauert3D_1=(-pi()/4)*BS_factor3*(t_1.^2)*dp;
dD_B_Glauert3D_2=(-pi()/4)*BS_factor3*(t_2.^2)*dp;
dD_B_Glauert3D_3=(-pi()/4)*BS_factor3*(t_3.^2)*dp;
% Solid Blockage 2D
e_sb_Pope_1=((pi()^2)/3)*((a_1^2)/(h^2));
e_sb_Pope_2=((pi()^2)/3)*((a_2^2)/(h^2));
e_sb_Pope_3=((pi()^2)/3)*((a_3^2)/(h^2));
e_sb_Glauert_1=((pi()^2)/3)*(BS_factor2/4)*((t_1.^2)/(h^2));
e_sb_Glauert_2=((pi()^2)/3)*(BS_factor2/4)*((t_2.^2)/(h^2));
e_sb_Glauert_3=((pi()^2)/3)*(BS_factor2/4)*((t_3.^2)/(h^2));
e_sb_Allen_1=G*Sigma_1; % For Airfoils
e_sb_Allen_2=G*Sigma_2; % For Airfoils
e_sb_Allen_3=G*Sigma_3; % For Airfoils
e_sb_Thom_1=K1_Thom*(Volume_1/(C^(3/2)));
e_sb_Thom_2=K1_Thom*(Volume_2/(C^(3/2)));
e_sb_Thom_3=K1_Thom*(Volume_3/(C^(3/2)));
% Solid Blockage 3D
e_sb_Herriot_wing_1=(K1_Herriot*T1_Herriot)*(Volume_1/(C^(3/2))); % For Wings
e_sb_Herriot_wing_2=(K1_Herriot*T1_Herriot)*(Volume_2/(C^(3/2))); % For Wings
e_sb_Herriot_wing_3=(K1_Herriot*T1_Herriot)*(Volume_3/(C^(3/2))); % For Wings
e_sb_Herriot_G_1=(K3_Herriot*T1_Herriot)*(Volume_1/(C^(3/2)));
e_sb_Herriot_G_2=(K3_Herriot*T1_Herriot)*(Volume_2/(C^(3/2)));
e_sb_Herriot_G_3=(K3_Herriot*T1_Herriot)*(Volume_3/(C^(3/2)));
e_sb_Thom3D_1=K*(Volume_1/(C^(3/2)));
e_sb_Thom3D_2=K*(Volume_2/(C^(3/2)));
e_sb_Thom3D_3=K*(Volume_3/(C^(3/2)));
% Wake Blockage 2D
e_wb_Thom_1=(c_1/(4*h))*cds_1;
e_wb_Thom_2=(c_1/(4*h))*cds_2;
e_wb_Thom_3=(c_1/(4*h))*cds_3;
e_wb_Maskell_1=(c_1/(2*h))*cds_1;
e_wb_Maskell_2=(c_2/(2*h))*cds_2;
e_wb_Maskell_3=(c_3/(2*h))*cds_3;
e_wb_Allen_1=G*Sigma_1; % For Airfoils
e_wb_Allen_2=G*Sigma_2; % For Airfoils
e_wb_Allen_3=G*Sigma_3; % For Airfoils
% Wake Blockage 3D
e_wb_Maskell3D_1=1/(1-m_1);
e_wb_Maskell3D_2=1/(1-m_2);
e_wb_Maskell3D_3=1/(1-m_3);
e_wb_highv_1=(1/4)*(A_1./C).*cds_1;
e_wb_highv_2=(1/4)*(A_2./C).*cds_2;
e_wb_highv_3=(1/4)*(A_3./C).*cds_3;
% Streamline Curvature Blockage
dCl_sc_1=Sigma_1*cls_1;
dCl_sc_2=Sigma_2*cls_2;
dCl_sc_3=Sigma_3*cls_3;
% dCl_sc_1=0;
% dCl_sc_2=0;
% dCl_sc_3=0;
%% Plotting Blockage Factors
% Horizontal Buoyancy
figure
plot(alphas,dD_B_Glauert_1,alphas,dD_B_Glauert3D_1)
yline(dD_B_Allen_1)
xlabel('alpha'), ylabel('Change in Drag Force'), title('Horizontal Buoyancy Blockage Calculations - Wing 1')
legend('Glauert 2D','Glauert3D','Allen 2D')
figure
plot(alphas,dD_B_Glauert_2,alphas,dD_B_Glauert3D_2)
yline(dD_B_Allen_2)
xlabel('alpha'), ylabel('Change in Drag Force'), title('Horizontal Buoyancy Blockage Calculations - Wing 2')
legend('Glauert 2D','Glauert3D','Allen 2D')
figure
plot(alphas,dD_B_Glauert_2,alphas,dD_B_Glauert3D_3)
yline(dD_B_Allen_3)
xlabel('alpha'), ylabel('Change in Drag Force'), title('Horizontal Buoyancy Blockage Calculations - Wing 3')
legend('Glauert 2D','Glauert3D','Allen 2D')
% Solid Blockage
figure
plot(alphas,e_sb_Glauert_1)
yyaxis left
yline(e_sb_Allen_1,'g'), yline(e_sb_Thom_1,'r'),yline(e_sb_Herriot_wing_1,'m--'), yline(e_sb_Herriot_G_1,'m'), yline(e_sb_Thom3D_1,'r--')
yyaxis right
yline(e_sb_Pope_1)
xlabel('alpha'), ylabel('Solid Blockage Component'), title('Solid Blockage Calculations - Wing 1')
legend('Glauert 2D','Allen 2D','Thom 2D','Herriot Wing 3D','Herriot Generic 3D','Thom 3D','Pope 2D')
figure
plot(alphas,e_sb_Glauert_2)
yyaxis left
yline(e_sb_Allen_2,'g'), yline(e_sb_Thom_2,'r'),yline(e_sb_Herriot_wing_2,'m--'), yline(e_sb_Herriot_G_2,'m'), yline(e_sb_Thom3D_2,'r--')
yyaxis right
yline(e_sb_Pope_2)
xlabel('alpha'), ylabel('Solid Blockage Component'), title('Solid Blockage Calculations - Wing 2')
legend('Glauert 2D','Allen 2D','Thom 2D','Herriot Wing 3D','Herriot Generic 3D','Thom 3D','Pope 2D')
figure
plot(alphas,e_sb_Glauert_3)
yyaxis left
yline(e_sb_Allen_3,'g'), yline(e_sb_Thom_3,'r'),yline(e_sb_Herriot_wing_3,'m--'), yline(e_sb_Herriot_G_3,'m'), yline(e_sb_Thom3D_3,'r--')
yyaxis right
yline(e_sb_Pope_3)
xlabel('alpha'), ylabel('Solid Blockage Component'), title('Solid Blockage Calculations - Wing 3')
legend('Glauert 2D','Allen 2D','Thom 2D','Herriot Wing 3D','Herriot Generic 3D','Thom 3D','Pope 2D')
% Wake Blockage
figure
plot(alphas,e_wb_Thom_1,alphas,e_wb_Maskell_1,alphas,e_wb_highv_1)
yyaxis left
yline(e_wb_Allen_1)
yyaxis right
yline(e_wb_Maskell3D_1)
xlabel('alpha'), ylabel('Wake Blockage Component'), title('Wake Blockage Calculations - Wing 1')
legend('Thom 2D','Maskell 2D','High Velocity 3D','Allen 2D','Maskell 3D')
figure
plot(alphas,e_wb_Thom_2,alphas,e_wb_Maskell_2,alphas,e_wb_highv_2)
yyaxis left
yline(e_wb_Allen_2)
yyaxis right
yline(e_wb_Maskell3D_2)
xlabel('alpha'), ylabel('Wake Blockage Component'), title('Wake Blockage Calculations - Wing 2')
legend('Thom 2D','Maskell 2D','High Velocity 3D','Allen 2D','Maskell 3D')
figure
plot(alphas,e_wb_Thom_3,alphas,e_wb_Maskell_3,alphas,e_wb_highv_3)
yyaxis left
yline(e_wb_Allen_3)
yyaxis right
yline(e_wb_Maskell3D_3)
xlabel('alpha'), ylabel('Wake Blockage Component'), title('Wake Blockage Calculations - Wing 3')
legend('Thom 2D','Maskell 2D','High Velocity 3D','Allen 2D','Maskell 3D')
% Streamline Curvature
figure
plot(alphas,dCl_sc_1,alphas,dCl_sc_2,alphas,dCl_sc_3)
xlabel('alpha'), ylabel('Change in Lift Component'), title('Streamline Curvature Calculations')
legend('Wing 1','Wing 2','Wing 3')
%% Choose Blockage Values
dD_B_1=dD_B_Glauert_1;
dD_B_2=dD_B_Glauert_2;
dD_B_3=dD_B_Glauert_3;
% dD_B_1=0;
% dD_B_2=0;
% dD_B_3=0;
e_sb_1=e_sb_Glauert_1;
e_sb_2=e_sb_Glauert_2;
e_sb_3=e_sb_Glauert_3;
% e_sb_1=0;
% e_sb_2=0;
% e_sb_3=0;
e_wb_1=e_wb_Thom_1;
e_wb_2=e_wb_Thom_2;
e_wb_3=e_wb_Thom_3;
% e_wb_1=0;
% e_wb_2=0;
% e_wb_3=0;
e_t_1=e_sb_1+e_wb_1;
e_t_2=e_sb_2+e_wb_2;
e_t_3=e_sb_3+e_wb_3;
%% Apply Blockages to Forces and Coefficients
drag_c_1=drags_1-dD_B_1;
drag_c_2=drags_2-dD_B_2;
drag_c_3=drags_3-dD_B_3;
Cd_u_1= drag_c_1./(q_1*s_1);
Cd_u_2= drag_c_2./(q_2*s_2);
Cd_u_3= drag_c_3./(q_3*s_3);
Cd_0_c_1=Cd_u_1(1,1)*(1-(3*e_sb_1)-(2*e_wb_1));
Cd_0_c_2=Cd_u_2(1,1)*(1-(3*e_sb_2)-(2*e_wb_2));
Cd_0_c_3=Cd_u_3(1,1)*(1-(3*e_sb_3)-(2*e_wb_3));
Cl_u_1=cls_1.*(1-Sigma_1-(2.*e_t_1));
Cl_u_2=cls_2.*(1-Sigma_2-(2.*e_t_2));
Cl_u_3=cls_3.*(1-Sigma_3-(2.*e_t_3));
Cl_c_1=Cl_u_1-dCl_sc_1;
Cl_c_2=Cl_u_2-dCl_sc_2;
Cl_c_3=Cl_u_3-dCl_sc_3;
Cd_i_1=(Cl_c_1.^2)./(pi()*rho*2);
Cd_i_2=(Cl_c_2.^2)./(pi()*rho*2);
Cd_i_3=(Cl_c_3.^2)./(pi()*rho*2);
Cd_c_1=Cd_0_c_1+Cd_i_1;
Cd_c_2=Cd_0_c_2+Cd_i_2;
Cd_c_3=Cd_0_c_3+Cd_i_3;
velocity_c_1=velocity_1*(1+e_t_1);
velocity_c_2=velocity_2*(1+e_t_2);
velocity_c_3=velocity_3*(1+e_t_3);
q_c_1=q_1*(1+(2*e_t_1));
q_c_2=q_2*(1+(2*e_t_2));
q_c_3=q_3*(1+(2*e_t_3));
Re_c_1=Re*(1+e_t_1);
Re_c_2=Re*(1+e_t_2);
Re_c_3=Re*(1+e_t_3);
%% Plot Final Values
figure
semilogx(BR_1,Cd_c_1,'r--',BR_1,cds_1,'r',BR_2,Cd_c_2,'b--',BR_2,cds_2,'b',BR_3,Cd_c_3,'k--',BR_3,cds_3,'k')
xlabel('Body Ratio'), ylabel('Coefficient of Drag'), title('Corrected vs Uncorrected Cd')
legend('Corrected Cd - Wing 1','Uncorrected Cd - Wing 1','Corrected Cd - Wing 2','Uncorrected Cd - Wing 2','Corrected Cd - Wing 3','Uncorrected Cd - Wing 3')
figure
semilogx(BR_1,Cl_c_1,'r--',BR_1,cls_1,'r',BR_2,Cl_c_2,'b--',BR_2,cls_2,'b',BR_3,Cl_c_3,'k--',BR_3,cls_3,'k')
xlabel('Body Ratio'), ylabel('Coefficient of Lift'), title('Corrected vs Uncorrected Cl')
legend('Corrected Cl - Wing 1','Uncorrected Cl - Wing 1','Corrected Cl - Wing 2','Uncorrected Cl - Wing 2','Corrected Cl - Wing 3','Uncorrected Cl - Wing 3')
figure
semilogx(BR_1,velocity_c_1,'r--',BR_2,velocity_c_2,'b--',BR_3,velocity_c_3,'k--')
yline(velocity_1,'r'),yline(velocity_2,'b'),yline(velocity_3,'k')
xlabel('Body Ratio'), ylabel('Velocity (m/s)'), title('Corrected vs Uncorrected Velocity')
legend('Corrected Velocity - Wing 1','Corrected Velocity - Wing 2','Corrected Velocity - Wing 3','Uncorrected Velocity - Wing 1','Uncorrected Velocity - Wing 2','Uncorrected Velocity - Wing 3')
figure
semilogx(BR_1,q_c_1,'r--',BR_2,q_c_2,'b--',BR_3,q_c_3,'k--')
yline(q_1,'r'),yline(q_2,'b'),yline(q_3,'k')
xlabel('Body Ratio'), ylabel('Dynamic Pressure'), title('Corrected vs Uncorrected Dynamic Pressure')
legend('Corrected Dynamic Pressure - Wing 1','Corrected Dynamic Pressure - Wing 2','Corrected Dynamic Pressure - Wing 3','Uncorrected Dynamic Pressure - Wing 1','Uncorrected Dynamic Pressure - Wing 2','Uncorrected Dynamic Pressure - Wing 3')
figure
semilogx(BR_1,Re_c_1,'r--',BR_2,Re_c_2,'b--',BR_3,Re_c_3,'k--')
yline(Re,'m')
xlabel('Body Ratio'), ylabel('Reynolds Number'), title('Corrected vs Uncorrected Reynolds Number')
legend('Corrected - Wing 1','Corrected - Wing 2','Corrected - Wing 3','Uncorrected - All Wings')
figure
plot(alphas,Cd_c_1,'r--',alphas,cds_1,'r',alphas,Cd_c_2,'b--',alphas,cds_2,'b',alphas,Cd_c_3,'k--',alphas,cds_3,'k')
xlabel('Alpha'), ylabel('Coefficient of Drag'), title('Corrected vs Uncorrected Cd')
legend('Corrected Cd - Wing 1','Uncorrected Cd - Wing 1','Corrected Cd - Wing 2','Uncorrected Cd - Wing 2','Corrected Cd - Wing 3','Uncorrected Cd - Wing 3')
%% Percent Contribution
Cd_hb_1=[0.0219;0.0190;0.0198;0.0253;0.0352;0.0472;0.0602;0.0744;0.0875;0.0675;0.0546];
Cd_hb_2=[0.0175;0.0155;0.0164;0.0207;0.0309;0.0421;0.0551;0.0691;0.0837;0.0643;0.0575];
Cd_hb_3=[0.0300;0.0252;0.0262;0.0331;0.0552;0.0735;0.1044;0.1441;0.1786;0.1460;0.1047];
Cd_sb_1=[0.0590;0.0561;0.0570;0.0626;0.0728;0.0853;0.0990;0.1143;0.1288;0.1077;0.0943];
Cd_sb_2=[0.0302;0.0282;0.0292;0.0335;0.0438;0.0553;0.0687;0.0831;0.0983;0.0785;0.0716];
Cd_sb_3=[0.0324;0.0276;0.0286;0.0356;0.0578;0.0762;0.1073;0.1474;0.1824;0.1498;0.1081];
Cd_wb_1=[0.0589;0.0559;0.0567;0.0621;0.0721;0.0843;0.0976;0.1125;0.1265;0.1077;0.0942];
Cd_wb_2=[0.0302;0.0282;0.0292;0.0335;0.0439;0.0555;0.0691;0.0839;0.0997;0.0819;0.0751];
Cd_wb_3=[0.0326;0.0277;0.0287;0.0359;0.0587;0.0780;0.1109;0.1540;0.1929;0.1641;0.1217];
Cd_sc_1=[0.0589;0.0557;0.0565;0.0620;0.0722;0.0845;0.0979;0.1125;0.1260;0.1041;0.0901];
Cd_sc_2=[0.0302;0.0281;0.0291;0.0334;0.0438;0.0552;0.0685;0.0827;0.0975;0.0775;0.0703];
Cd_sc_3=[0.0324;0.0276;0.0286;0.0356;0.0578;0.0762;0.1073;0.1473;0.1820;0.1492;0.1075];
percent_hb_1=(Cd_c_1-Cd_hb_1)./Cd_c_1;
percent_hb_2=(Cd_c_2-Cd_hb_2)./Cd_c_2;
percent_hb_3=(Cd_c_3-Cd_hb_3)./Cd_c_3;
percent_hb=[mean(percent_hb_1) mean(percent_hb_2) mean(percent_hb_3)];
percent_sb_1=(Cd_c_1-Cd_sb_1)./Cd_c_1;
percent_sb_2=(Cd_c_2-Cd_sb_2)./Cd_c_2;
percent_sb_3=(Cd_c_3-Cd_sb_3)./Cd_c_3;
percent_sb=[mean(percent_sb_1) mean(percent_sb_2) mean(percent_sb_3)];
percent_wb_1=(Cd_c_1-Cd_wb_1)./Cd_c_1;
percent_wb_2=(Cd_c_2-Cd_wb_2)./Cd_c_2;
percent_wb_3=(Cd_c_3-Cd_wb_3)./Cd_c_3;
percent_wb=[mean(percent_wb_1) mean(percent_wb_2) mean(percent_wb_3)];
percent_sc_1=(Cd_c_1-Cd_sc_1)./Cd_c_1;
percent_sc_2=(Cd_c_2-Cd_sc_2)./Cd_c_2;
percent_sc_3=(Cd_c_3-Cd_sc_3)./Cd_c_3;
percent_sc=[mean(percent_sc_1) mean(percent_sc_2) mean(percent_sc_3)];