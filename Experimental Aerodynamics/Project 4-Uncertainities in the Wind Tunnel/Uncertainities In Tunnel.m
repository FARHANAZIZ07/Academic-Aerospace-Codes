clc; close all; clc;

%% Load All Data %%

Data_Z_100g=readmatrix('force.txt');
Data_Z_500g=readmatrix('force1.txt');
Data_Z_1kg=readmatrix('force2.txt');
Data_Z_2kg=readmatrix('force3.txt');

Data_X_100g=readmatrix('force4.txt');
Data_X_500g=readmatrix('force5.txt');
Data_X_1kg=readmatrix('force6.txt');
Data_X_2kg=readmatrix('force7.txt');

Data_Y_100g_2in=readmatrix('force8.txt');
Data_Y_100g_1in=readmatrix('force9.txt');
Data_Y_200g_1in=readmatrix('force10.txt');
Data_Y_200g_2in=readmatrix('force11.txt');

Data_RC_10_Neg_u=readmatrix('force12.txt');
Data_RC_10_None_u=readmatrix('force13.txt');
Data_RC_10_Pos_u=readmatrix('force14.txt');

Data_RC_15_Neg_u=readmatrix('force15.txt');
Data_RC_15_None_u=readmatrix('force16.txt');
Data_RC_15_Pos_u=readmatrix('force17.txt');

Data_RC_10_Mount=readmatrix('force18.txt');
Data_RC_15_Mount=readmatrix('force19.txt');

Data_RC_10_Neg=Data_RC_10_Neg_u-Data_RC_10_Mount;
Data_RC_10_None=Data_RC_10_None_u-Data_RC_10_Mount;
Data_RC_10_Pos=Data_RC_10_Pos_u-Data_RC_10_Mount;

Data_RC_15_Neg=Data_RC_15_Neg_u-Data_RC_15_Mount;
Data_RC_15_None=Data_RC_15_None_u-Data_RC_15_Mount;
Data_RC_15_Pos=Data_RC_15_Pos_u-Data_RC_15_Mount;

%% Part 1: Calibration %%

% Average the Data, Break Into Components %

Fx_m_Z_100g=mean(Data_Z_100g(:,2));
Fy_m_Z_100g=mean(Data_Z_100g(:,3));
Fz_m_Z_100g=mean(Data_Z_100g(:,4));
Mx_m_Z_100g=mean(Data_Z_100g(:,5));
My_m_Z_100g=mean(Data_Z_100g(:,6));
Mz_m_Z_100g=mean(Data_Z_100g(:,7));

Fx_m_Z_500g=mean(Data_Z_500g(:,2));
Fy_m_Z_500g=mean(Data_Z_500g(:,3));
Fz_m_Z_500g=mean(Data_Z_500g(:,4));
Mx_m_Z_500g=mean(Data_Z_500g(:,5));
My_m_Z_500g=mean(Data_Z_500g(:,6));
Mz_m_Z_500g=mean(Data_Z_500g(:,7));

Fx_m_Z_1kg=mean(Data_Z_1kg(:,2));
Fy_m_Z_1kg=mean(Data_Z_1kg(:,3));
Fz_m_Z_1kg=mean(Data_Z_1kg(:,4));
Mx_m_Z_1kg=mean(Data_Z_1kg(:,5));
My_m_Z_1kg=mean(Data_Z_1kg(:,6));
Mz_m_Z_1kg=mean(Data_Z_1kg(:,7));

Fx_m_Z_2kg=mean(Data_Z_2kg(:,2));
Fy_m_Z_2kg=mean(Data_Z_2kg(:,3));
Fz_m_Z_2kg=mean(Data_Z_2kg(:,4));
Mx_m_Z_2kg=mean(Data_Z_2kg(:,5));
My_m_Z_2kg=mean(Data_Z_2kg(:,6));
Mz_m_Z_2kg=mean(Data_Z_2kg(:,7));

Fx_m_X_100g=mean(Data_X_100g(:,2));
Fy_m_X_100g=mean(Data_X_100g(:,3));
Fz_m_X_100g=mean(Data_X_100g(:,4));
Mx_m_X_100g=mean(Data_X_100g(:,5));
My_m_X_100g=mean(Data_X_100g(:,6));
Mz_m_X_100g=mean(Data_X_100g(:,7));

Fx_m_X_500g=mean(Data_X_500g(:,2));
Fy_m_X_500g=mean(Data_X_500g(:,3));
Fz_m_X_500g=mean(Data_X_500g(:,4));
Mx_m_X_500g=mean(Data_X_500g(:,5));
My_m_X_500g=mean(Data_X_500g(:,6));
Mz_m_X_500g=mean(Data_X_500g(:,7));

Fx_m_X_1kg=mean(Data_X_1kg(:,2));
Fy_m_X_1kg=mean(Data_X_1kg(:,3));
Fz_m_X_1kg=mean(Data_X_1kg(:,4));
Mx_m_X_1kg=mean(Data_X_1kg(:,5));
My_m_X_1kg=mean(Data_X_1kg(:,6));
Mz_m_X_1kg=mean(Data_X_1kg(:,7));

Fx_m_X_2kg=mean(Data_X_2kg(:,2));
Fy_m_X_2kg=mean(Data_X_2kg(:,3));
Fz_m_X_2kg=mean(Data_X_2kg(:,4));
Mx_m_X_2kg=mean(Data_X_2kg(:,5));
My_m_X_2kg=mean(Data_X_2kg(:,6));
Mz_m_X_2kg=mean(Data_X_2kg(:,7));

Fx_m_Y_100g_2in=mean(Data_Y_100g_2in(:,2));
Fy_m_Y_100g_2in=mean(Data_Y_100g_2in(:,3));
Fz_m_Y_100g_2in=mean(Data_Y_100g_2in(:,4));
Mx_m_Y_100g_2in=mean(Data_Y_100g_2in(:,5));
My_m_Y_100g_2in=mean(Data_Y_100g_2in(:,6));
Mz_m_Y_100g_2in=mean(Data_Y_100g_2in(:,7));

Fx_m_Y_100g_1in=mean(Data_Y_100g_1in(:,2));
Fy_m_Y_100g_1in=mean(Data_Y_100g_1in(:,3));
Fz_m_Y_100g_1in=mean(Data_Y_100g_1in(:,4));
Mx_m_Y_100g_1in=mean(Data_Y_100g_1in(:,5));
My_m_Y_100g_1in=mean(Data_Y_100g_1in(:,6));
Mz_m_Y_100g_1in=mean(Data_Y_100g_1in(:,7));

Fx_m_Y_200g_2in=mean(Data_Y_200g_2in(:,2));
Fy_m_Y_200g_2in=mean(Data_Y_200g_2in(:,3));
Fz_m_Y_200g_2in=mean(Data_Y_200g_2in(:,4));
Mx_m_Y_200g_2in=mean(Data_Y_200g_2in(:,5));
My_m_Y_200g_2in=mean(Data_Y_200g_2in(:,6));
Mz_m_Y_200g_2in=mean(Data_Y_200g_2in(:,7));

Fx_m_Y_200g_1in=mean(Data_Y_200g_1in(:,2));
Fy_m_Y_200g_1in=mean(Data_Y_200g_1in(:,3));
Fz_m_Y_200g_1in=mean(Data_Y_200g_1in(:,4));
Mx_m_Y_200g_1in=mean(Data_Y_200g_1in(:,5));
My_m_Y_200g_1in=mean(Data_Y_200g_1in(:,6));
Mz_m_Y_200g_1in=mean(Data_Y_200g_1in(:,7));
%%


% Applied Forces and Moments %

Fg_Container=0.1;
Fg_1=(0.1+Fg_Container)*9.81;
Fg_2=(0.5+Fg_Container)*9.81;
Fg_3=(1+Fg_Container)*9.81;
Fg_4=(2+Fg_Container)*9.81;
Fg_5=(0.2+Fg_Container)*9.81;
F_a=[Fg_1 Fg_2 Fg_3 Fg_4];

d_1=2*0.0254;
d_2=1.25*0.0254;

M_1=2*(Fg_1*d_1);
M_2=2*(Fg_1*d_2);
M_3=2*(Fg_5*d_1);
M_4=2*(Fg_5*d_2);
M_a=[M_2 M_4 M_1 M_3];

% Combine Data %

Fx_m_X=[Fx_m_X_100g Fx_m_X_500g Fx_m_X_1kg Fx_m_X_2kg];
Fz_m_X=[Fz_m_X_100g Fz_m_X_500g Fz_m_X_1kg Fz_m_X_2kg];
My_m_X=[My_m_X_100g My_m_X_500g My_m_X_1kg My_m_X_2kg];

Fx_m_Z=[Fx_m_Z_100g Fx_m_Z_500g Fx_m_Z_1kg Fx_m_Z_2kg];
Fz_m_Z=[Fz_m_Z_100g Fz_m_Z_500g Fz_m_Z_1kg Fz_m_Z_2kg];
My_m_Z=[My_m_Z_100g My_m_Z_500g My_m_Z_1kg My_m_Z_2kg];

Fx_m_Y=[Fx_m_Y_100g_1in Fx_m_Y_200g_1in Fx_m_Y_100g_2in Fx_m_Y_200g_2in];
Fz_m_Y=[Fz_m_Y_100g_1in Fz_m_Y_200g_1in Fz_m_Y_100g_2in Fz_m_Y_200g_2in];
My_m_Y=[My_m_Y_100g_1in My_m_Y_200g_1in My_m_Y_100g_2in My_m_Y_200g_2in];
%%


% Plots and Slopes %

p_Fx_X=polyfit(F_a,Fx_m_X,1);
Fx_Calc_X=polyval(p_Fx_X,F_a);
dFx_dFx=p_Fx_X(1,1);

p_Fz_X=polyfit(F_a,Fz_m_X,1);
Fz_Calc_X=polyval(p_Fz_X,F_a);
dFz_dFx=p_Fz_X(1,1);

p_My_X=polyfit(F_a,My_m_X,1);
My_Calc_X=polyval(p_My_X,F_a);
dMy_dFx=p_My_X(1,1);

p_Fx_Z=polyfit(-F_a,Fx_m_Z,1);
Fx_Calc_Z=polyval(p_Fx_Z,-F_a);
dFx_dFz=p_Fx_Z(1,1);

p_Fz_Z=polyfit(-F_a,Fz_m_Z,1);
Fz_Calc_Z=polyval(p_Fz_Z,-F_a);
dFz_dFz=p_Fz_Z(1,1);

p_My_Z=polyfit(-F_a,My_m_Z,1);
My_Calc_Z=polyval(p_My_Z,-F_a);
dMy_dFz=p_My_Z(1,1);

p_Fx_Y=polyfit(M_a,Fx_m_Y,1);
Fx_Calc_Y=polyval(p_Fx_Y,M_a);
dFx_dMy=p_Fx_Y(1,1);

p_Fz_Y=polyfit(M_a,Fz_m_Y,1);
Fz_Calc_Y=polyval(p_Fz_Y,M_a);
dFz_dMy=p_Fz_Y(1,1);

p_My_Y=polyfit(M_a,My_m_Y,1);
My_Calc_Y=polyval(p_My_Y,M_a);
dMy_dMy=p_My_Y(1,1);

figure (1)
plot(F_a,Fx_Calc_X,'k--',F_a,Fx_m_X,'k*',F_a,Fz_Calc_X,'r--',F_a,Fz_m_X,'r*',F_a,My_Calc_X,'b--',F_a,My_m_X,'b*')
xlabel('Applied Force in x-direction'),ylabel('Measured Force'),title('Applied vs Measured Forces - X-axis')
legend('Calculated Fx','Measured Fx','Calculated Fz','Measured Fz','Calculated My','Measured My')

figure (2)
plot(-F_a,Fx_Calc_Z,'k--',-F_a,Fx_m_Z,'k*',-F_a,Fz_Calc_Z,'r--',-F_a,Fz_m_Z,'r*',-F_a,My_Calc_Z,'b--',-F_a,My_m_Z,'b*')
xlabel('Applied Force in Z-direction'),ylabel('Measured Force'),title('Applied vs Measured Forces - Z-axis')
legend('Calculated Fx','Measured Fx','Calculated Fz','Measured Fz','Calculated My','Measured My')

figure (3)
plot(M_a,Fx_Calc_Y,'k--',M_a,Fx_m_Y,'k*',M_a,Fz_Calc_Y,'r--',M_a,Fz_m_Y,'r*',M_a,My_Calc_Y,'b--',M_a,My_m_Y,'b*')
xlabel('Applied Moment in Y-direction'),ylabel('Measured Force'),title('Applied vs Measured Forces - Y-axis')
legend('Calculated Fx','Measured Fx','Calculated Fz','Measured Fz','Calculated My','Measured My')

Slopes=[dFx_dFx dFz_dFx dMy_dFx;dFx_dFz dFz_dFz dMy_dFz;dFx_dMy dFz_dMy dMy_dMy];
Calibration=inv(Slopes);
%%


%% Part 2: RC Plane %%

% Average the Data, Break Into Components %

Fx_m_RC_10_Neg=mean(Data_RC_10_Neg(:,2));
Fy_m_RC_10_Neg=mean(Data_RC_10_Neg(:,3));
Fz_m_RC_10_Neg=mean(Data_RC_10_Neg(:,4));
Mx_m_RC_10_Neg=mean(Data_RC_10_Neg(:,5));
My_m_RC_10_Neg=mean(Data_RC_10_Neg(:,6));
Mz_m_RC_10_Neg=mean(Data_RC_10_Neg(:,7));

Fx_m_RC_10_None=mean(Data_RC_10_None(:,2));
Fy_m_RC_10_None=mean(Data_RC_10_None(:,3));
Fz_m_RC_10_None=mean(Data_RC_10_None(:,4));
Mx_m_RC_10_None=mean(Data_RC_10_None(:,5));
My_m_RC_10_None=mean(Data_RC_10_None(:,6));
Mz_m_RC_10_None=mean(Data_RC_10_None(:,7));

Fx_m_RC_10_Pos=mean(Data_RC_10_Pos(:,2));
Fy_m_RC_10_Pos=mean(Data_RC_10_Pos(:,3));
Fz_m_RC_10_Pos=mean(Data_RC_10_Pos(:,4));
Mx_m_RC_10_Pos=mean(Data_RC_10_Pos(:,5));
My_m_RC_10_Pos=mean(Data_RC_10_Pos(:,6));
Mz_m_RC_10_Pos=mean(Data_RC_10_Pos(:,7));

Fx_m_RC_15_Neg=mean(Data_RC_15_Neg(:,2));
Fy_m_RC_15_Neg=mean(Data_RC_15_Neg(:,3));
Fz_m_RC_15_Neg=mean(Data_RC_15_Neg(:,4));
Mx_m_RC_15_Neg=mean(Data_RC_15_Neg(:,5));
My_m_RC_15_Neg=mean(Data_RC_15_Neg(:,6));
Mz_m_RC_15_Neg=mean(Data_RC_15_Neg(:,7));

Fx_m_RC_15_None=mean(Data_RC_15_None(:,2));
Fy_m_RC_15_None=mean(Data_RC_15_None(:,3));
Fz_m_RC_15_None=mean(Data_RC_15_None(:,4));
Mx_m_RC_15_None=mean(Data_RC_15_None(:,5));
My_m_RC_15_None=mean(Data_RC_15_None(:,6));
Mz_m_RC_15_None=mean(Data_RC_15_None(:,7));

Fx_m_RC_15_Pos=mean(Data_RC_15_Pos(:,2));
Fy_m_RC_15_Pos=mean(Data_RC_15_Pos(:,3));
Fz_m_RC_15_Pos=mean(Data_RC_15_Pos(:,4));
Mx_m_RC_15_Pos=mean(Data_RC_15_Pos(:,5));
My_m_RC_15_Pos=mean(Data_RC_15_Pos(:,6));
Mz_m_RC_15_Pos=mean(Data_RC_15_Pos(:,7));

% Knowns %

p=1.225;
v_1=10;
v_2=15;
q_1=0.5*p*v_1^2;
q_2=0.5*p*v_2^2;

c_root=0.07387;
c_tip=0.04844;
b_root=0.11168;
b_tip=0.14168;
a_in=c_root*b_root;
a_out= b_tip*0.5*(c_root+c_tip);
S=2*(a_in+a_out);
%%

% Calculating Applied Forces vs Measured %

Fl_10_m=[Fz_m_RC_10_Neg Fz_m_RC_10_None Fz_m_RC_10_Pos];
Fl_15_m=[Fz_m_RC_15_Neg Fz_m_RC_15_None Fz_m_RC_15_Pos];

Fd_10_m=[Fx_m_RC_10_Neg Fx_m_RC_10_None Fx_m_RC_10_Pos];
Fd_15_m=[Fx_m_RC_15_Neg Fx_m_RC_15_None Fx_m_RC_15_Pos];

My_10_m=[My_m_RC_10_Neg My_m_RC_10_None My_m_RC_10_Pos];
My_15_m=[My_m_RC_15_Neg My_m_RC_15_None My_m_RC_15_Pos];

Measured_10=[Fd_10_m;Fl_10_m;My_10_m];
Measured_15=[Fd_15_m;Fl_15_m;My_15_m];

Applied_10=Calibration*Measured_10;
Applied_15=Calibration*Measured_15;

Fl_10_a=Applied_10(2,:);
Fl_15_a=Applied_15(2,:);
Fd_10_a=Applied_10(1,:);
Fd_15_a=Applied_15(1,:);
My_10_a=Applied_10(3,:);
My_15_a=Applied_15(3,:);

% Calculating Coefficients %

Cl_10_m=abs(Fl_10_m)/(q_1*S);
Cl_15_m=abs(Fl_15_m)/(q_2*S);
Cd_10_m=Fd_10_m/(q_1*S);
Cd_15_m=Fd_15_m/(q_2*S);
Cm_10_m=My_10_m/(q_1*S*c_root);
Cm_15_m=My_15_m/(q_2*S*c_root);

Cl_10=abs(Fl_10_a)/(q_1*S);
Cl_15=abs(Fl_15_a)/(q_2*S);
Cd_10=Fd_10_a/(q_1*S);
Cd_15=Fd_15_a/(q_2*S);
Cm_10=My_10_a/(q_1*S*c_root);
Cm_15=My_15_a/(q_2*S*c_root);

Position=[-1 0 1];
%%


% Plotting Coefficients %

figure (4)
plot(Position,Cl_10,'b',Position,Cd_10,'r',Position,Cl_10_m,'b*',Position,Cd_10_m,'r*')
xlabel('Position'),ylabel('Force Coefficient'),title('Force Coefficient vs Elevator Position at 10 m/s')
legend 'Coefficient of Lift-Calibrated' 'Coefficient of Drag-Calibrated' 'Coefficient of Lift-Measured' 'Coefficient of Drag-Measured'

figure (5)
plot(Position,Cl_15,'b',Position,Cd_15,'r',Position,Cl_15_m,'b*',Position,Cd_15_m,'r*')
xlabel('Position'),ylabel('Force Coefficient'),title('Force Coefficient vs Elevator Position at 15 m/s')
legend 'Coefficient of Lift-Calibrated' 'Coefficient of Drag-Calibrated' 'Coefficient of Lift-Measured' 'Coefficient of Drag-Measured'

figure (6)
plot(Position,Cm_10,'b',Position,Cm_15,'r',Position,Cm_10_m,'b*',Position,Cm_15_m,'r*')
xlabel('Position'),ylabel('Moment Coefficient'),title('Moment Coefficient vs Elevator Position')
legend 'Calibrated 10 m/s' 'Calibrated 15 m/s' 'Measured 10 m/s' 'Measured 15 m/s'

% Moving Moments to New Location %

x_mount=4*0.0254;
x_CG=3*0.0254;

My_disp_10=((x_CG-x_mount).*Fl_10_a)+My_10_a;
My_disp_15=((x_CG-x_mount).*Fl_15_a)+My_15_a;

Cm_disp_10=((x_CG-x_mount).*Cl_10)/c_root+Cm_10;
Cm_disp_15=((x_CG-x_mount).*Cl_15)/c_root+Cm_15;

figure (7)
plot(Position,My_10_a,'r',Position,My_disp_10,'b')
xlabel('Position'),ylabel('Moment (Nm)'),title('Moment at Mount vs Center of Gravity at 10 m/s')
legend 'Mount' 'CG'

figure (8)
plot(Position,My_15_a,'r',Position,My_disp_15,'b')
xlabel('Position'),ylabel('Moment (Nm)'),title('Moment at Mount vs Center of Gravity at 15 m/s')
legend 'Mount' 'CG'
%%


%% Calculating Uncertainty %%

n=2;
Wy_r=0.005;

df_dF_10=1/(q_1*S);
df_dF_15=1/(q_2*S);
df_dM_10=1/(q_1*S*c_root);
df_dM_15=1/(q_2*S*c_root);

W_Fd=1/(2*80);
W_Fl=1/(2*40);
W_My=10/(2*13333);
W_v=0.005;

df_dv_10_Fl=-((p*v_1*S)/(0.5*v_1^2*S)^2)*Fl_10_a;
df_dv_15_Fl=-((p*v_2*S)/(0.5*v_2^2*S)^2)*Fl_15_a;
df_dv_10_Fd=-((p*v_1*S)/(0.5*v_1^2*S)^2)*Fd_10_a;
df_dv_15_Fd=-((p*v_2*S)/(0.5*v_2^2*S)^2)*Fd_15_a;

W_total_10=sqrt((df_dF_10*W_Fd).^2+(df_dF_10*W_Fl).^2+(df_dv_10_Fl*W_v).^2+(df_dv_10_Fd*W_v).^2+(df_dM_10*W_My).^2+(df_dM_15*W_My).^2);
W_total_15=sqrt((df_dF_15*W_Fd).^2+(df_dF_15*W_Fl).^2+(df_dv_15_Fl*W_v).^2+(df_dv_15_Fd*W_v).^2);

W_Need_Fl_10=(df_dF_10).^(-1)*(Wy_r/sqrt(n));
W_Need_Fl_15=(df_dF_15).^(-1)*(Wy_r/sqrt(n));

W_Need_v_10=(df_dv_10_Fl).^(-1)*(Wy_r/sqrt(n));
W_Need_v_15=(df_dv_15_Fl).^(-1)*(Wy_r/sqrt(n));