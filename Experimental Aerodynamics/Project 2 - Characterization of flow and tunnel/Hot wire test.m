clear; clc;
%% Calibration %%
load('Calibration.mat')
load('Test_Data.mat')
V_Cal=[0.05 1.1 2.34 4.91 7.46 10.07 12.78 15.11 17.6 20.12 22.62 25.15];
%%

% Average Voltage for Each Speed %
E_0rpm=mean(Cal_0rpm);
E_0rpm=table2array(E_0rpm);
E_25rpm=mean(Cal_25rpm);
E_25rpm=table2array(E_25rpm);
E_50rpm=mean(Cal_50rpm);
E_50rpm=table2array(E_50rpm);
E_100rpm=mean(Cal_100rpm);
E_100rpm=table2array(E_100rpm);
E_150rpm=mean(Cal_150rpm);
E_150rpm=table2array(E_150rpm);
E_200rpm=mean(Cal_200rpm);
E_200rpm=table2array(E_200rpm);
E_250rpm=mean(Cal_250rpm);
E_250rpm=table2array(E_250rpm);
E_300rpm=mean(Cal_300rpm);
E_300rpm=table2array(E_300rpm);
E_350rpm=mean(Cal_350rpm);
E_350rpm=table2array(E_350rpm);
E_400rpm=mean(Cal_400rpm);
E_400rpm=table2array(E_400rpm);
E_450rpm=mean(Cal_450rpm);
E_450rpm=table2array(E_450rpm);
E_500rpm=mean(Cal_500rpm);
E_500rpm=table2array(E_500rpm);
E_Cal=[E_0rpm E_25rpm E_50rpm E_100rpm E_150rpm E_200rpm E_250rpm E_300rpm E_350rpm E_400rpm E_450rpm E_500rpm];
p=polyfit(E_Cal,V_Cal,4);
E_fit=linspace(1.3,3,12);
V_fit=polyval(p,E_fit);
%%

% Plot Calibration %
figure (1)
plot(E_Cal,V_Cal,'b',E_fit,V_fit,'o')
xlabel('Voltage (V)');
ylabel('Velocity (m/s)');
title 
%%
'Velocity vs Voltage for Calibration of Hot Wire'

%% Flow Turbulance Intensity  %%
t=linspace(0,5,100000);
Test_TI_10ms=table2array(Test_TI_10ms);
Test_TI_10ms_E=Test_TI_10ms;
V_10ms_TI_calc=polyval(p,Test_TI_10ms);
V_10ms_TI=V_10ms_TI_calc';
V_10ms_avg=mean(V_10ms_TI);
V_10ms_diff=(V_10ms_TI-V_10ms_avg).^2;
V_10ms_diff_avg=mean(V_10ms_diff);
TI_10ms=((sqrt(V_10ms_diff_avg))/V_10ms_avg)*100;
Test_TI_20ms=table2array(Test_TI_20ms);
Test_TI_20ms_E=Test_TI_20ms;
V_20ms_TI_calc=polyval(p,Test_TI_20ms);
V_20ms_TI=V_20ms_TI_calc';
V_20ms_avg=mean(V_20ms_TI);
V_20ms_diff=(V_20ms_TI-V_20ms_avg).^2;
V_20ms_diff_avg=mean(V_20ms_diff);
TI_20ms=((sqrt(V_20ms_diff_avg))/V_20ms_avg)*100;
figure (2)
plot(t,V_10ms_TI,'b')
hold on
yline(V_10ms_avg,'r')
xlabel('time (s)');
ylabel('Velocity (m/s)');
title 'Velocity vs Time for 10 m/s Using Hot Wire'
hold off
figure (3)
plot(t,V_20ms_TI,'b')
hold on
yline(V_20ms_avg,'r')
xlabel('time (s)');
ylabel('Velocity (m/s)');
title 'Velocity vs Time for 20 m/s Using Hot Wire'
hold off
%%

%% Flow Uniformity %%
% Flow Uniformity at 5 m/s %
Test_FU_5ms_0in=table2array(Test_FU_5ms_0in);
Test_FU_5ms_0in_E=Test_FU_5ms_0in;
V_FU_5ms_0in_calc=polyval(p,Test_FU_5ms_0in_E);
V_FU_5ms_0in=V_FU_5ms_0in_calc';
V_FU_5ms_0in_avg=mean(V_FU_5ms_0in);
Test_FU_5ms_5in=table2array(Test_FU_5ms_5in);
Test_FU_5ms_5in_E=Test_FU_5ms_5in;
V_FU_5ms_5in_calc=polyval(p,Test_FU_5ms_5in_E);
V_FU_5ms_5in=V_FU_5ms_5in_calc';
V_FU_5ms_5in_avg=mean(V_FU_5ms_5in);
Test_FU_5ms_10in=table2array(Test_FU_5ms_10in);
Test_FU_5ms_10in_E=Test_FU_5ms_10in;
V_FU_5ms_10in_calc=polyval(p,Test_FU_5ms_10in_E);
V_FU_5ms_10in=V_FU_5ms_10in_calc';
V_FU_5ms_10in_avg=mean(V_FU_5ms_10in);
Test_FU_5ms_13in=table2array(Test_FU_5ms_13in);
Test_FU_5ms_13in_E=Test_FU_5ms_13in;
V_FU_5ms_13in_calc=polyval(p,Test_FU_5ms_13in_E);
V_FU_5ms_13in=V_FU_5ms_13in_calc';
V_FU_5ms_13in_avg=mean(V_FU_5ms_13in);
Test_FU_5ms_14in=table2array(Test_FU_5ms_14in);
Test_FU_5ms_14in_E=Test_FU_5ms_14in;
V_FU_5ms_14in_calc=polyval(p,Test_FU_5ms_14in_E);
V_FU_5ms_14in=V_FU_5ms_14in_calc';
V_FU_5ms_14in_avg=mean(V_FU_5ms_14in);
Test_FU_5ms_15in=table2array(Test_FU_5ms_15in);
Test_FU_5ms_15in_E=Test_FU_5ms_15in;
V_FU_5ms_15in_calc=polyval(p,Test_FU_5ms_15in_E);
V_FU_5ms_15in=V_FU_5ms_15in_calc';
V_FU_5ms_15in_avg=mean(V_FU_5ms_15in);
x=[0 5 10 13 14 15];
V_FU_total=[V_FU_5ms_0in_avg V_FU_5ms_5in_avg V_FU_5ms_10in_avg V_FU_5ms_13in_avg V_FU_5ms_14in_avg V_FU_5ms_15in_avg];
figure (4)
plot(x,V_FU_total,'b')
xlabel('Distance From Center (in)');
ylabel('Velocity (m/s)');
title 'Velocity vs Distance From Center for 5 m/s Using Hot Wire'
%%

% Flow Uniformity at 10 m/s %
Test_FU_10ms_0in=table2array(Test_FU_10ms_0in);
Test_FU_10ms_0in_E=Test_FU_10ms_0in;
V_FU_10ms_0in_calc=polyval(p,Test_FU_10ms_0in_E);
V_FU_10ms_0in=V_FU_10ms_0in_calc';
V_FU_10ms_0in_avg=mean(V_FU_10ms_0in);
Test_FU_10ms_5in=table2array(Test_FU_10ms_5in);
Test_FU_10ms_5in_E=Test_FU_10ms_5in;
V_FU_10ms_5in_calc=polyval(p,Test_FU_10ms_5in_E);
V_FU_10ms_5in=V_FU_10ms_5in_calc';
V_FU_10ms_5in_avg=mean(V_FU_10ms_5in);
Test_FU_10ms_10in=table2array(Test_FU_10ms_10in);
Test_FU_10ms_10in_E=Test_FU_10ms_10in;
V_FU_10ms_10in_calc=polyval(p,Test_FU_10ms_10in_E);
V_FU_10ms_10in=V_FU_10ms_10in_calc';
V_FU_10ms_10in_avg=mean(V_FU_10ms_10in);
Test_FU_10ms_13in=table2array(Test_FU_10ms_13in);
Test_FU_10ms_13in_E=Test_FU_10ms_13in;
V_FU_10ms_13in_calc=polyval(p,Test_FU_10ms_13in_E);
V_FU_10ms_13in=V_FU_10ms_13in_calc';
V_FU_10ms_13in_avg=mean(V_FU_10ms_13in);
Test_FU_10ms_14in=table2array(Test_FU_10ms_14in);
Test_FU_10ms_14in_E=Test_FU_10ms_14in;
V_FU_10ms_14in_calc=polyval(p,Test_FU_10ms_14in_E);
V_FU_10ms_14in=V_FU_10ms_14in_calc';
V_FU_10ms_14in_avg=mean(V_FU_10ms_14in);
Test_FU_10ms_15in=table2array(Test_FU_10ms_15in);
Test_FU_10ms_15in_E=Test_FU_10ms_15in;
V_FU_10ms_15in_calc=polyval(p,Test_FU_10ms_15in_E);
V_FU_10ms_15in=V_FU_10ms_15in_calc';
V_FU_10ms_15in_avg=mean(V_FU_10ms_15in);
x=[0 5 10 13 14 15];
V_FU_total_10ms=[V_FU_10ms_0in_avg V_FU_10ms_5in_avg V_FU_10ms_10in_avg V_FU_10ms_13in_avg V_FU_10ms_14in_avg V_FU_10ms_15in_avg];
figure (5)
plot(x,V_FU_total_10ms,'b')
xlabel('Distance From Center (in)');
ylabel('Velocity (m/s)');
title 'Velocity vs Distance From Center for 10 m/s Using Hot Wire'
%%

%% Spectral Content %%
figure (6)
Fs_10 = 20000;            % Sampling frequency                   
T_10 = 1/Fs_10;             % Sampling period      
L_10 = 100000;             % Length of signal
t_10 = (0:L_10-1)*T_10;        % Time vector
X_10 = V_10ms_TI - 10;     % Signal from Hotwire
Y_10 = fft(X_10);
P2_10 = abs(Y_10/L_10);
P1_10 = P2_10(1:L_10/2+1);
P1_10(2:end-1) = 2*P1_10(2:end-1);
f_10 = Fs_10/L_10*(0:(L_10/2));
plot(f_10,P1_10,"LineWidth",3)
title("Single-Sided Amplitude Spectrum of X(t) for 10 m/s")
xlabel("f (Hz)")
ylabel("|P1(f)|")
figure (7)
Fs_20 = 20000;            % Sampling frequency                   
T_20 = 1/Fs_20;             % Sampling period      
L_20 = 100000;             % Length of signal
t_20 = (0:L_20-1)*T_20;        % Time vector
X_20 = V_20ms_TI - 20;     % Signal from Hotwire
Y_20 = fft(X_20);
P2_20 = abs(Y_20/L_20);
P1_20 = P2_20(1:L_20/2+1);
P1_20(2:end-1) = 2*P1_20(2:end-1);
f_20 = Fs_20/L_20*(0:(L_20/2));
plot(f_20,P1_20,"LineWidth",3)
title("Single-Sided Amplitude Spectrum of X(t) for 20 m/s")
xlabel("f (Hz)")
ylabel("|P1(f)|")