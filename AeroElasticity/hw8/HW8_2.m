clc
clear all
syms Omega k b m  a w_h_sq w_theta_sq x_theta_bar mu g_prime k Lh La Mh r_theta_bar_sq Ma 
b= 3;  %ft
m = 269.0; %slugs/ft
a = -0.2;
w_h = 10  %wh^2
w_h = 10
w_theta = 25; %wtheta^2
x_theta_bar = 0.1;
mu = 20.0;
% g_prime = 0;
k = 1./linspace(0.001,5.0,100)
1./k
r_theta_bar = 0.5; %rtheta_bar^2
Mh = 0.5;
Ma = 3./8 -1i.*(1./k);
for n = 1:length(k)
ck = 1- ((0.165.*k(n))./(k(n)-0.0455.*1i))-((0.355.*k(n))./(k(n)-0.3.*1i));
Ma = 3./8 - 1i.*(1./k(n));
Lh = 1 - ((1i.*(2.*ck))./k(n));
La = 0.5 - ((1i.*(1+2.*ck))./k(n))-(2.*ck./(k(n)).^2);
A_b =  mu.*((1-Omega.*((w_h.^2)./(w_theta.^2))))+Lh;
B_b = (mu.*x_theta_bar)+ La - Lh.*(0.5+a);
D_b = (mu.*x_theta_bar)+ Mh - Lh.*(0.5+a);
E_b = (mu.*(r_theta_bar.^2)).*(1-Omega) - (Mh.*(0.5+a)) + Ma - La.*(0.5+a) + Lh.*((0.5+a).^2);
Eq = A_b.*E_b - B_b.*D_b;
omega(:,n) = eval(solve(Eq,Omega));
end
omega
w = w_theta./sqrt(real(omega))
g = real(((w_theta.^2).*1i - omega.*(w.^2).*1i)./(w_theta.^2))
V = w.*b./k
plot(V(2,:),g(2,:),"ro")
hold on
plot(V(1,:),g(1,:)," bo ")
grid on
Vf = interp1([-0.00495264 0.00425538],[165.721 166.692],0)
Vf = 166.2433
xline(Vf)
yline(0)
xlabel('Velocity - ft/s')
ylabel('Damping - g')
title('V-G diagram')
hold off  
plot(V(1,:),w(1,:),"ro")
hold on
plot(V(2,:),w(2,:),"bo")
title('V-\omega diagram')
xlabel('Velocity - ft/s')
ylabel('frequency - \omega')
grid on