clc
clear all
syms R xtheta rtheta r0 b wtheta
format default
wh = 12
wtheta = 30
b = 40
R = wh./wtheta
rtheta = 0.45.*b.^2
xbar = 0.09
xtheta = 0.09.*b
bxtheta = b.*xtheta
Bxtheta = linspace(-1,1,100).*b
r0 = (rtheta.^2- bxtheta.^2).^0.5
Ro = (rtheta.^2- Bxtheta.^2).^0.5
omegasq1 = (1 + R^2 + sqrt((1-R^2)^2+4*R^2*(bxtheta/rtheta)^2)) /(2*(r0/rtheta)^2)
omegasq2 = (1 + R^2 - sqrt((1-R^2)^2+4*R^2*(bxtheta/rtheta)^2)) /(2*(r0/rtheta)^2)
for i= 1:length(Bxtheta)
Omega1(i) = ((1 + R.^2) - (sqrt((1 - (R.^2)).^2 + 4.*(R^2).*(Bxtheta(i).*b ./rtheta).^2)))./(2.*(Ro(i)./rtheta).^2);
Omega2(i) = ((1 + R.^2) + (sqrt((1 - (R.^2)).^2 + 4.*(R^2).*(Bxtheta(i).*b ./rtheta).^2)))./(2.*(Ro(i)./rtheta).^2);
end
plot(Bxtheta,Omega2)
hold on
plot(Bxtheta,Omega1)
yline(0)
xlabel('b*x\theta')
ylabel('omega')
hold off
omega1 = sqrt(omegasq1)
omega2 = sqrt(omegasq2)
A = [omega1.^2 + R.^2, -omega1.^2 .*xtheta; -omega1.^2 .*xtheta, (-omega1.^2 + 1).*rtheta.^2]
B = [omega2.^2 + R.^2, -omega2.^2 .*xtheta; -omega2.^2 .*xtheta, (-omega2.^2 + 1).*rtheta.^2]
phi1 = (omega1^2+R^2-omega1^2*bxtheta)/ (omega1^2*bxtheta - ((-omega1^2+1)*rtheta^2))
phi_ratio1 = omega1^2*xbar/(R^2-omega1^2)
phi_ratio2 = omega2^2*xbar/(R^2-omega2^2)
ubar_12 = 1/phi_ratio2
Phi1 = [1;ubar_12]
Phi2 = [phi_ratio1;1]
rtheta_b = 0.45
Cl = 6.5
mu = 15
a = 0.2
bxa = (b - 0.25.*(2*b))
eb = (bxa + b*(a))./(b)
db = eb + xbar
Vb = linspace(0,2.5,100)
A = rtheta_b.^2 - xbar.^2
B = rtheta_b.^2.*(1+R.^2) - ((db.*Cl)./(pi.*mu)).*Vb.^2
C = rtheta_b.^2.*R.^2 - ((R.^2.*Cl.*eb)./(mu.*pi))*Vb.^2
O1 = (B + sqrt((B.^2) - (4.*A.*C)))./(2.*A)
O2 = (B - sqrt((B.^2) - (4.*A.*C)))./(2.*A)
B1 = rtheta_b.^2.*(1+R.^2) - ((db.*Cl)./(pi.*mu)).*Vb.^2
C1 = rtheta_b.^2.*R.^2 - ((R.^2.*Cl.*eb)./(mu.*pi)).*Vb.^2
O11 = (B1 + sqrt((B1.^2) - (4.*A.*C1)))./(2.*A)
O22 = (B1 - sqrt((B1.^2) - (4.*A.*C1)))./(2.*A)
plot(Vb, real(O1), 'o')
hold on
plot (Vb,real(O2), 'o')
grid on
hold off
ylabel ("Omega")
xlabel ("V_{bar}^2") 
title ('Real Omega values vs V_{bar}^2')
plot(Vb, imag(O11),'o')
hold on
plot (Vb,imag(O22), 'o')
grid on
hold off
xlim([0.936 1.812])
ylim([-0.242 0.389])
ylabel ("Omega")
xlabel ("V_{bar}^2") 
title ('Imaginary Omega values vs V_{bar}^2')
D2 = B1.^2 - 4.*A.*C1
plot(Vb,real(D2), 'ro')
grid on
yline(0)
xlim([0.746 2.031])
ylim([-0.0242 0.0475])
ylabel (" B^2 - 4AC")
xlabel ("V_{bar}^2") 
title ('B^2 - 4AC vs V_{bar}^2')