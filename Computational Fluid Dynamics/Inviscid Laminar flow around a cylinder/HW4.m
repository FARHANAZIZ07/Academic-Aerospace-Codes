clc; clear all; clf; 
%% Dimensions 
Rc = 0.5; 
Rm = 10; 
%%Inlet conditions
Uinf = 1; 
%%Mesh Size
Nr = 100; 
Nt = 100; 
th =(linspace(0,2*pi,Nt))';
r = linspace(Rc,Rm,Nr); 
dr = (Rm-Rc)./(Nr-1); 
dt = (2.*pi)./(Nt-1);
tol = 1e-8; 
Psi=zeros(Nr,Nt);
Psi(Nr,:) = (Uinf.*Rm.*sin(th));
%discretization 
for j = 1:Nr 
    for k=1:Nt 
        rj(j,k)=Rc+(j-1)*dr; 
        tk(j,k)=(k-1)*dt; 
        x(j,k)=rj(j,k)*cos(tk(j,k));
        y(j,k)=rj(j,k)*sin(tk(j,k));
    end 
end
%% Calculating for stream function 
omega=2./(1 + sqrt(1-(cos(pi./(Nr-1))*cos(pi./(Nt-1))))); 
Res(1)=1; 
itr=1; 
while(Res(itr)>tol)
    for j=2:(Nr-1) 
        for k=2:(Nt-1) 
            den=((((2./dr.^2))+2./(rj(j,k).^2.*dt.^2)));
dPsi(j,k)=((1/(rj(j,k)).*((Psi(j+1,k)-Psi(j-1,k))/(2*dr))+((1/dr^2)*(Psi(j+1,k)+Psi(j-1,k)))+((1.0/(rj(j,k)^2*dt^2))*(Psi(j,k+1)+Psi(j,k-1)))))/den-Psi(j,k); 
Psi(j,k)=Psi(j,k)+omega*dPsi(j,k); 
        end
    end 
    err=norm(dPsi); 
    itr=itr+1; 
    Res(itr)=err; 
    disp(Res(itr)); 
end 
figure(1) 
pcolor(x,y,Psi) 
colorbar 
figure(2) 
plot(1:itr,log10(Res)) 
grid on 
title('Convergence history') 
xlabel('No. of Iterations') 
ylabel('residual') 
Cp=zeros(Nr,Nt); 
ur=zeros(Nr,Nt); 
ut=zeros(Nr,Nt); 
%% Calculating Velocity in interior 
for j=2:Nr-1 
    for k=2:Nt-1 
        ur(j,k)=((Psi(j,k+1)-Psi(j,k-1))/2/dr)./r(j);
        ut(j,k)=(-Psi(j+1,k)+Psi(j-1,k))/2/dr; 
    end 
end 
% Calculating Velocity on boundary 
j=1; 
for k=1:Nt 
    ur(j,k)=0; 
    ut(j,k)=3*Psi(j,k)-4*Psi(j+1,k)+Psi(j+2,k)./(2*dr); 
end
j=Nr-2; 
for k=1:Nt 
    ur(j,k)=Uinf*cos(th(k)); 
    ut(j,k)=-Uinf*sin(th(k)); 
end 
k=1; 
for j=1:Nr
    ur(j,k)=(-3*Psi(j,k)+4*Psi(j,k+1)-Psi(j,k+2)./(2*dr))./r(j); 
    ut(j,k)=0;
end 
k=Nt; 
for j=1:Nr 
    ur(j,k)=(3*Psi(j,k)-4*Psi(j,k-1)+Psi(j,k-2)./(2*dr))./r(j);
    ut(j,k)=0; 
end 
%% Calculating Cpressure numerical and exact 
% Calculate Cp 
for j=1:Nr 
    for k=1:Nt 
        Cp(j,k)=1.0-(ur(j,k)^2+ut(j,k)^2)/Uinf; 
    end 
end 
for i = 1:Nt 
    Cpex(i)= 1-(4*(sin(th(i)).^2)); 
end 
figure(3) 
pcolor(x,y,Cp) 
colorbar 
title('Pressure distribution') 
%%Error calculation 
error=(norm ((Cpex-Cp(2,(1:Nr)))))./(sqrt(Nr.*Nt)); 
disp(error); 
figure(4) 
plot(th,Cp(2,1:Nt)); 
hold on; 
plot(th,Cpex); 
legend({'Cpnum','Cpexact'}) 
title('Cp vs thetha') 
xlabel('thetha') 
ylabel('Cp')