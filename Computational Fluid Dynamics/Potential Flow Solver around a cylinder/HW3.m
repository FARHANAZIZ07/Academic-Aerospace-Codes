clc
close all

U_inf = 1.0;                 %velocity at inlet
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−Mesh Size
Nr=50; % grid count in R−dimension
Nth=50; % grid count in theta−dimension
%−−−−−−−−−−−−−−−−−−−−−−−−−SOR relaxation parameter
w_opt = 2/(1+sqrt(1-cos(pi/(Nr+1))*cos(pi/(Nth+1))));
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−Other Options
tol=10^(-8); % machine tolerance
%iter =10000; % maximum iterations
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−Dimensions
Rcyl = 0.5;
Rmax = 10;
dr=(Rmax-Rcyl)/(Nr-1);
dth=2*pi/(Nth-1);
%−−−−−−−−−−−−−−−−−−−−−−−−−Boundary and initial Conditions
Psi = zeros(Nr,Nth);
for j = 1:Nr
    for k=1:Nth
rj(j,k)=(j-1)*dr+Rcyl;
theta_k(j,k)=(k-1)*dth;
x(j,k)=rj(j,k)*cos(theta_k(j,k));
y(j,k)=rj(j,k)*sin(theta_k(j,k));
        if (j==Nr)
            Psi (j,k)=U_inf*y(j,k);
        end
    end
end
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−SOR Method
res(1)=1;
iter=1;
54
while res(iter)>tol
iter=iter+1;
itcount=[iter]
    for j=2:(Nr-1)
        for k=2:(Nth-1)
        num=1/(2/dr^2+2/(rj(j,k)^2*dth^2));
Delta_Psi(j,k)=num*((Psi(j+1,k)-Psi(j-1,k))/(2*rj(j,k)*dr)+(1/dr^2)*(Psi(j+1,k)+Psi(j-1,k))+(1/(rj(j,k)^2*dth^2))*(Psi(j,k+1)+Psi(j,k-1)))-Psi(j,k );
Psi(j,k)=Psi(j,k)+w_opt*Delta_Psi(j,k);
        end
    end
res(iter)=norm(Delta_Psi);
disp(res(iter));
end
disp(iter);
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−Ploting the residual

figure(1)
semilogy(1:iter,res,'-k','LineWidth',2)
xlabel('Number of Iterations')
ylabel('residual')
title('Convergence History')
grid on
% pause
% set(ya,'FontSize',20);
set(gca,'FontSize',30);
ax=gca ;
c=ax.Color;
ax.Color='w';
ax.Box='on';
ax.LineWidth=2;
%−−−−−−−−−−−−−−−−−−−−−Calculating the radial and tangnetial velocities

for j=1:Nr
    for k=1:Nth
        if ((j==1)|( j==Nr))
            if j==1
%Forward Difference
%utheta(j,k)=(psi(j,k)−Psi(j+1,k))/(dr);
utheta(j,k)=(Psi(j+2,k)-4*Psi(j+1,k)-3*Psi(j,k))/(2*dr);
                else
%Backward Difference
%utheta(j,k)=(Psi(j−1,k)−Psi(j,k))/(dr);
utheta(j,k)=(4*Psi(j-1,k)-3*Psi(j,k)-Psi(j-2,k))/(2*dr);
                end
            else
%Centered Difference
% Comment utheta=−dpsi/dr Centered difference
utheta(j,k)=(Psi(j-1,k)-Psi(j+1,k))/(2*dr);
%utheta(j,k)=(Psi(j+2,k)−8∗Psi(j+1,k)+8∗Psi(j−1,k)−Psi(j−2,k))/(12∗dr);% My Heavy mi s t ake
            end
            if((k==1)|(k==Nth))
                if k==1
        % ur( j , k )=(P si ( j , k+1)−P si ( j , k ) ) / ( dth∗ r j ( j , k ) ) ;
                  else
    %Centered D i f f e r e n c e
ur(j,k)=(Psi(j,k+1)-Psi(j,k-1))/(2*dth*rj(j,k));
%ur( j , k )=(−P si ( j , k+2)+8∗P si ( j , k+1)−8∗P si ( j , k−1)+P si ( j , k−2))/(12∗dth∗rj(j,k));% My Heavy mi s t ake
                end
            end
    end
end
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−Changing from polor to cartesian co o r di n a t e −−−−−−−−−−−−−−−−−−−−−−
[u1,u2]=pol2cart(ur,utheta);
Cp_num=1-((u1.^2)+(u2.^2)/U_inf^2);
%−−−−−−−−−−−−−−−−−−−−−−−−Calculating the exact numerical pressure c o e f f i c i e n t s −−−−−−−−−−−−−−−−−
theta=theta_k(1,(1:Nth));
for n=1:Nth
    Cp_ex(n)=1-4*sin(theta(n)).^2;
end
%−−−−−−−−−−−−−−−−−−−−−−−−Pl o t ti n g the e x a c t and n um e ri c al p r e s s u r e
%coefficients −−−−−−−−−−−−−−−−−%
 figure ( 2 )
 plot (theta , Cp_num (1,1:Nth) , '-k' , 'LineWidth' ,2 )
 hold on
 plot (theta , Cp_ex , '−−k ' , 'LineWidth ' , 2 )
 xlabel ( '\ t h e t a ' )
 ylabel ( 'Cp ' )
 title ( 'Cp vs . \ theta ' )
 legend ( 'Numerical ' , 'Exact ' )
 grid on
 set( gca , 'FontSize ' , 30 ) ;
 ax = gca ;
 c = ax . Color ;
 ax . Color = 'w ' ;
 ax . Box = 'on ' ;
 ax . LineWidth = 2 ;
 % pause
 figure ( 3 )
 pcolor( Psi )

 %−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−C al c ul a ti n g the e r r o r
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

 err2=(norm (( Cp_ex - Cp_num(1,(1:Nr)))))./( sqrt(Nr.*Nth));
 %−−−−−−−−−−−−−−−−−−−−−−−−−−−P r e s s u r e D i s t r i b u t i o n around the e n d t i r e
% domain−−−−−−−−−−−−−−−−−−−−

 Contours =91; % number o f s t r e a m l i n e s t o pl o t
 Grid='on' ; % tu rn c f d g ri d p oi n t s on/ o f f
 colordef white % Add c y l i n d e r t o s t r e amli n e pl o t
 plotcirc=linspace ( 0 , 2.*pi ,100 ) ;
 xcirc=Rcyl .* cos ( plotcirc ) ;
 ycirc=Rcyl .* sin ( plotcirc ) ;
 plot ( xcirc , ycirc , 'b' ) ;
 hold on
 xlabel ('itcount')
 ylabel ('y')
 title ( 'CFD−Potential Flow' )
 if Grid=='on'
 for j =1:Nr
 for k=1:Nth
 plot (x(j,k),y(j,k),'−')
 hold on
 end
 end
 end
 %
 xplot=x ;
 yplot=y ;
 Prplot=Psi ;
 nplot=Nr;
 for j =1:Nr
 Prplot(j,Nth+1)=Psi(j,1); 
 xplot (j,Nth+1) = x(j ,1) ;
 yplot (j,Nth+1) = y(j,1) ;
 end
 for j =1:(Nr)
 for k=1:Nth+1
 Prplot2 (j,k)=Prplot (j,k) ;
 xplot2 (j,k)=xplot (j,k) ;
 yplot2 (j,k)=yplot (j,k) ;
 end
 end
 hold on
 contour ( xplot2 , yplot2 , Prplot2 , Contours)
 contour ( xplot , yplot , Prplot , Contours)
 title ( 'Streamlines' )
 grid off
 fill ( xcirc , ycirc , 'w' )
 axis ([-Rcyl.*5.0 Rcyl.*5.0 -Rcyl.*5 Rcyl.*5 ] )
 axis equal
 % pause
 set ( gca , 'F on tSize ' , 30 ) ;
 ax = gca ;
 c = ax . Color ;
 ax . Color = 'w ' ;
 ax . Box = 'on ' ;
 ax . LineWidth = 2 ;
 hold off
 %
 figure (4)
 wplot=Cp_num ;
 for j =1:Nr
 wplot ( j , Nth+1)=Cp_num( j , 1) ;
 end
 pcolor ( xplot , yplot , wplot )
 shading interp
 colorbar
 axis ([-Rcyl.*5.0  Rcyl .*5.0 -Rcyl.*5 Rcyl*5] )
 axis equal
 hold on
 grid off
 title ( 'PressureCoefficient' )
 fill( xcirc , ycirc , 'w')
 set( gca , 'FontSize' , 30 ) ;
 ax = gca ;
 c = ax . Color ;
 ax . Color = 'w ' ;
 ax . Box = 'on' ;
 ax . LineWidth = 2 
    