clc
clear all
M=39
dx= 1/(M+1)
A= zeros(M,M)
B= 6.*dx.^3.*(1:39)'
coeff=[1 -2 1]
cc=1
A(1,1:2)= [-2 1]
A(end,end-1:end) = [1 -2]


B(1)=B(1)-0;
B(end)=B(end)-1;

for r = 2:M-1
    A(r,cc : cc+2)= coeff;
    cc= cc+1;
end

u= linsolve(A,B)
x= dx.*(1:39)
tol=1e-8
err= inf
%% POINT-JACOBI METHOD



itr1=0
while(err<tol)
    x1=x
    for i= 1:size(A)  
        sum=0
        for j=1:size(A)
            if j~=i 
        sum=sum + A(i,j).*x1(j)
        end
        x(i)=(1./A(i,i)).*(B(i)-sum)
        end
    end
    itr1=itr1+1
    err=norm(B-A.*x1)
  Res1=err./norm(B)
end
         
%% GAUSS SEIDAL METHOD



itr2=0
while(err>tol)
    x1=x
    for i= 1:size(A)  
        sum=0
        for j=1:(i-1)
        sum=sum + A(i,j).*x(j)
        end
        for j=i+1:size(A)
        sum=sum+A(i,j).*x1(j)
        end
        x(i)=(1./A(i,i)).*(B(i)-sum)
       
    end
    itr2=itr2+1
    err=norm(B-A.*x1)
    Res2=err./norm(B)
    
end

%% SOR METHOD

w=2 ./( 1 + sin(pi./(M+1)) )
itr3=0
while(err>tol)
    x1=x
    for i= 1:size(A)
        sum=0
        for j=1:(i-1)
        sum=sum + A(i,j).*x(j)
        end
        for j=i+1:size(A)
        sum=sum+A(i,j).*x1(j)
        end
        x(i)=((1-w).*x(i))+(w./A(i,i)).*(B(i)-sum)
       
    end
    itr3=itr3+1
    err=norm(B-A.*x1)
    Res3=err./norm(B)
    
end



plot(x,u,x,x.^3,x,itr1,x,itr,2,x,itr3)
legend('numerical','exact')