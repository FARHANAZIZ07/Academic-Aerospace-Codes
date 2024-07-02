clc;
clear all;
w = [0.0 0.0 0.0 -2.49479E-02	-2.49473E-02	-2.49473E-02	1.68267E-02	9.60586E-03	1.64595E-03	7.90909E-02	6.75550E-02	5.62843E-02];
q = [3.90467E+0 1.10640E+02	8.85940E+01	3.82591E+02	1.07253E+02	8.56026E+01	3.63625E+02	9.90872E+01	7.84376E+01	6.32280E+02	1.59865E+02	1.23717E+02];
x = [6.0 24.0 42.0 6.0 24.0 42.0 6.0 24.0 42.0];
y = [0.0 0.0 0.0 60.0 60.0 60.0 120.0 120.0 120.0];
N = 9;
D = 1;
for i = 1:N
    for j = 1:N
        r(i) = (x(j) - x(i)).^2 + (y(j) - y(i)).^2;
        %             K((i),(j)) = 1./(16.*pi.*D).*(r(i)).*log(r(i));
        K((i),(j)) = (r(i)).*log(r(i));
    end
end
K(isnan(K))= 0

format shortG
for i = 1:12
    for j = 1:12
        CMAT(i,j) = 0;
        if i == 1 && j >= 4
            CMAT(i,j) = 1;
        end
        if  j == 1 && i >= 4
            CMAT(i,j) = 1;
        end
        if i ==2 && j>=4
            CMAT(i,j) = x(j-3);
        end
        if i == 3 && j>=4
            CMAT(i,j) = y(j-3);
        end
        if j ==2 && i>=4
            CMAT(i,j) = x(i-3);
        end
        if j == 3 && i>=4
            CMAT(i,j) = y(i-3);
        end
        if (i > 3 && j > 3)
            CMAT(i,j) = K(i-3,j-3);
        end
    end
end

CMAT


xa = [6.0 18.0 36.0 6.0 18.0 36.0 6.0 18.0 36.0 6.0 18.0 36.0];
ya = [12.0 12.0 12.0 36.0 36.0 36.0 60.0 60.0 60. 96.0 96.0 96.0];

for i = 1:12
    for j = 1:9
        ra(i) = (x(j) - xa(i)).^2 + (y(j) - ya(i)).^2;
        %    k((i),(j)) = 1./(16.*pi.*D).*(r(i)).*log(r(i));
        Ka((i),(j)) = (ra(i)).*log(ra(i));
    end
end
Ka

for i = 1:12
    for j = 1:12
        %         KMAT(i,j) = 0;
        if i >= 0 && j == 1
            KMAT(i,j) = 1;
        end
        if  j == 2 && i >= 0
            KMAT(i,j) = xa(i);
        end
        if i >= 0 && j == 3
            KMAT(i,j) = ya(i);
        end
        if (i >= 0 && j >= 4)
            KMAT(i,j) = Ka(i,j-3);
        end
    end
end
KMAT(isnan(KMAT))= 0
SPLINEload = transpose(KMAT*(inv(CMAT)))*transpose(q)

SPLINEload(isnan(SPLINEload))= 0





AERO_LOADS = sum(q)

AERO_LOADs_struct = sum(SPLINEload(4:12),1)

%% Displacments;

splinedisp = (KMAT*(inv(CMAT))*transpose(w))
Struct_dis = sum(w)
aer_disp = sum(splinedisp)



%% SLOPE ;
%% slope;






for i = 1:12
    for j = 1:9
        ra(i) = (x(j) - xa(i)).^2 + (y(j) - ya(i)).^2;
        k_a((i),(j)) = 2.*(x(j)-xa(i)).*(1+log(ra(i)));   %% Dk = dk/xa = 2.*(x(i)-xa(j)).*(1+log(ra(i)));
    end
end
k_a

for i = 1:12
    for j = 1:12
        %         KMAT(i,j) = 0;
        if i >= 0 && j == 1
            K_MAT(i,j) = 0;
        end
        if  j == 2 && i >= 0
            K_MAT(i,j) = 1;
        end
        if i >= 0 && j == 3
            K_MAT(i,j) = 0;
        end
        if (i >= 0 && j >= 4)
            K_MAT(i,j) = k_a(i,j-3);
        end
    end
end

K_MAT


slope_D = K_MAT*inv(CMAT)*transpose(w)

a= rad2deg(slope_D)