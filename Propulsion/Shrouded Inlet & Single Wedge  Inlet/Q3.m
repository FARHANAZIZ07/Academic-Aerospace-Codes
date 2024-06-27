G=1.4
finalM=-100
Mos=1.82
Gmos=((G+1)/(G-1))^0.5*atan((((G-1)/(G+1))*(Mos^2-1))^0.5)-atan((Mos^2-1)^0.5)
theta=11
Gm1 = Gmos - (theta*3.14)/180
rhs=[];
lhs=Gm1;
Ms=[0:0.001:4];
for M=0:0.001:4
    temp=((G+1)/(G-1))^0.5*atan((((G-1)/(G+1))*(M^2-1))^0.5)-atan((M^2-1) ^ 0.5);
    
    rhs=[rhs temp];
end

threshold=abs(rhs-lhs);
[~,I]=min(threshold)

finalM=Ms(I)