%% 
% QUESTION1
% 
% a)

%p1po=(((((G+1)/2)*(M^2))/(1+(((G-1)/2)*(M^2))))^(G/(G-1)))/((((2*G)/(G+1))*(M^2))-((G-1)/(G+1)))^(1/(G-1))
%AoAot=(1/M)*((2/(G+1))*(1+((G-1)/2)*M^2))^((G+1)/(2*(G-1)))
% Define the equation with variables M and G
eqn = @(M,G) ((((((G+1)/2)*(M^2))/(1+(((G-1)/2)*(M^2))))^(G/(G-1)))/((((2*G)/(G+1))*(M^2))-((G-1)/(G+1)))^(1/(G-1))*(1/M)*((2/(G+1))*(1 + ((G-1)/2)*M^2))^((G+1)/(2*(G-1))));
AcAt = 1.17;%  target result of the equation
G = 1.4;%  fixed value for G
M = [1.5 4]; % initial guess for M
Mo = fzero(@(M) eqn(M, G) - AcAt, M);
% Display the final value of M
disp(['value of Mo: ' num2str(Mo)]);

%% 
% b)

AoAos=(1/Mo)*((2/(G+1))*(1 + ((G-1)/2)*Mo^2))^((G+1)/(2*(G-1)))% Ao/Aost
AtAts=AoAos/AcAt
% Define the equation with variables M and G
target = AtAts;%  target result of the equation
G = 1.4;%  fixed value for G
eqn = @(M) ((1/M)*((2/(G+1))*(1 + ((G-1)/2)*M^2))^((G+1)/(2*(G-1))))- target;


%M = [-1 4]; % initial guess for M
Mt = fzero(eqn,1.6);
% Display the final value of M
disp([' value of Mt: ' num2str(Mt)]);
%% 
% c)

poypoO=(((((G+1)/2)*(Mt^2))/(1+(((G-1)/2)*(Mt^2))))^(G/(G-1)))/((((2*G)/(G+1))*(Mt^2))-((G-1)/(G+1)))^(1/(G-1))%Poy/PoO
%% 
% d)

AoAos=AcAt;%Ao/Aos=Ao/At
% Define the equation with variables M and G
eqn = @(M,G) ((1/M)*((2/(G+1))*(1 + ((G-1)/2)*M^2))^((G+1)/(2*(G-1))));
target = AoAos;%  target result of the equation
G = 1.4;%  fixed value for G
M = [-1 4]; % initial guess for M
Mu= fzero(@(M) eqn(M, G) - target, M);
% Display the final value of M
disp([' value of Mu: ' num2str(Mu)]);
%% 
% e)

poypoO=(((((G+1)/2)*(Mu^2))/(1+(((G-1)/2)*(Mu^2))))^(G/(G-1)))/((((2*G)/(G+1))*(Mu^2))-((G-1)/(G+1)))^(1/(G-1))% Poy/PoO when mu= 1.49
AAs=(1/Mu)*((2/(G+1))*(1+((G-1)/2)*Mu^2))^((G+1)/(2*(G-1)))%A/A* when Mu=1.49
AoAt=AAs*poypoO%Ao/At=A/A*xPoy/PoO
FS=(AcAt-AoAt)/AcAt%Fraction Spilled=(Ac/At-Ao/At)/Ac/At