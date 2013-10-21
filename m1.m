% MODEL1 in Barillas - Fernandez Villaverde (2007)
% Solved with Endogenous Grid Method by Carroll (2005)
%
% Created:
% 9 Oct 2013, Antti Ripatti
%
% (c) Antti Ripatti, 2013-
%
% Special thanks to Gero Dolfus in patiently explaining me details
%
%% Constants 
clear global all; close all;
global tau delta beta alpha;
gSize = 1000;
options = optimset('Display','off');
relativeKss = 0.75; % Range limits relative to capital stock steady-state
DispEveryIter = 50;
MaxIt = 10000; % Maximum number of iterations in EGM
% parameters
delta  = 0.0196; % Depreciation
alpha  = 0.4; % Capital Share 
tau    = 2; % risk aversion
beta   = 0.9896; % Discount factor
rho    = 0.95; % Autoregressive 
sigma  = 0.007; % Variance
z1 = 0;
%% Steady states
k_ss = ((1/beta-1+delta)/alpha)^(1/(alpha-1));
y_ss = k_ss^alpha;
c_ss = k_ss^alpha-delta*k_ss;
disp('Steady state values');
fprintf(1,'Consumption: %4.2f\n',c_ss);
fprintf(1,'Output:      %4.2f\n',y_ss);
fprintf(1,'Capital:     %4.2f\n',k_ss);
%% Grid choice
kMin = (1-relativeKss)*k_ss;
kMax = (1+relativeKss)*k_ss;
kStep = (kMax-kMin)/gSize;
zStates = 2;
piz = [0.8 0.2;
       0.2 0.8];
z = ones(gSize,zStates)*[0 0; 0 -1];
tolerance = 1e-5; % tolerance for iteration
%% Initialize grids
Gk1 = kMin:kStep:kMax;
Gk1 = Gk1(1:gSize); % truncate
GY1 = exp(z1).*Gk1 + (1-delta).*Gk1;
Gz = 1:2;
% Initialize Vtilde0, 
Vn = zeros(gSize,zStates);
for i = 1:gSize;
  temp = i/gSize + (1/(1-beta))*c_ss^(1-tau)/(1-tau);
  Vn(i,:) = temp.*ones(1,zStates);
end;
% cash in hand
cih = zeros(gSize,zStates);
for i = 1:zStates;
  cih(:,i) = (exp(z(:,i)').*(Gk1.^alpha)+(1-delta).*Gk1)';
end;
% Action starts with initialization
normi = 1;
iter = 0;
Yend = zeros(gSize,zStates);
Vend1 = Yend;
Vn1 = Vn;
while normi>tolerance && iter<MaxIt ;
  iter = iter + 1;
  % Derivatives at grid points; glumsy implementation
  D = zeros(gSize,zStates);
  D(1,:)=(Vn(2,:)-Vn(1,:))/kStep;
	D(end,:) = (Vn(end,:)-Vn(end-1,:))/kStep;
	for i = 2:gSize-1  % this should be vectorized
		D(i,:) = (Vn(i+1,:)-Vn(i-1,:))/(2*kStep);
  end;
  Cstar = D.^(-1/tau);
  Vend = (1/(1-tau))*Cstar.^(1-tau) + Vn;
  for i = 1:zStates
    Yend(:,i) = Cstar(:,i) + Gk1';
    Vend1(:,i) = interpol(Yend(:,i),Vend(:,i),cih(:,i));
  end;
  for i = 1:gSize % this could be vectorized too
    for j = 1:zStates
      Vn1(i,j) = beta.*piz(j,:)*Vend1(i,:)';
    end;
  end;
  normi = norm(Vn1-Vn);
  Vn = Vn1;
  if mod(iter,DispEveryIter) == 0
    fprintf(1,'Iteration %4.0f with norm %12.6e\n',iter,normi);
  end;
end;
fprintf(1,'Iteration converged with %4.0f iterations and the norm %12.6e\n',iter,normi);
%% calculate k_t 
kend = Yend; % just initialize
k0 = [Gk1(1) Gk1(1)];
for i=1:gSize
  for j=1:zStates
    kend(i,j) = fsolve(@(x) (Yend(i,j)-exp(z(i,j))*x^alpha-(1-delta)*x),k0(j),options);
  end;
  k0 = kend(i,:);
end;
%% Interpolate
for i = 1:zStates
    Vend(:,i) = interpol(kend(:,i),Vend(:,i),Gk1);
    Gk1 = interpol(kend(:,i),Gk1,Gk1);
    Cstar(:,i) = interpol(kend(:,i),Cstar(:,i),Gk1);
end;
%% Plotting results
fig = figure;
fig = plot(kend(:,1),[kend(:,1) Gk1']);
fig = title('K(t)and K(t+1) in boom');
xlabel('k(t)');
ylabel('k(t+1)');
legend({'45°','k(t+1)'});
fig = figure;
fig = plot(kend(:,2),[kend(:,2) Gk1']);
title('K(t)and K(t+1) in recession');
xlabel('k(t)');
ylabel('k(t+1)');
legend({'45°','k(t+1)'});