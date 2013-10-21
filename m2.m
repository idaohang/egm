% MODEL2 in Barillas - Fernandez Villaverde (2007)
% Solved with Endogenous Grid Method by Carroll (2005)
%
% Created:
% 15 Oct 2013, Antti Ripatti
%
% (c) Antti Ripatti, 2013-
%
% Special thanks to Gero Dolfus in patiently explaining me details
%
%% Constants 
clear('all'); clear('global','all'); close('all');
global tau delta beta alpha theta MaxIt;
gSize = 10000;
options = optimset('Display','off');
relativeKss = 0.25; % Range limits [(1-r)*Kss, (1+r)*Kss] relative to capital stock steady-state
DispEveryIter = 100; % how often iteration results are displayed
MaxIt = 10000; % Maximum number of iterations in EGM
tolerance = 1e-8; % tolerance for iteration
nWorkers = 4; % number of workers in parallell computing
% parameters
delta  = 0.0196; % Depreciation
alpha  = 0.4; % Capital Share 
tau    = 2; % risk aversion
beta   = 0.9896; % Discount factor
theta  = 0.357; % CD parameter in labour
rho    = 0.95; % Autoregressive 
sigma  = 0.007; % Variance
%% Steady states
psi = (1/alpha*(1/beta-1+delta))^(1/(1-alpha));
Psi = theta/(1-theta)*(1-alpha)*psi^(-alpha);
%Omega = psi^(1/alpha) - delta; % as in the paper
Omega = psi^(1-alpha) - delta; % as in the code
k_ss = Psi/Omega + psi*Psi;
l_ss = psi*k_ss;
c_ss = Omega*k_ss;
y_ss = k_ss^alpha*l_ss^(1-alpha);
disp('Steady state values');
fprintf(1,'Consumption: %4.2f\n',c_ss);
fprintf(1,'Output:      %4.2f\n',y_ss);
fprintf(1,'Capital:     %4.2f\n',k_ss);
fprintf(1,'Labour:      %4.2f\n',l_ss);
%% Grid choice
kMin = (1-relativeKss)*k_ss;
kMax = (1+relativeKss)*k_ss;
kStep = (kMax-kMin)/gSize;
zStates = 2;
piz = [0.8 0.2;
       0.2 0.8];
z = ones(gSize,zStates)*[0 0; 0 -1];
%% Initialize grids
Gk1 = kMin:kStep:kMax;
Gk1 = Gk1(1:gSize); % truncate
% Initialize Vtilde0, 
Vn = zeros(gSize,zStates);
for i = 1:gSize;
  temp = i/gSize + (1/(1-beta))*(c_ss^theta*(1-l_ss)^(1-theta))^(1-tau)/(1-tau);
  Vn(i,:) = temp.*ones(1,zStates);
end;
% cash in hand; This is the G_Y(t+1)
cih = zeros(gSize,zStates);
for i = 1:zStates;
  cih(:,i) = (exp(z(:,i)').*(Gk1.^alpha).*(l_ss^(1-alpha))+(1-delta).*Gk1)';
end;
% Action starts with initialization
normi = 1;
iter = 0;
Yend = zeros(gSize,zStates);
Vend1 = Yend;
Vn1 = Vn;
%% Step1: Labour fixed to l_ss. Simple EGM iterations
disp('STEP1: Labour fixed to l_ss. Simple EGM iterations');
while normi>tolerance && iter<MaxIt ;
  iter = iter + 1;
  % Derivatives at grid points; glumsy implementation
  D = zeros(gSize,zStates);
  D(1,:)=(Vn(2,:)-Vn(1,:))/kStep;
	D(end,:) = (Vn(end,:)-Vn(end-1,:))/kStep;
	for i = 2:gSize-1  % this should be vectorized
		D(i,:) = (Vn(i+1,:)-Vn(i-1,:))/(2*kStep);
  end;
  Cstar = (D.*(1-l_ss)^((theta-1)*(1-tau))/theta).^(1/(theta*(1-tau)-1));
  Vend = (Cstar.^theta*(1-l_ss).^(1-theta)).^(1-tau)./(1-tau) + Vn;
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
if iter == MaxIt
  warning('Iteration did not converge with %4.0f iterations (=MaxIt). The norm %12.6e\nConsider increasing MaxIt\n',iter,normi);
else
  fprintf(1,'Iteration converged with %4.0f iterations and the norm %12.6e\n',iter,normi);
end;
%% calculate k_t
disp('Calculating k_t, interpolating Cstar, Gk1, etc...');
% lets parallelize this
if matlabpool('size')==0
  matlabpool('open','local',nWorkers);
end;
kend = Yend; % just initialize
k0 = [Gk1(1) Gk1(1)];
parfor i=1:gSize % note parfor
  for j=1:zStates
    kend(i,j) = fsolve(@(x) (Yend(i,j)-exp(z(i,j))*x^alpha-(1-delta)*x),k0(j),options);
  end;
  %k0 = kend(i,:); % this cannot be used in parfor loop
end;
%% Interpolate
for i = 1:zStates
    Vend(:,i) = interpol(kend(:,i),Vend(:,i),Gk1);
    Gk1 = interpol(kend(:,i),Gk1,Gk1);
    Cstar(:,i) = interpol(kend(:,i),Cstar(:,i),Gk1);
end;
Gc = Cstar;
Gl = ones(size(Gc)).*l_ss;
Gk = kend;
%% Step 2: Standard value function iteration
disp('STEP2: standard value function iteration');
[Vn,Gk,Gc,Gl,diff] = valfuniter(Gk1,z,piz,Vn,Gk,Gc,Gl);
%% Step 3: Cycle between EGM step and standard value function iteration
disp('STEP3: Cycle between EGM step and standard value function iteration');
if diff > tolerance
  iter = 1;
	while diff < tolerance
    % The first step is take the input from the standard algorithm and compute the values of Kend such that
    % g_k(end)=k_grid
    fprintf(1,'Iteration %5.0f diff in value function iteration %8.3e',iter,diff);
    for ik = 1:gSize
      for iz = 1:nStates
          kend = golden(Gk1(ik),Gk1, Gk(:,iz));
      end;
    end;
    % Now do the endogenous grid for a few more extra steps
    Vn = Vn/(1-beta);
    Vn = EGMiter(Vn, Gk1, z, piz, kend, Gl, tolerance);
    % One standard iteration to recover policy functions.
    [Vn,Gk,Gc,Gl,diff] = valfuniter(Gk1,zt,piz,Vn,Gk,Gc,Gl);
	end;
end;


%% calculate k_t
kend = Yend; % just initialize
k0 = [Gk1(1) Gk1(1)];
parfor i=1:gSize % note parfor
  for j=1:zStates
    kend(i,j) = fsolve(@(x) (Yend(i,j)-exp(z(i,j))*x^alpha-(1-delta)*x),k0(j),options);
  end;
  %k0 = kend(i,:); % this cannot be used in parfor loop
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