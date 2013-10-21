function vn = EGMiter(vn, Gk1, z, piz, kend, Gl, diffstandard)
% EGMiter Endogenous Grid Method step
%call sub_valueend2(value,grid_k,length_grid_k,y,transition,kend,g_l,diff)
% Created:
% 21.10.2013, Antti Ripatti
%
% based on the fortran code by Barillas and Fernandez Villaverde (2006)
%
global alpha beta delta theta tau MaxIt;
gSize = length(Gk1);
zStates = size(z,2);
L1 = zeros(size(kend));
Vn1 = zeros(size(L1));
iter = 1;
diff = 1000;
% Need to interpolate the policy function and the value function at the new Kend values
for iz = 1:zStates;
  L1(:,iz) = interpol(Gk1,Gl(:,iz),kend(:,iz));
end;
Vend1 = vn;
Vn = vn;
normi = 10;
while normi > 0.1*diffstardard && iter < MaxIt ;
  iter = iter + 1;
  for ik = 1:gSize
    for iz = 1:nStates
      Vn(ik,iz) = beta*piz(iz,:)*Vend(ik,:);
    end;
  end;
    % Derivatives at grid points; glumsy implementation
  D = zeros(gSize,zStates);
  D(1,:)=(Vn(2,:)-Vn(1,:))/kStep;
	D(end,:) = (Vn(end,:)-Vn(end-1,:))/kStep;
	for i = 2:gSize-1  % this should be vectorized
		D(i,:) = (Vn(i+1,:)-Vn(i-1,:))/(2*kStep);
  end;
  Cstar = (D.*(1-L1)^((theta-1)*(1-tau))/theta).^(1/(theta*(1-tau)-1));
  Vend = (Cstar.^theta*(1-L1).^(1-theta)).^(1-tau)./(1-tau) + Vn;
  for ik = 1:gSize
    for iz = 1:nStates
      kend(ik,iz) = kendogenousnewton(z(ik,iz),L1(ik,iz),Cstar(ik,iz),Gk1(ik));
    end;
  end;
  for i = 1:zStates
    Vn1(:,i) = interpol(kend(:,i),Vend(:,i),Gk1);
    L1(:,i) = interpol(Gk1(:,i),Gl(:,i),kend(:,i));
  end;
  normi = norm(Vn1-Vn);
  Vn = Vn1;
  if mod(iter,DispEveryIter) == 0
    fprintf(1,'EMG2 iteration %4.0f with norm %12.6e\n',iter,normi);
  end;
end;
vn = Vn*(1-beta);
if iter == MaxIt
  fprintf(1,'Iteration did not converge with %4.0f iterations (=MaxIt). Thehe norm %12.6e\n',iter,normi);
else
  fprintf(1,'Iteration converged with %4.0f iterations and the norm %12.6e\n',iter,normi);
end;
