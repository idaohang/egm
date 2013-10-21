function [vn,Gk,Gc,Gl,diff] = valfuniter(Gk1,zt,piz,vn,Gk,Gc,Gl)
% Standard value function iteration
%
% Created:
% 18.10.2013, Antti Ripatti
%
% based on the fortran code by Barillas and Fernandez Villaverde (2006)
%
global alpha beta delta theta tau;
gSize = length(Gk1);
zStates = size(zt,2);
% allocating matrices
Glprov = Gl;
Gkprov = Gk;
vnprov = zeros(size(vn));
diff = 100;
for igrid = 1:gSize
  k = Gk1(igrid);
  for iz = 1:zStates
    if iz == 1
      prev_kp = 1;
    else
      prev_kp = index_max_k_so_far;
    end;
    z = zt(igrid,iz);
    value_max_so_far = -1e8;
    flag = 1;
    for ikp = prev_kp:gSize
fprintf(1,'.');
      kp = Gk1(ikp);
      lp = mynewton(kp,k,z);
      c = exp(z)*(k^alpha)*(lp^(1-alpha))+(1-delta)*k - kp;
      if c > 0
        value_comp = (1-beta)*(((c^theta)*((1-lp)^(1-theta)))^(1-tau))/(1-tau) ...
          +beta.*piz(iz,:)*vn(ikp,:)';
        if value_comp > value_max_so_far
          value_max_so_far = value_comp;
					Glprov(ikp,iz) = lp;
					Gkprov(ikp,iz) = kp;
					index_max_k_so_far = ikp;
        else
          if flag == 1 && gSize > 20 
            [kp, lp, value_refinement] = refinement(k,z,piz,vn,Gk1,ikp,iz);
            flag = flag+1;
						if value_max_so_far < value_refinement 
							value_max_so_far = value_refinement;
							Glprov(ikp,iz) = lp;
							Gkprov(ikp,iz) = kp;
            end;
          end;
        end;
      else
        disp('hypataan');
        return;
      end;
    end;
    vnprov(ikp,iz) = value_max_so_far;
  end;
end;
disp('tultiin');
diff  = max(abs(vnprov - vn));
diff2 = max(abs(g_lprov - g_l));
vn = vnprov;
Gl = Glprov;
Gk = Gkprov;
for ik = 1:gSize;
	for iz = 1:zStates;
    Gc(ik,iz)=exp(z(iz))*Gk1(ik)^alpha*Gl(ik,iz)^(1-alpha) + ...
      (1-delta)*Gk1(ik) - Gk(ik,iz);
	end;
end;
fprintf(1,'valfuniter; sup norm %g, labor change %g',diff, diff2);