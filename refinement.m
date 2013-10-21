function [kp, lp, value_max_so_far] = refinement(k,z,piz,vn,Gk1,ikp,iz)
% refinement, value function refinement
%
% Created:
% 21.10.2013, Antti Ripatti
%
% based on the fortran code by Barillas and Fernandez Villaverde (2006)
%
global alpha beta delta theta tau;
gSize = length(Gk1);
% allocating matrices 
% algorithm
ikp_now = max(ikp-1,1);
kp = Gk1(ikp_now);
lp = mynewton(kp,k,z);
c = exp(z)*(k^alpha)*(lp^(1-alpha))+(1-delta)*k-kp;
value_max_so_far = (1-beta)*(((c^theta)*((1-lp)^(1-theta)))^(1-tau))/(1-tau) + ...
											beta*piz(iz,:)*vn(ikp_now,:)';
ikp_low  = max(1,ikp-2);
ikp_high = min(gSize,ikp);
%length_index = ikp_high-ikp_low+1; % should be used by the interpolator
kp_low  = Gk1(ikp_low);
kp_high = Gk1(ikp_high);
kp_guess = kp_low+0.6*(kp_high-kp_low);
while abs(kp_high-kp_low)<0.000001
  myvn = interpol(Gk1(ikp_low:ikp_high),vn(ikp_low:ikp_high,:),kp_guess);
	lp_guess = mynewton(kp_guess,k,z);
  c = exp(z)*(k^alpha)*(lp_guess^(1-alpha))+(1-delta)*k-kp_guess;
	if (c<0.0)
			c = 0.00000001;
	end;
  vn_comp = (1-beta)*(((c^theta)*((1-lp_guess)^(1-theta)))^(1-tau))/(1-tau) + ...
            beta*piz(iz,:)*myvn';
  if vn_comp>value_max_so_far
		if kp_guess < kp
			kp_high = kp;
		else
			kp_low = kp;
		end;
		kp_guess_new = kp_low+0.6*(kp_high-kp_low);
    kp = kp_guess;
		value_max_so_far = vn_comp;
	else
		if kp_guess < kp
			kp_low = kp_guess;
		else
			kp_high = kp_guess;
		end;
		kp_guess_new = kp_low+0.6*(kp_high-kp_low);
  end;
  kp_guess = kp_guess_new;
	kp = kp_guess;
	lp = mylabor_guess;
end;
	
