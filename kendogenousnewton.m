function kend =  kendogenousnewton(z,labor,Cstar,kprime)
% kendogenousnewton Newton step when k is endogenous?
%
% Created:
% 21.10.2013, Antti Ripatti
%
% based on the fortran code by Barillas and Fernandez Villaverde (2006)
%
kend = kprime;
kp = 1e5;
while abs(kp-kend)<0.0000001
  f_k  = kprime + Cstar-exp(z)*kend^alpha*labor^(1-alpha)-(1-delta)*kend;
	f_kp = -alpha*exp(z)*kend^(alpha-1)*labor^(1-alpha)-(1-delta);
	kp = kend-(f_k/f_kp);
	kend = kp;
end;
