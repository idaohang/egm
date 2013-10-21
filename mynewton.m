function l = mynewton(kp,k,z)
% newton iteration within value function iteration
%
% Created:
% 18.10.2013, Antti Ripatti
%
% based on the fortran code by Barillas and Fernandez Villaverde (2006)
%
global alpha delta theta;
lp = 0.3;
l = 0;
while abs(lp-l)<0.0000001 || lp>0.95 || lp<0.01
		l = lp;
		f_l  = kp-exp(z)*(k^alpha)*(l^(1-alpha))-(1-delta)*k+(theta/(1-theta))*(1-alpha)*exp(z)*(k^alpha)*(l^(-alpha))*(1-l);
		f_lp = -(1-alpha)*exp(z)*(k^alpha)*(l^(-alpha))+(-alpha)*(theta/(1-theta))*(1-alpha)*exp(z)*(k^alpha)*(l^(-1-alpha))*(1-l)-(theta/(1-theta))*(1-alpha)*exp(z)*(k^alpha)*(l^(-alpha));
		lp = l-(f_l/f_lp);
end;
l = min(0.95,lp);
l = max(0.01,l);
