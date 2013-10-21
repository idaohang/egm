function u = util(Y,k1)
% utility
global tau; % parameters as globals
u = ((Y - k1).^(1-tau))/(1-tau);