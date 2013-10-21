function nul = resconst(z,k)
% resource constraint function to be solved
global alpha delta;
nul = exp(z)*k^alpha + (1-delta)*k;