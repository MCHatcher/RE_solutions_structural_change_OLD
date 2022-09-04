% Application 3 - Pension reform
% Diamond (1965) model with CES utility; see Fedotenkov (2016,EL) and
% Hatcher(2019, EL) for log utility version of the model

%Calibration
alfa = 0.4;
betta = 1/(1+0.01)^35;
n = 0;
sigma = 0.90;  
tau = 0.20;
tau1 = 0.15;
betta_tild = betta/(1 + (1-alfa)*tau/alfa);
betta_tild1 = betta/(1 + (1-alfa)*tau1/alfa);

%Find steady state capital
Kmax = 0.15;
Kmin = 0.00001;
Nstep = 200000;
kstack = linspace(Kmin,Kmax,Nstep);

for j=1:Nstep

k(j) = kstack(j);
%Orig steady state
Resid(j) = abs( k(j)^(1-alfa)* ( (alfa + (1-alfa)*tau)*k(j)^(alfa-1)  + (alfa*betta)^(1/sigma)*k(j)^((alfa-1)/sigma) ) - (alfa*betta)^(1\sigma)*(1-tau)*(1-alfa)/(1+n) * k(j)^((alfa-1)/sigma) );
%Resid(j) = abs( (alfa + (1-alfa)*tau)*(k(j))^((1-alfa)/sigma) + (alfa*betta)^(1/sigma)*(k(j))^(1-alfa) - (alfa*betta)^(1\sigma)*(1-tau)*(1-alfa)/(1+n) );
%New steady state
Resid1(j) = abs( k(j)^(1-alfa)* ( (alfa + (1-alfa)*tau1)*k(j)^(alfa-1)  + (alfa*betta)^(1/sigma)*k(j)^((alfa-1)/sigma) ) - (alfa*betta)^(1\sigma)*(1-tau1)*(1-alfa)/(1+n) * k(j)^((alfa-1)/sigma) );
%Resid1(j) = abs( (alfa + (1-alfa)*tau1)*(k(j))^((1-alfa)/sigma) + (alfa*betta)^(1/sigma)*(k(j))^(1-alfa) - (alfa*betta)^(1\sigma)*(1-tau1)*(1-alfa)/(1+n) );
end

%Log utility - for comparison
klog = (betta_tild/(1+betta_tild)*(1-tau)*(1-alfa)/(1+n))^(1/(1-alfa));
klog1 = (betta_tild1/(1+betta_tild1)*(1-tau1)*(1-alfa)/(1+n))^(1/(1-alfa));

[MN,Ind] = min(Resid);
[MN1,Ind1] = min(Resid1);

%Orig steady state
kstar = k(Ind);
rstar = alfa*kstar^(alfa-1);
wstar = (1-alfa)*kstar^alfa;
cystar = (1-tau)*(1-alfa)*kstar^alfa - (1+n)*kstar;
costar = (alfa + (1-alfa)*tau)*kstar^alfa;

%New steady state
kstar1 = k(Ind1);
rstar1 = alfa*kstar1^(alfa-1);
wstar1 = (1-alfa)*kstar1^alfa;
cystar1 = (1-tau1)*(1-alfa)*kstar1^alfa - (1+n)*kstar1;
costar1 = (alfa + (1-alfa)*tau1)*kstar1^alfa;

%Reference regime
B1 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 1 0 (1-tau)*(-1)*(1-alfa)*kstar^alfa/cystar (1+n)*kstar/cystar 0; 0 0 0 0 1];
B2 = [0 -1/sigma 0 0 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
B3 = [0 0 0 0 0; 0 0 0 alfa-1 0; 0 0 0 alfa 0; 0 0 0 0 0; 0 0 0 alfa 0];
B4  = zeros(5,1);
B5 = zeros(5,1);

%Alternative (terminal) regime
xbar = [log(cystar); log(rstar); log(wstar); log(kstar); log(costar)];
xstar = [log(cystar1); log(rstar1); log(wstar1); log(kstar1); log(costar1)];
B1_tild = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 1 0 (1-tau1)*(-1)*(1-alfa)*kstar1^alfa/cystar1 (1+n)*kstar1/cystar1 0; 0 0 0 0 1];
B2_tild = [0 -1/sigma 0 0 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
B3_tild = [0 0 0 0 0; 0 0 0 alfa-1 0; 0 0 0 alfa 0; 0 0 0 0 0; 0 0 0 alfa 0];
B4_tild  = zeros(5,1);
B5_tild = (B1_tild - B2_tild - B3_tild)*(xstar-xbar);
