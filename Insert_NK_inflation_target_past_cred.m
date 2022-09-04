% Application 1 - NK Shift in Inflation Target (imperfect credibility)
% Model in Cagliarini and Kulish (2013)

%Calibration
alfa = 0.25; betta = 0.9925; psie = 0.1;
alfa_tild = alfa/(1+betta*alfa);
betta_tild = betta/(1+betta*alfa);
psie_tild = psie/(1+betta*alfa);

sigma = 1;
theta_pi = 0.5; theta_y = 0.1; theta_dy = 0.2;
rho_i = 0.65;

pistar = 0.05/4;  pistar_new = pistar/2;
i_ss = pistar - log(betta); i_ss_new = pistar_new - log(betta);

rho_a = 0.9; rho_g = 0.9; rho_mu = 0.9;
sigma_a = 0.007; sigma_g = 0.02; sigma_mu = 0.001;
VCOV = eye(3);

X_init = [pistar; 0; i_ss; 0; 0; 0];  %Set at SS as default
X_init_new = [pistar_new; 0; i_ss_new; 0; 0; 0];  %Set at SS as default
%e_vec = randn(length(VCOV),T);  %Stochastic simulation 
e_vec(:,1) = [0; 1; 0];  e_vec(:,2:T_sim) = zeros(3,T_sim-1); %Impulse responses

%Imperfect credibility matrices
Dummy = 0;   %Set to 1 (0) for expectations that track past inflation (original inflation target) 
F0 = zeros(length(X_init),length(X_init));
F1 = [Dummy 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
F2 = zeros(length(X_init),length(VCOV));
F3 = [(1-Dummy)*pistar; 0; 0; 0; 0; 0];

%Reference regime
B1 = [1 -sigma*psie_tild 0 psie_tild 0 1/(1+betta*alfa); 0 1 1/sigma 0 -(1-rho_g)/sigma 0; -theta_pi -(theta_y+theta_dy)  1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
B2 = [betta_tild 0 0 0 0 0; 1/sigma 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
B3 = [alfa_tild 0 0 0 0 0; 0 0 0 0 0 0; 0 -theta_dy  rho_i 0 0 0; 0 0 0 rho_a 0 0; 0 0 0 0 rho_g 0; 0 0 0 0 0 rho_mu];
B4  = [0 0 0; 0 0 0; 0 0 0; sigma_a 0 0; 0 sigma_g 0; 0 0 sigma_mu];
B5 = [(1+betta*(alfa-1)-alfa)*pistar/(1+betta*alfa); -(1/sigma)*log(betta); (1-rho_i)*(pistar-log(betta))-theta_pi*pistar; 0; 0; 0];

%Alternative regime
B1_tild = B1;
B2_tild = B2;
B3_tild = B3;
B4_tild = B4;
B5_tild = [(1+betta*(alfa-1)-alfa)*pistar_new/(1+betta*alfa); -(1/sigma)*log(betta); (1-rho_i)*(pistar_new-log(betta))-theta_pi*pistar_new; 0; 0; 0];



