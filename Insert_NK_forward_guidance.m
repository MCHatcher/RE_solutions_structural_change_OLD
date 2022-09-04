% Application 2 - NK Forward Guidance
% Model in Cagliarini and Kulish (2013)

%Calibration
alfa = 0.25; betta = 0.9925; psie = 0.1;
alfa_tild = alfa/(1+betta*alfa);
betta_tild = betta/(1+betta*alfa);
psie_tild = psie/(1+betta*alfa);

sigma = 1;
rho_i = 0.65;
theta_pi = 0.5; theta_y = 0.1; theta_dy = 0.2;

pistar = (0.05/4) / 2;  
i_ss = pistar - log(betta);
i_zlb = 0;

rho_a = 0.9; rho_g = 0.9; rho_mu = 0.9;
sigma_a = 0.007; sigma_g = 0.02; sigma_mu = 0.001;
VCOV = eye(3);

X_init = [pistar; 0; i_ss; 0; 0; 0];  %Set at SS as 
e_vec = zeros(length(VCOV),T_sim);  %e_vec = randn(length(VCOV),T_sim); %Stochastic simulation
e_vec(:,1) = [0; -4; 0];  %Impulse responses

%Reference regime
B1 = [1 -sigma*psie_tild 0 psie_tild 0 1/(1+betta*alfa); 0 1 1/sigma 0 -(1-rho_g)/sigma 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
B2 = [betta_tild 0 0 0 0 0; 1/sigma 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
B3 = [alfa_tild 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 rho_a 0 0; 0 0 0 0 rho_g 0; 0 0 0 0 0 rho_mu];
B4 = [0 0 0; 0 0 0; 0 0 0; sigma_a 0 0; 0 sigma_g 0; 0 0 sigma_mu];
B5 = [(1+betta*(alfa-1)-alfa)*pistar/(1+betta*alfa); -(1/sigma)*log(betta); i_zlb; 0; 0; 0];

%Alternative regime
B1_tild = [1 -sigma*psie_tild 0 psie_tild 0 1/(1+betta*alfa); 0 1 1/sigma 0 -(1-rho_g)/sigma 0; -theta_pi -(theta_y+theta_dy)  1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
B2_tild = [betta_tild 0 0 0 0 0; 1/sigma 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
B3_tild = [alfa_tild 0 0 0 0 0; 0 0 0 0 0 0; 0 -theta_dy  rho_i 0 0 0; 0 0 0 rho_a 0 0; 0 0 0 0 rho_g 0; 0 0 0 0 0 rho_mu];
B4_tild  = [0 0 0; 0 0 0; 0 0 0; sigma_a 0 0; 0 sigma_g 0; 0 0 sigma_mu];
B5_tild = [(1+betta*(alfa-1)-alfa)*pistar/(1+betta*alfa); -(1/sigma)*log(betta); (1-rho_i)*(pistar-log(betta))-theta_pi*pistar; 0; 0; 0];






