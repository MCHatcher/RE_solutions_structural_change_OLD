% Applications 2 and 3 - NK Forward Guidance 
% Model in Cagliarini and Kulish (2013)

%Calibration
alfa = 0.25; betta = 0.9925; psie = 0.1;
alfa_tild = alfa/(1+betta*alfa);
betta_tild = betta/(1+betta*alfa);
psie_tild = psie/(1+betta*alfa);

sigma = 1;
rho_i = 0.65;
theta_pi = 0.5; theta_y = 0.1; theta_dy = 0.2;
%theta_pi = 0.9*(1-rho_i);

pistar = (0.05/4) / 2;  
i_ss = pistar - log(betta);

rho_a = 0.9; rho_g = 0.9; rho_mu = 0.9;
sigma_a = 0.007; sigma_g = 0.02; sigma_mu = 0.001;

%Write model in the form in Binder and Pesaran (1997) and Cho and Moreno
%(2011), i.e. with exogenous processes not included in x(t)

%Terminal regime
B1_tild = [1 -sigma*psie_tild 0; 0 1 1/sigma; -theta_pi -(theta_y+theta_dy)  1 ];
B2_tild = [betta_tild 0 0; 1/sigma 1 0; 0 0 0];
B3_tild = [alfa_tild 0 0; 0 0 0; 0 -theta_dy rho_i];
B4_tild  = [-psie_tild -1/(1+betta*alfa) 0; 0 0 (1-rho_g)/sigma; 0 0 0];
B5_tild = [(1+betta*(alfa-1)-alfa)*pistar/(1+betta*alfa); -(1/sigma)*log(betta); (1-rho_i)*(pistar-log(betta))-theta_pi*pistar];






