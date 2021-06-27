% Application 2 - NK Forward Guidance

%Calibration
alfa = 0.25; betta = 0.9925; psi = 0.1;
alfa_tild = alfa/(1+betta*alfa);
betta_tild = betta/(1+betta*alfa);
psi_tild = psi/(1+betta*alfa);

sigma = 1; theta_pi = 0.5;
theta_y = 0.1; theta_g = 0.2;
rho_i = 0.65;

pistar = 0.05/4;
i_ss = pistar - log(betta);
i_star = 0.0001;

rho_a = 0.9; rho_z = 0.9; rho_e = 0.9;
sigma_a = 0; sigma_z = sigma_a; sigma_e = sigma_a;
R = [rho_z 0 0; 0 rho_e 0; 0 0 rho_a];

%Reference regime
B1 = [1 -sigma*psi_tild 0; 0 1 1/sigma; -theta_pi -(theta_y+theta_g)  1 ];
B2 = [betta_tild 0 0; 1/sigma 1 0; 0 0 0];
B3 = [alfa_tild 0 0; 0 0 0; 0 -theta_g rho_i];
B4  = [-psi_tild -1/(1+betta*alfa) 0; 0 0 (1-rho_a)/sigma; 0 0 0];
B5 = [(1+betta*(alfa-1)-alfa)*pistar/(1+betta*alfa); -(1/sigma)*log(betta); (1-rho_i)*(pistar-log(betta))-theta_pi*pistar];

A = B1 \ B2; B = B1 \ B3;
C = B1 \ B4; D = B1 \ B5;

%Alternative regime
B1_tild = [1 -sigma*psi_tild 0; 0 1 1/sigma; 0 0  1 ];
B2_tild = [betta_tild 0 0; 1/sigma 1 0; 0 0 0];
B3_tild = [alfa_tild 0 0; 0 0 0; 0 0 0];
B4_tild  = [-psi_tild -1/(1+betta*alfa) 0; 0 0 (1-rho_a)/sigma; 0 0 0];
B5_tild = [(1+betta*(alfa-1)-alfa)*pistar/(1+betta*alfa); -(1/sigma)*log(betta); i_star];

A_tild = B1_tild \ B2_tild; B_tild = B1_tild \ B3_tild;
C_tild = B1_tild \ B4_tild; D_tild = B1_tild \ B5_tild;
