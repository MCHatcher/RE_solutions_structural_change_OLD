%Date t announcement applied to:
%Ireland (2007) NK Model. Replicates a result in Fig 3 of Hatcher (2021).
%To study a different example, simply change the parameters and matrices
%A,B,C,D, A_tild, B_tild, C_tild, D_tild, R in the 'Insert' files

clc; clear; close all;

% Final date before terminal structure
T_tild = 4;
T = 12; %Simulation length

% Model and calibration
run Insert_NK_forward_guidance
%run Insert_NK_inflation_target

% Fixed structure solutions (Cho and Moreno 2011, JEDC)
run Cho_and_Moreno

% Solution matrices
Omega_bar = BT; Gama_bar = CT; Psi_bar = DT;
Omega_bar1 = BT_tild; Gama_bar1 = CT_tild; Psi_bar1 = DT_tild;

%Indicator variable    
%To pick a random structure use: ind = randi([0,1],T_tild,1)
ind = zeros(T_tild,1); ind(1:T_tild-2,1) = 1;

%Computation of matrix recursion
 for j=1:T_tild       
            
    Aj = ind(T_tild+1-j,1)*A + (1-ind(T_tild+1-j,1))*A_tild;
    Bj = ind(T_tild+1-j,1)*B + (1-ind(T_tild+1-j,1))*B_tild;
    Cj = ind(T_tild+1-j,1)*C + (1-ind(T_tild+1-j,1))*C_tild;
    Dj = ind(T_tild+1-j,1)*D + (1-ind(T_tild+1-j,1))*D_tild;
        
    if j == 1
         
        % If alternative regime is terminal regime
        %Omeg = (I - Aj*Omega_bar1) \ Bj; 
        %Gama = (I - Aj*Omega_bar1) \ (Aj*Gama_bar1*R + Cj); 
        %Psi = (I - Aj*Omega_bar1) \ (Aj*Psi_bar1 + Dj);
        
        %%If reference regime is terminal regime
        Omeg = (I - Aj*Omega_bar) \ Bj;
        Gama = (I - Aj*Omega_bar) \ (Aj*Gama_bar*R + Cj); 
        Psi = (I - Aj*Omega_bar) \ (Aj*Psi_bar + Dj);
        
    end    
        
    if j > 1 
    
        Omeg = (I - Aj*Omeg) \ Bj;
        Gama = (I - Aj*Omeg) \ (Aj*Gama*R + Cj);
        Psi = (I - Aj*Omeg) \ (Aj*Psi + Dj);
    
    end
              
    Omeg_j(:,:,j) = Omeg;
    Gama_j(:,:,j) = Gama;
    Psi_j(:,:,j) = Psi;
    
 end

%Prepare for simulations 
pi = zeros(T,1); y = pi; int = pi;
rng(10); %Set random seed to replicate results
eps_a = sigma_a*randn(T,1);
eps_z = sigma_z*randn(T,1);
eps_e = sigma_e*randn(T,1);
 
%Initial values
Z_init(1:length(R),1) = 0; Z = Z_init;
X_init = [pistar; 0; i_ss]; X = X_init;
%Deterministic case
eps_z(2:T) = 0; eps_e(2:T) = 0; eps_a(2:T) = 0;
eps_a(1) = -0.125;
        
% Original structure
X_orig = X_init; 
% Final structure
X_init2 = [pistar; 0; i_ss]; 
X_fin = X_init2;

%Simulation results        
for t=1:T 
         
        Z = R*Z + [eps_z(t); eps_e(t); eps_a(t)];
        
        if t <= T_tild
            X = Omeg_j(:,:,T_tild+1-t)*X + Gama_j(:,:,T_tild+1-t)*Z + Psi_j(:,:,T_tild+1-t);
        end
        
        %Terminal structure
        if t > T_tild 
            %X = Omega_bar1*X + Gama_bar1*Z + Psi_bar1;
            %%If terminal structure reference regime
            X = Omega_bar*X + Gama_bar*Z + Psi_bar;
        end
        
        pi(t) = X(1); y(t) = X(2); int(t) = X(3);
    
        Periods(t) = t-1;
    
end       

hold on, 
subplot(1,3,1), plot(Periods, 4*100*pi, 'g'), hold on,  title('Inflation'), xlabel('Periods'), ylabel('% (annualized)')
subplot(1,3,2), plot(Periods, 100*y,'g'), hold on, title('Output'), xlabel('Periods'), ylabel('%')
subplot(1,3,3), plot(Periods, 4*100*int,'g'), hold on, title('Nominal interest rate'), 
xlabel('Periods'), ylabel('% (annualized)')
