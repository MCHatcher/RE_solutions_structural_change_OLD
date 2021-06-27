%Arbitrary announcement date T_ann, applied to:
%Ireland (2007) NK Model, replicates Cagliarini and Kulish (2013, Fig 3)
%To study a different example, simply change the parameters and matrices
%A,B,C,D, A_tild, B_tild, C_tild, D_tild, R in the 'Insert' files

clc; clear; close all;

% Announcement date and final date before terminal structure
T_ann = 4; T_tild = 7;
T = 16; %Simulation length

% Model and calibration
%run Insert_NK_forward_guidance
run Insert_NK_inflation_target

% Fixed structure solutions (Cho and Moreno 2011, JEDC)
run Cho_and_Moreno

Check1 = AT; Check2 = AT_tild;
% Solution matrices
Omega_bar = BT; Gama_bar = CT; Psi_bar = DT;
Omega_bar1 = BT_tild; Gama_bar1 = CT_tild; Psi_bar1 = DT_tild;

%Indicator variable    
%To pick a random structure use: ind = randi([0,1],T_tild,1)
ind = zeros(T_tild,1); ind(1:T_tild,1) = 1;

%Computation of matrix recursion
 for j=1:T_tild       
            
    Aj = ind(T_tild+1-j,1)*A + (1-ind(T_tild+1-j,1))*A_tild;
    Bj = ind(T_tild+1-j,1)*B + (1-ind(T_tild+1-j,1))*B_tild;
    Cj = ind(T_tild+1-j,1)*C + (1-ind(T_tild+1-j,1))*C_tild;
    Dj = ind(T_tild+1-j,1)*D + (1-ind(T_tild+1-j,1))*D_tild;
        
    if j == 1
         
        Omeg = (I - Aj*Omega_bar1) \ Bj; 
        Gama = (I - Aj*Omega_bar1) \ (Aj*Gama_bar1*R + Cj); 
        Psi = (I - Aj*Omega_bar1) \ (Aj*Psi_bar1 + Dj);
        
        %%If reference regime is terminal regime
        %Omeg = (I - Aj*Omega_bar) \ Bj; Gama = (I - Aj*Omega_bar) \ (Aj*Gama_bar*R + Cj); 
        %Psi = (I - Aj*Omega_bar) \ (Aj*Psi_bar + Dj)
        
    end    
        
    if j > 1 && j <= T_tild+1-T_ann
    
        Omeg = (I - Aj*Omeg) \ Bj;
        Gama = (I - Aj*Omeg) \ (Aj*Gama*R + Cj);
        Psi = (I - Aj*Omeg) \ (Aj*Psi + Dj);
    
    end
      
    if j > T_tild+1-T_ann       
        Omeg = Omega_bar; Gama = Gama_bar; Psi = Psi_bar; 
    end
              
    Omeg_j(:,:,j) = Omeg;
    Gama_j(:,:,j) = Gama;
    Psi_j(:,:,j) = Psi;
    
 end

% Prepare for simulations 
pi = zeros(T,1); y = pi; int = pi;
rng(10); %Set random seed to replicate results
eps_a = sigma_a*randn(T,1);
eps_z = sigma_z*randn(T,1);
eps_e = sigma_e*randn(T,1);
 
%Initial values
Z_init(1:length(R),1) = 0; Z = Z_init;
X_init = [pistar; 0; i_ss]; X = X_init; 
eps_z(1) = 0; eps_e(1) = 0; eps_a(1) = 0.02;        
        
% Original structure
X_orig = X_init; 
% Final structure
X_init2 = [pistar_new; 0; i_ss_new]; 
X_fin = X_init2;

%Deterministic case
eps_z(2:T+1) = 0; eps_e(2:T+1) = 0; eps_a(2:T+1) = 0;

%Simulation results        
for t=1:T 
         
        Z = R*Z + [eps_z(t); eps_e(t); eps_a(t)];
        
        if t <= T_tild
            X = Omeg_j(:,:,T_tild+1-t)*X + Gama_j(:,:,T_tild+1-t)*Z + Psi_j(:,:,T_tild+1-t);
        end
        
        %Terminal structure
        if t > T_tild 
            X = Omega_bar1*X + Gama_bar1*Z + Psi_bar1;
            %%If terminal structure reference regime
            %X = Omega_bar*X + Gama_bar*Z + Psi_bar;
        end
        
        pi(t) = X(1); y(t) = X(2); int(t) = X(3);
        
        X_orig = Omega_bar*X_orig + Gama_bar*Z + Psi_bar;
        pi_orig(t) = X_orig(1); y_orig(t) = X_orig(2); 
        int_orig(t) = X_orig(3);
        
        X_fin = Omega_bar1*X_fin + Gama_bar1*Z + Psi_bar1;
        pi_fin(t) = X_fin(1); y_fin(t) = X_fin(2);
        int_fin(t) = X_fin(3);
    
        Periods(t) = t;
    
end       

hold on, 
subplot(1,3,1), plot(Periods, 4*100*pi, 'k'), hold on, plot(Periods, 4*100*pi_orig,'--k'), 
hold on, plot(Periods, 4*100*pi_fin,'g'), title('Inflation'), 
axis([-inf,inf,2,6]), xlabel('Periods'), ylabel('% (annualized)')
subplot(1,3,2), plot(Periods, 100*y,'k'), hold on, plot(Periods, 100*y_orig, '--k'), title('Output'),
xlabel('Periods'), ylabel('%')
subplot(1,3,3), plot(Periods, 4*100*int,'k'), hold on, plot(Periods, 4*100*int_orig,'--k'), axis([-inf,inf,-0.002,0.006]) 
hold on, plot(Periods, 4*100*int_fin,'g'), title('Nominal interest rate'), axis([-inf,inf,5,9.01]),
xlabel('Periods'), ylabel('% (annualized)')
