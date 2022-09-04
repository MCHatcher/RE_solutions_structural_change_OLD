%Arbitrary announcement date T_ann, imperfect credibility applied to:
%Ireland (2007) NK Model, replicates Cagliarini and Kulish (2013, Fig 3)
%To study a different example, simply change the parameters and matrices
%Model structures are defined in the 'Insert' files
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

clc; clear; %close all;

% Announcement date and final date before terminal structure
T_ann = 4; T = 7; T_tild = 10;
T_sim = 16; %Simulation length

% Model and calibration
%run Insert_NK_forward_guidance
%run Insert_NK_inflation_target
run Insert_NK_inflation_target_past_cred

% Fixed structure solutions (Cho and Moreno 2011, JEDC)
run Cho_and_Moreno

%Indicator variable    
ind = ones(T_sim,1); ind(T+1:T_sim,1) = 0;  %permanent structural change
%Initial values for recursion
Omeg = Omega_tild; Gama = Gama_tild; Psi = Psi_tild;  

%Computation of matrix recursion
 for t=T_sim:-1:1       
            
    B1t = ind(t)*B1 + (1-ind(t))*B1_tild;
    B2t = ind(t)*B2 + (1-ind(t))*B2_tild;
    B3t = ind(t)*B3 + (1-ind(t))*B3_tild;
    B4t = ind(t)*B4 + (1-ind(t))*B4_tild;
    B5t = ind(t)*B5 + (1-ind(t))*B5_tild;
    
    %for solution refinements
    
    Lambda = [0.70 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
    %Lambda = eye(6);  %for perfect credibility case
    
    if t >= T_tild+1
        Lambda = eye(6);
    end
    
    B1t_hat = B1t - B2t*(eye(length(B1t))-Lambda)*F0;
    B2t_hat = B2t*Lambda;
    B3t_hat = B3t + B2t*(eye(length(B1t))-Lambda)*F1;
    B4t_hat = B4t + B2t*(eye(length(B1t))-Lambda)*F2;
    B5t_hat = B5t + B2t*(eye(length(B1t))-Lambda)*F3;
         
   Omeg = (B1t_hat - B2t_hat*Omeg) \ B3t_hat; 
   Gama = (B1t_hat - B2t_hat*Omeg) \ B4t_hat; 
   Psi = (B1t_hat - B2t_hat*Omeg) \ (B2t_hat*Psi + B5t_hat);    
        
    if t >= T_tild+1
        Omeg = Omega_tild;
        Gama = Gama_tild;
        Psi = Psi_tild;   
    end
    
    if t <= T_ann-1
        Omeg = Omega_bar;
        Gama = Gama_bar;
        Psi = Psi_bar;
    end
                      
    Omeg_t(:,:,t) = Omeg;
    Gama_t(:,:,t) = Gama;
    Psi_t(:,:,t) = Psi;
    
 end
 
 %Prepare for simulations
 X = X_init;
 
 %For comparison
 X_orig = X; X_fin = X_init_new; 

%Simulation results        
for t=1:T_sim-1 
        
        X = Omeg_t(:,:,t)*X + Gama_t(:,:,t)*e_vec(:,t) + Psi_t(:,:,t);
        Xe = Omeg_t(:,:,t+1)*X + Psi_t(:,:,t+1);
        
        %Store for later
        X_stack(:,t) = X;  Xe_stack(:,t) = Xe;
        
        %Under original structure
        X_orig = Omega_bar*X_orig + Gama_bar*e_vec(:,t) + Psi_bar;
        X_stack_orig(:,t) = X_orig;
        
        %Under terminal structure
        X_fin = Omega_tild*X_fin + Gama_tild*e_vec(:,t) + Psi_tild;
        X_stack_fin(:,t) = X_fin;
    
        Periods(t) = t;
    
end 

NK_inflation_target_plotter_cred
