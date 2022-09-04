%Announcement in the first period, Type 2 Credibility, forward guidance 
%Ireland (2007) NK Model, replicates Cagliarini and Kulish (2013, Fig 3)
%To study a different example, simply change the parameters and matrices
%Model structures are defined in the 'Insert' files
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

clc; clear; %close all;

T = 3; % Start of implementation period
T1 = 6; % Final period before terminal structure (B_tild)
T_tild = 6; % Final period of imperfect credibility
T_sim = 13; % Simulation length

% Model and calibration
%run Insert_NK_inflation_target
run Insert_NK_forward_guidance

% Fixed structure solutions (Cho and Moreno 2011, JEDC)
run Cho_and_Moreno

%Indicator variable    
%ind = ones(T_sim,1); ind(T_tild+1:T_sim,1) = 0;  %permanent structural change
ind = zeros(T_sim,1); ind(T:T_tild) = 1; % temporary structural change
ind_IC = zeros(T_sim,1); ind(T:T_tild-2) = 1;  %alternative structure under imperfect credibility

%Credibility parameter
p = 0.5;

%Initial values for recursions
Omeg_star = Omega_tild; Gama_star = Gama_tild; Psi_star = Psi_tild;
Omeg_IC = Omega_tild; Gama_IC = Gama_tild; Psi_IC = Psi_tild;

%Computation of matrix recursion
 for t=T_sim:-1:1     
     
    B1t = ind(t)*B1 + (1-ind(t))*B1_tild;
    B2t = ind(t)*B2 + (1-ind(t))*B2_tild;
    B3t = ind(t)*B3 + (1-ind(t))*B3_tild;
    B4t = ind(t)*B4 + (1-ind(t))*B4_tild;
    B5t = ind(t)*B5 + (1-ind(t))*B5_tild;
            
    B1t_IC = ind_IC(t)*B1 + (1-ind_IC(t))*B1_tild;
    B2t_IC = ind_IC(t)*B2 + (1-ind_IC(t))*B2_tild;
    B3t_IC = ind_IC(t)*B3 + (1-ind_IC(t))*B3_tild;
    B4t_IC = ind_IC(t)*B4 + (1-ind_IC(t))*B4_tild;
    B5t_IC = ind_IC(t)*B5 + (1-ind_IC(t))*B5_tild;
         
   Omeg_star = (B1t - B2t*Omeg_star) \ B3t; 
   Gama_star = (B1t - B2t*Omeg_star) \ B4t; 
   Psi_star = (B1t - B2t*Omeg_star) \ (B2t*Psi_star + B5t); 
   
   Omeg_IC = (B1t_IC - B2t_IC*Omeg_IC) \ B3t_IC; 
   Gama_IC  = (B1t_IC - B2t_IC*Omeg_IC) \ B4t_IC; 
   Psi_IC  = (B1t_IC - B2t_IC*Omeg_IC) \ (B2t_IC*Psi_IC + B5t_IC);
        
    if t >= T_tild+1
        Omeg_star = Omega_tild;
        Gama_star = Gama_tild;
        Psi_star = Psi_tild;   
        
        Omeg_IC = Omega_tild; 
        Gama_IC = Gama_tild; 
        Psi_IC = Psi_tild;
    end
                      
    Omeg_star_t(:,:,t) = Omeg_star;
    Gama_star_t(:,:,t) = Gama_star;
    Psi_star_t(:,:,t) = Psi_star;
    
    Omeg_IC_t(:,:,t) = Omeg_IC;
    Gama_IC_t(:,:,t) = Gama_IC;
    Psi_IC_t(:,:,t) = Psi_IC;
    
    Omeg_tild_IC_t(:,:,t) = p*Omeg_star_t(:,:,t) + (1-p)*Omeg_IC_t(:,:,t);
    Gama_tild_IC_t(:,:,t) = p*Gama_star_t(:,:,t) + (1-p)*Gama_IC_t(:,:,t);
    Psi_tild_IC_t(:,:,t) = p*Psi_star_t(:,:,t) + (1-p)*Psi_IC_t(:,:,t);
        
 end

 %Initial values for final recursion
Omeg = Omega_tild; Gama = Gama_tild; Psi = Psi_tild;
Lambda = 0.7;

%Computation of matrix recursion
 for t=T_sim:-1:1       
            
    B1t = ind(t)*B1 + (1-ind(t))*B1_tild;
    B2t = ind(t)*B2 + (1-ind(t))*B2_tild;
    B3t = ind(t)*B3 + (1-ind(t))*B3_tild;
    B4t = ind(t)*B4 + (1-ind(t))*B4_tild;
    B5t = ind(t)*B5 + (1-ind(t))*B5_tild;
    
    if t <= T_tild
        B1t = B1t - B2t*(1-Lambda)*Omeg_tild_IC_t(:,:,t+1);
        B2t = B2t*Lambda;
        B5t = B5t + B2t*(1-Lambda)*Psi_tild_IC_t(:,:,t+1);
    end
         
   Omeg = (B1t - B2t*Omeg) \ B3t; 
   Gama = (B1t - B2t*Omeg) \ B4t; 
   Psi = (B1t - B2t*Omeg) \ (B2t*Psi + B5t);    
        
    if t >= T_tild+1
        Omeg = Omega_tild;
        Gama = Gama_tild;
        Psi = Psi_tild;   
    end
                      
    Omeg_t(:,:,t) = Omeg;
    Gama_t(:,:,t) = Gama;
    Psi_t(:,:,t) = Psi;
    
 end
 
 %Prepare for simulations
 X = X_init;
 
 %For comparison
 X_fin = X_init;

%Simulation results        
for t=1:T_sim 
        
        X = Omeg_t(:,:,t)*X + Gama_t(:,:,t)*e_vec(:,t) + Psi_t(:,:,t);
        
        %Store for later
        X_stack(:,t) = X;
        
        %Under terminal structure (for comparison if wanted)
        X_fin = Omega_tild*X_fin + Gama_tild*e_vec(:,t) + Psi_tild;
        X_stack_fin(:,t) = X_fin;
    
        Periods(t) = t-1;  %To plot from period 0
    
end 

%NK_inflation_target_plotter
NK_forward_guidance_plotter
