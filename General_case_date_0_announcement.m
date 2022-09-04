%Announcement in the first period (t=1 in Matlab), forward guidance (baseline) 
%Ireland (2007) NK Model, replicates Cagliarini and Kulish (2013, Fig 3)
%To study a different example, simply change the parameters and matrices
%Model structures are defined in the 'Insert' files
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

clc; clear; %close all;

T_tild = 6; % Final date before terminal structure
T = 3; % Start of implementation period
T_sim = 15; % Simulation length

% Model and calibration
%run Insert_NK_inflation_target
run Insert_NK_forward_guidance

% Fixed structure solutions (Cho and Moreno 2011, JEDC)
run Cho_and_Moreno

%Indicator variable    
%ind = ones(T_sim,1); ind(T_tild+1:T_sim,1) = 0;  %permanent structural change
ind = zeros(T_sim,1); ind(T:T_tild) = 1; % temporary structural change

%Initial values for recursion
Omeg = Omega_tild; Gama = Gama_tild; Psi = Psi_tild;  

%Computation of matrix recursion
 for t=T_sim:-1:1       
            
    B1t = ind(t)*B1 + (1-ind(t))*B1_tild;
    B2t = ind(t)*B2 + (1-ind(t))*B2_tild;
    B3t = ind(t)*B3 + (1-ind(t))*B3_tild;
    B4t = ind(t)*B4 + (1-ind(t))*B4_tild;
    B5t = ind(t)*B5 + (1-ind(t))*B5_tild;
         
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
