% Arbitrary announcement date T_ann, applied to Application 3 - Pension reform
% Diamond (1965) model with CES utility; see Fedotenkov (2016,Econ Lett) and
% Hatcher(2019, Econ Lett) for log utility version of the model
%Model structures are defined in the 'Insert' files
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

clc; clear;

% Announcement date and final date before terminal structure
T_ann = 4; T_tild = 4;
T_sim = 80; %Simulation length

% Model and calibration
run Insert_pension_reform

% Fixed structure solutions (Cho and Moreno 2011, JEDC)
run Cho_and_Moreno

%Indicator variable    
ind = ones(T_sim,1); ind(T_tild+1:T_sim,1) = 0;  %permanent structural change
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
    
    if t <= T_ann-1
        Omeg = Omega_bar;
        Gama = Gama_bar;
        Psi = Psi_bar;
    end
                      
    Omeg_t(:,:,t) = Omeg;
    Gama_t(:,:,t) = Gama;
    Psi_t(:,:,t) = Psi;
    
 end

%Initial values
X_init = [0; 0; 0; 0; 0]; 

%Prepare for simulations
X = X_init;
 
%For comparison
X_u = X_init;  %Unannounced reform

%Simulation results        
for t=1:T_sim 
    
        if t < T_ann
            X = Omega_bar*X + Psi_bar;
        end
        
        if t >= T_ann
            X = Omeg_t(:,:,t)*X + Psi_t(:,:,t);
        end
        
        %Store for later
        X_stack(:,t) = X;
       
        
        if t <= T_ann
            X_u = Omega_bar*X_u + Psi_bar;
        end
        
        if t > T_ann
            X_u = Omega_tild*X_u + Psi_tild;
        end
        
        Xu_stack(:,t) = X_u;
    
        Periods(t) = t-T_ann;  %To plot from period 0
    
end 

Pension_reform_plotter



