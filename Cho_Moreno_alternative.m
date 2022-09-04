% Alternative algorithm for computing the unique fixed-structure fundamental rational 
% expectations solution using time iteration (along similar lines to Cho and Moreno (2011, JEDC))
% Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

J = 10000;  %no. iterations

%% Warning: DO NOT AMEND. Structures are specified in the main files.

%Taken from main file
I = eye(length(B1));

%Initial values (run main file with structure)
OmegT = B3;  GamaT = B4;  PsiT = B5;
OmegT_tild = B3_tild; GamaT_tild = B4_tild;  PsiT_tild = B5_tild;

for i=1:J
    
    OmegT = (B1 - B2*OmegT) \ B3;
    GamaT = (B1 - B2*OmegT) \ B4;
    PsiT =  (B1 - B2*OmegT) \ (B2*PsiT + B5);
    
    OmegT_tild = (B1_tild - B2_tild*OmegT_tild) \ B3_tild;
    GamaT_tild = (B1_tild - B2_tild*OmegT_tild) \ B4_tild;
    PsiT_tild =  (B1_tild - B2_tild*OmegT_tild) \ (B2_tild*PsiT_tild + B5_tild);
      
end

%Solution matrices
Omega_bar = OmegT; Gama_bar = GamaT; Psi_bar = PsiT;
Omega_tild = OmegT_tild; Gama_tild = GamaT_tild; Psi_tild = PsiT_tild;