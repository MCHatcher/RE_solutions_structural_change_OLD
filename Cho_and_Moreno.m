% Based on Cho and Moreno (2011, Journal of Econ Dynamics and Control)
% Recursive algorithm for computing the unique fixed-structure fundamental
% rational expectations solution that satisfies the no-bubbles condition
% Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

%Cho and Moreno (2011) method for computing fixed structure solutions
J = 10000;  %no. iterations

%% Warning: DO NOT AMEND. Structures are specified in the main files.

%Taken from main file
I = eye(length(A));

%Initial values (run main file with structure)
AT = A; BT = B; CT = C; DT = D;
AT_tild = A_tild; BT_tild = B_tild;
CT_tild = C_tild; DT_tild = D_tild;

for i=1:J
    
    BT = (I - A*BT) \ B;     
    AT = (I - A*BT) \  A*AT;
    CT = (I - A*BT) \ (C + A*CT*R);
    DT = (I - A*BT) \ (D + A*DT);
    
    BT_tild = (I - A_tild*BT_tild) \ B_tild;     
    AT_tild = (I - A_tild*BT_tild) \  A_tild*AT_tild;
    CT_tild = (I - A_tild*BT_tild) \ (C_tild + A_tild*CT_tild*R);
    DT_tild = (I - A_tild*BT_tild) \ (D_tild + A_tild*DT_tild);
      
end

Check1 = AT, Check2 = AT_tild  %Can use these to check no-bubbles condition
%Solution matrices
Omega_bar = BT; Gama_bar = CT; Psi_bar = DT;
Omega_bar1 = BT_tild; Gama_bar1 = CT_tild; Psi_bar1 = DT_tild;
