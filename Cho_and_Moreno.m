% Based on Cho and Moreno (2011, Journal of Econ Dynamics and Control)
% Recursive algorithm for computing the unique fixed-structure fundamental
% rational expectations solution that satisfies the no-bubbles condition
% Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

%Cho and Moreno (2011) method for computing fixed structure solutions
J = 10000;  %no. iterations

%% Warning: DO NOT AMEND. Structures are specified in the main files.

if abs(det(B1)) >  0 && abs(det(B1_tild)) > 0
     
    A = B1 \ B2;  B = B1 \ B3;  C = B1 \ B4;  D = B1 \ B5;
    A_tild = B1_tild \ B2_tild;  B_tild = B1_tild \ B3_tild;  
    C_tild = B1_tild \ B4_tild;  D_tild = B1_tild \ B5_tild;

%Note R = 0_{m \times m} as shocks in e(t) are white noise. Persistent shocks can be included in x(t).
%See Kulish and Pagan (2017, Journal of Applied Econometrics)
    
%Taken from main file
I = eye(length(B1));

%Initial values (run main file with structure)
AT = A; BT = B; CT = C; DT = D;
AT_tild = A_tild; BT_tild = B_tild; CT_tild = C_tild; DT_tild = D_tild;

for i=1:J
    
    decomp = decomposition(I - A*BT);
    if isIllConditioned(decomp)==1
        break
    end
    
    BT = (I - A*BT) \ B;     
    AT = (I - A*BT) \  A*AT;
    CT = (I - A*BT) \ C;            %Since R = 0_{ m \times m}
    DT = (I - A*BT) \ (D + A*DT);
    
end

for i = 1:J
    
    decomp = decomposition(I - A_tild*BT_tild);
    if isIllConditioned(decomp)==1
        break
    end
    
    BT_tild = (I - A_tild*BT_tild) \ B_tild;     
    AT_tild = (I - A_tild*BT_tild) \  A_tild*AT_tild;
    CT_tild = (I - A_tild*BT_tild) \ C_tild;            %Since R = 0_{ m \times m}
    DT_tild = (I - A_tild*BT_tild) \ (D_tild + A_tild*DT_tild);
      
end

%Check no-bubbles condition 
Check1 = AT, Check2 = AT_tild 

%Solution matrices
Omega_bar = BT; Gama_bar = CT; Psi_bar = DT;
Omega_tild = BT_tild; Gama_tild = CT_tild; Psi_tild = DT_tild;

end

if abs(det(B1)) ==  0    ||   abs(det(B1_tild)) == 0  
    disp('Matrix A and/or A_tild non-invertible')
    run Cho_Moreno_alternative
end
