%Code to check whether the Blanchard and Kahn (1980) condition for
%determinacy is satisfied by the terminal structure (Bi_tild, i=1,...,5) in
%the forward guidance example (Application 3)

clear; clc;

% Model and calibration
run Insert_NK_forward_guidance_BK

%Write model in the BK form: A*y(t) = B*E(t)y(t+1) + other(t)
%--> E(t)y(t+1) = *y(t) - B^(-1)*other(t), where C := B^(-1)*A
%Variable is predetermined if E(t)z(t+1) = z(t+1)

n = length(B1_tild);
B = [-B2_tild B1_tild; zeros(n) eye(n)];
A = [zeros(n) B3_tild; eye(n) zeros(n)];

%C = B \ A;
%if det(B) == 0
%    C = pinv(B)*A;
%end

%Degrees of indeterminacy

if det(B) ~= 0
    EIGEN = abs(eig(B\A));
end
    
if det(B) == 0
    [V,EIGEN] = eig(A,B);  %Generalized eigenvalues
end

EIGEN = abs(EIGEN);

%No. of non-predetermined vars
for i=1:n
    count(i) = 0;
    if nnz(B2_tild(:,i) > 0) || nnz(B1_tild(i,1:i-1) ~= 0)  || nnz(B1_tild(i,i+1:n) ~= 0)
        count(i) = 1;
    end
end
    
nf = sum(count)
k = nf - nnz(EIGEN>1) 

disp('No. of degrees of indeterminacy is')
k = nf - nnz(EIGEN>1) 

