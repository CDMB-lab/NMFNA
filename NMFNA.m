% Algorithm for NMFNA
function [S11,S22,G1,G2,res,errorx] = NMFNA(X1, X2, X3,W1,W2, K, maxiter,lamda1,lamda2)

% Input argument:
% X1 (N1,N1): N (dimensionallity) x M1 (dimensionallity) non negative input matrix
% X2 (N2,N2): N (dimensionallity) x M2 (dimensionallity) non negative input matrix
% X3 (N1,N2): N (dimensionallity) x M3 (dimensionallity) non negative input matrix
% K - reduced dimension of matrix factorization, that is the number of identified modules;
% maxiter - the maximum number of iterations;
% lamda1 and lamda2: graph regularized constraint parameters 

% Output:
% G1: N1 x K matrix
% G2: N2 x K matrix
% S11: K x K matrix
% S22: K x K matrix

% input matrix
R11=X1;
R12=X3;
R22=X2;

%size of matrix
[nr1,nc1] = size(R11);
[nr,nc] = size(R12);
[nr2,nc2] = size(R22);

DCol1 = full(sum(W1,2));
D1 = spdiags(DCol1,0,nc1,nc1);
DCol2 = full(sum(W2,2));
D2 = spdiags(DCol2,0,nc2,nc2);

if((nr1 ~= nc1) || (nr2 ~= nc2) || (nr ~= nr1) || (nc ~= nr2))
    error('The input data matrices can not be applied to the algorithm!')
end

n1 = nr1;
n2 = nr2;

k = K;
if((k >= n1) || (k >= n2))
    error('The number of clusters should be set smaller!')
end

tol =10^(-5);

% Initialization of G1,G2,S11,S22.
[U1,V11,~]=svds(R11,k);
[U2,V22,~]=svds(R22,k);

G1=abs(U1);
G2=abs(U2);

S11=abs(V11);
S22=abs(V22);

% Initialization of parameters
lambda12 = n1/n2;
lambda22 = (n1*n1)/(n2*n2);
iter = 1;
residual = zeros(maxiter,4);
res = Inf;
forRes = Inf;

while ((res > tol) && (iter <= maxiter))
    fprintf('Iter: %d\n', iter);
    % Fix G1,G2, update S11,S12.
    num = G1'*R11*G1;
    sG1 = G1'*G1;
    den = sG1*S11*sG1 + eps;
    S11 = S11.*(num./den);
    clear num den
    
    num = G2'*R22*G2;
    sG2 = G2'*G2;
    den = sG2*S22*sG2 + eps;
    S22 = S22.*(num./den);
    
    % Fix S11,S22, update G1,G2 alternatively.
    % Fix G2, update G1.
    num = 2*R11*(G1*S11) + lambda12*R12*G2+lamda1*W1*G1;
    den = 2*G1*(S11*sG1*S11) + lambda12*G1*sG2 +2*lamda2*D1*G1 + eps;
    
    G1 = G1.*(num./den);
    clear numr den sG1
    
    % Fix G1, update G2.
    num = 2*lambda22*R22*(G2*S22) + lambda12*R12'*G1++lamda2*W2*G2;
    sG1 = G1'*G1;
    den = 2*lambda22*G2*(S22*sG2*S22) + lambda12*G2*sG1+lamda2*D2*G2  + eps;
    G2 = G2.*(num./den);
    clear num den sG1 sG2
    
    residual(iter,1) = norm(R11 - G1*S11*G1','fro')^2;
    residual(iter,2) = lambda12*norm(R12 - G1*G2','fro')^2;
    residual(iter,3) = lambda22*norm(R22 - G2*S22*G2','fro')^2;
    residual(iter,4) = sum(residual(iter,1:3));
    
    res = abs(residual(iter,4) - forRes);
    forRes = residual(iter,4);
    iter = iter + 1;
    
    val = max(max(S11));
    sqrtval = sqrt(val);
    S11 = S11/val;
    G1 = G1*sqrtval;
    G2 = G2/sqrtval;
    S22 = S22*val;
    
    clear val sqrtval
    
    errorx1 = mean(mean(abs(R11 - G1*S11*G1')))/mean(mean(R11));
    errorx2 = mean(mean(abs(R22 - G2*S22*G2')))/mean(mean(R22));
    errorx3 = mean(mean(abs(R12 - G1*G2')))/mean(mean(R12));
    errorx(iter-1) = errorx1 + errorx2 + errorx3;
   
end

end