function x = OMP(y, Phi, K)
% Othogonal matching pursuit of Tropp et Gilbert
% y : data
% Phi : sensing matrix
% K : sparsity
% SYM March 13 2013

[~,N] = size(Phi);
x = zeros(1,N);
S = []; % positions indexes of components of s
res = y; % first residual
PhiS = []; % Matrix of the columns used to represent y
for t=1:K;
    [~,j]=max(abs(Phi'*res));
    S = [S j];
    PhiS = [PhiS Phi(:,j)];
    x_est = pinv(PhiS)*y;
    res = y- PhiS*x_est;
    x(S) = x_est;
end;

