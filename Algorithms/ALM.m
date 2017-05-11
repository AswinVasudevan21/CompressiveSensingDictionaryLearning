% modified from Allen Yang's code
% SYM March 13 2013
function [x, nIter] = ALM(A, b)

% Initialize parameters
[m,n] = size(A) ;
tol = 1e-6 ;
tol_apg = 1e-6 ;
rho = 1.25 ;

At = A';
G = At*A ;
opts.disp = 0;
tau = eigs(G,1,'lm',opts);

nIter = 0 ;

mu = 1/tau ;
mubar = 1e10*mu ;
tauInv = 1*mu ;

lambda = ones(m,1) ;
x = zeros(n,1) ;

converged_main = 0 ;

maxIter = 200 ;
maxIter_apg = 50 ;


while ~converged_main
    nIter = nIter + 1 ;
    
    muInv = 1/mu ;
    lambdaScaled = muInv*lambda ;
    
    x_old_main = x ;
    
    
    temp = b + lambdaScaled ;
    temp = At*temp ;
    
    converged_apg = 0 ;
    
    nIter_apg = 0 ;
    
    t1 = 1 ; z = x ;
    
    muTauInv = muInv*tauInv ;
    
    while ~converged_apg
        
        nIter_apg = nIter_apg + 1 ;
        
        x_old_apg = x ;
        
        temp1 = z - tauInv*(G*z - temp) ;
        
        x = shrink(temp1, muTauInv) ;
        
        if norm(x_old_apg - x) < tol_apg*norm(x_old_apg)
            converged_apg = 1 ;
        end
        
        if nIter_apg >= maxIter_apg
            converged_apg = 1 ;
        end
        
        t2 = (1+sqrt(1+4*t1*t1))/2 ;
        z = x + ((t1-1)/t2)*(x-x_old_apg) ;
        t1 = t2 ;
        
    end
    
    lambda = lambda + mu*(b - A*x) ;
    mu = min(rho*mu,mubar) ;
    
    % Stopping criteria
    if norm(x_old_main - x) < tol*norm(x_old_main) % small change
        converged_main = 1 ;
    end
    
    if ~converged_main && norm(x_old_main-x) < 100*eps % small change
        converged_main = 1 ;
    end
    
    if ~converged_main && nIter >= maxIter % max Iter num
        converged_main = 1 ;
    end
    
end

function Y = shrink(X, alpha)

Y = sign(X).*max(abs(X)-alpha,0) ;
