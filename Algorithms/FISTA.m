function [x_hat,nIter] = SolveFISTA(A,b, varargin)
% FISTA
% Modified from Alen Yang's code at Berkeley
% b - m x 1 vector of observations/data (required input)
% A - m x n measurement matrix (required input)
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
% maxIter - maxilambdam number of iterations
%         - DEFAULT 10000, if omitted or -1.
% lineSearchFlag - 1 if line search is to be done every iteration
%                - DEFAULT 0, if omitted or -1.
% continuationFlag - 1 if a continuation is on lambda
%                  - DEFAULT 1, if omitted or -1.
% eta - line search parameter, should be in (0,1)
%     - ignored if lineSearchFlag is 0.
%     - DEFAULT 0.9, if omitted or -1.
% lambda - relaxation parameter
%    - ignored if continuationFlag is 1.
%    - DEFAULT 1e-3, if omitted or -1.
% outputFileName - Details of each iteration are dumped here
%
% x_hat - estimate of coeeficient vector
% numIter - number of iterations until convergence

DEBUG = 0 ;

STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_SUBGRADIENT;

stoppingCriterion = STOPPING_DEFAULT;
maxIter = 1000 ;
tolerance = 1e-3;
[m,n] = size(A) ;
x0 = zeros(n,1) ;
xG = [];

%% Initializing optimization variables
t_k = 1 ; 
t_km1 = 1 ;
L0 = 1 ;
G = A'*A ;
nIter = 0 ;
c = A'*b ;
lambda0 = 0.99*L0*norm(c,inf) ;
eta = 0.6 ;
lambda_bar = 1e-4*lambda0 ;
xk = zeros(n,1) ;
lambda = lambda0 ;
L = L0 ;
beta = 1.5 ;

% Parse the optional inputs.
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' lambdast be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue;
        case 'groundtruth'
            xG = parameterValue;
        case 'tolerance'
            tolerance = parameterValue;
        case 'linesearchflag'
            lineSearchFlag = parameterValue;
        case 'lambda'
            lambda_bar = parameterValue;
        case 'maxiteration'
            maxIter = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        case 'continuationflag'
            continuationFlag = parameterValue;
        case 'initialization'
            xk = parameterValue;
            if ~all(size(xk)==[n,1])
                error('The dimension of the initial xk does not match.');
            end
        case 'eta'
            eta = parameterValue;
            if ( eta <= 0 || eta >= 1 )
                disp('Line search parameter out of bounds, switching to default 0.9') ;
                eta = 0.9 ;
            end
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

if stoppingCriterion==STOPPING_GROUND_TRUTH && isempty(xG)
    error('The stopping criterion must provide the ground truth value of x.');
end

keep_going = 1 ;
nz_x = (abs(xk)> eps*10);
f = 0.5*norm(b-A*xk)^2 + lambda_bar * norm(xk,1);
xkm1 = xk;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;
    
    yk = xk + ((t_km1-1)/t_k)*(xk-xkm1) ;
    
    stop_backtrack = 0 ;
    
    temp = G*yk - c ; % gradient of f at yk
    
    while ~stop_backtrack
        
        gk = yk - (1/L)*temp ;
        
        xkp1 = soft(gk,lambda/L) ;
        
        temp1 = 0.5*norm(b-A*xkp1)^2 ;
        temp2 = 0.5*norm(b-A*yk)^2 + (xkp1-yk)'*temp + (L/2)*norm(xkp1-yk)^2 ;
        
        if temp1 <= temp2
            stop_backtrack = 1 ;
        else
            L = L*beta ;
        end
        
    end
    
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            keep_going = norm(xG-xkp1)>tolerance;
        case STOPPING_SUBGRADIENT
            sk = L*(yk-xkp1) + G*(xkp1-yk) ;
            keep_going = norm(sk) > tolerance*L*max(1,norm(xkp1));
        case STOPPING_SPARSE_SUPPORT
            % compute the stopping criterion based on the change
            % of the number of non-zero components of the estimate
            nz_x_prev = nz_x;
            nz_x = (abs(xkp1)>eps*10);
            num_nz_x = sum(nz_x(:));
            num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
            if num_nz_x >= 1
                criterionActiveSet = num_changes_active / num_nz_x;
                keep_going = (criterionActiveSet > tolerance);
            end
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            prev_f = f;
            f = 0.5*norm(b-A*xkp1)^2 + lambda_bar * norm(xk,1);
            criterionObjective = abs(f-prev_f)/(prev_f);
            keep_going =  (criterionObjective > tolerance);
        case STOPPING_DUALITY_GAP
            error('Duality gap is not a valid stopping criterion for PGBP.');
        otherwise
            error('Undefined stopping criterion.');
    end
    
    lambda = max(eta*lambda,lambda_bar) ;
    
    
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;
    
    t_km1 = t_k ;
    t_k = t_kp1 ;
    xkm1 = xk ;
    xk = xkp1 ;
end

x_hat = xk ;

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end