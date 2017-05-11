function x = SP(y, Phi, K)
% Subspace Pursuit by Wei Dai and Olgica Milenkovic
% y : data
% Phi : sensing matrix
% K : sparsity
% SYM March 13 2013

[~,N]=size(Phi);
y_r = y;
iter = 1;

cv = abs(y_r'*Phi);
[~, S] = sort(cv,'descend');
S = sort(S(1:K));
Phi_x = Phi(:,S);
S_save(iter,:) = S;

x_p = pinv(Phi_x) * y;
y_r = y - Phi_x*x_p;
norm_save(iter) = norm(y_r);

while 1
   iter = iter+1;

   % find set of size 2K
   cv = abs(y_r'*Phi );
   [~, S] = sort(cv,'descend');
   S = sort(S(1:K) );
   Snew = union(S_save(iter-1,:), S); % 2S now
   Phi_x = Phi(:,Snew);

   % find the most significant K indices
   x_p = pinv(Phi_x) * y;
   [~, S] = sort(abs(x_p) , 'descend' );
   S = Snew(S(1:K));
   S = sort(S);
   Phi_x = Phi(:,S);
   S_save(iter,:)=S;

   % calculate the residue
   x_p = pinv(Phi_x) * y;
   y_r = y - Phi_x*x_p;
   norm_save(iter) = norm(y_r);

   if ( norm_save(iter) == 0 ) || (norm_save(iter)/norm_save(iter-1) >= 1)
       break;
   end
end

x = zeros(N,1);
x(S_save(iter,:) ) = reshape(x_p,K,1);
