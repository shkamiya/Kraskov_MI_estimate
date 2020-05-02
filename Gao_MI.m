function I = Gao_MI(X,k)
%% Gao's discrete-continuous mixture; Extention of Kraskov's method (phase)
% +++ INPUT +++
% N: # of data points
% d: dimension of data: this must be two in this funciton
% X: = [x1, x2, ..., xd], column vectors of the same length
% k: the degree of k nearest neighbor method

N = size(X,1); d = size(X,2);
if d~= 2
    error('Dimension must be 2 in this function!')
end

%% Gao's method using KSG method
n_vec = zeros(N,d+1); % n_vec: i'th row is [n_{x,i} n_{y,i} k_i]
D_cheb = squareform(pdist(X,'chebychev'));
for ii = 1:N
    min_k_vec = mink(D_cheb(:,ii),k+1); % k+1 because one needs to count the ii'th row itself
    epsilon = 2 * min_k_vec(k+1); 

    n_vec(ii,1) = sum(abs(X(:,1) - X(ii,1)) < epsilon/2) - 1;
    n_vec(ii,2) = sum(abs(X(:,2) - X(ii,2)) < epsilon/2) - 1;
    
    if epsilon == 0
        n_vec(ii,3) = sum(D_cheb(:,ii)==0) - 1;
    else
        n_vec(ii,3) = k;
    end
    n_vec(n_vec==-1) = 0;
end

sigma_psi_minus_log = sum([-log(n_vec(:,1:2)+1) psi(n_vec(:,3))], 2);
I = mean(sigma_psi_minus_log) + log(N);

% case 3
%     n_vec1 = zeros(N,d);
%     [idx, D_cheb] = knnsearch(X,X,'Distance','chebychev','K',k+1); % k+1 because one needs to count the ii'th row itself
%     for ii = 1:N
%         min_k_vec = D_cheb(ii,k+1); 
%         epsilon = 2 * min_k_vec(k+1); 
%         
%         for jj = 1:d
%             n_vec1(ii,jj) = sum(abs(X(:,jj) - X(ii,jj)) < epsilon/2) - 1;
%         end
%     end
%     n_vec1 = n_vec1 + 1;
%     psi_n_vector = psi(n_vec1);
%     sum_psi_dim = sum(psi_n_vector,2);
%     I = psi(k) + (d-1) * psi(N) - mean(sum_psi_dim);
% 
% end
end

