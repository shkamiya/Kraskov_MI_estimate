function I = KSG_MI(X,k,varargin)
%% Kraskov's method (phase)
% +++ INPUT +++
% N: # of data points
% d: dimension of data
% X: = [x1, x2, ..., xd], column vectors of the same length
% k: the degree of k nearest neighbor method
% type: 1 for kraskov

N = size(X,1); d = size(X,2);
if nargin == 2 type = 1;
elseif nargin > 2 type = varargin{1}; end
type = 1;

switch type
    
%% Kraskov's method 1
    case 1
    n_vec1 = zeros(N,d);
    D_cheb = squareform(pdist(X,'chebychev'));
    for ii = 1:N
        min_k_vec = mink(D_cheb(:,ii),k+1); % k+1 because one needs to count the ii'th row itself
        epsilon = 2 * min_k_vec(k+1); 
        
        for jj = 1:d
            n_vec1(ii,jj) = sum(abs(X(:,jj) - X(ii,jj)) < epsilon/2) - 1;
        end
%        k = sum(D_cheb(:,n) < epsilon /2);
%         dist = abs(X - repmat(X(ii,:),N,1)); % coordinates' difference
%         [min_k_vec, Ind] = mink(max(dist,[],2), k+1); % k+1 because one needs to count the ii'th row itself
%         k_neighbor = dist(Ind(end),:); % coordinate of the k'th neighbor
%         n_vec1 = zeros(1,M);
%         for jj = 1:M
%             n_vec1(ii,jj) = sum(dist(:,jj) < k_neighbor(jj));
%         end
    end
    n_vec1 = n_vec1 + 1;
    psi_n_vector = psi(n_vec1);
    sum_psi_dim = sum(psi_n_vector,2);
    I = psi(k) + (d-1) * psi(N) - mean(sum_psi_dim);

%% Kraskov's method 2
%     case 2
%     n_vec2
%     for kk = 1:N
%         dist = abs(X - repmat(X(kk,:),N,1));
%         [min_k_mat, Ind] = mink(dist,k);
%         k_neighbor;
%         
%         %%%%% editing %%%%%
%         
%     end
% 
%% method 1 - another expression
case 3
    n_vec1 = zeros(N,d);
    [idx, D_cheb] = knnsearch(X,X,'Distance','chebychev','K',k+1); % k+1 because one needs to count the ii'th row itself
    for ii = 1:N
        min_k_vec = D_cheb(ii,k+1); 
        epsilon = 2 * min_k_vec(k+1); 
        
        for jj = 1:d
            n_vec1(ii,jj) = sum(abs(X(:,jj) - X(ii,jj)) < epsilon/2) - 1;
        end
    end
    n_vec1 = n_vec1 + 1;
    psi_n_vector = psi(n_vec1);
    sum_psi_dim = sum(psi_n_vector,2);
    I = psi(k) + (d-1) * psi(N) - mean(sum_psi_dim);

end
end

