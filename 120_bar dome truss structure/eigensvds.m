function [L, X] = eigensvds(K, M, dof_constrnt)
% PURPOSE
%  Solve the generalized eigenvalue problem
%  [K - L*M]X = 0, considering boundary conditions.
%
% INPUT:
%    K : global stiffness matrix, dim(K) = ele_nod x ele_nod
%    M : global mass matrix, dim(M) = ele_nod x ele_nod
%    dof_constrnt : boundary condition matrix
%        dim(dof_active) = dof_active x 1
% OUTPUT:
%    L : eigenvalues stored in a vector with length (ele_nod - dof_active) 
%    X : eigenvectors dim(X) = ele_nod x num_dof, num_dof : number of dof's

[ele_nod, ~] = size(K);
dof = (1:ele_nod)';

if nargin == 3
    dof_constrnt = dof_constrnt(:);
    dof(dof_constrnt) = [];
    if nargout == 2
        % Use svds to solve the generalized eigenvalue problem
        [U, S, V] = svds(M(dof, dof), length(dof));
        M_inv_sqrt = U * diag(1 ./ sqrt(diag(S))) * V';
        K_transformed = M_inv_sqrt' * K(dof, dof) * M_inv_sqrt;
        [X1, D] = eig(K_transformed);
        X1 = M_inv_sqrt * X1;
        d = diag(D);
        [L, ii] = sort(d);
        X2 = X1(:, ii);
        X = zeros(ele_nod, length(dof));
        X(dof, :) = X2;
        % Normalize eigenvectors
        for jj = 1:length(dof)
            mnorm = sqrt(X(:, jj)' * M * X(:, jj));
            X(:, jj) = X(:, jj) / mnorm;
        end
    else
        [U, S, V] = svds(M(dof, dof), length(dof));
        M_inv_sqrt = U * diag(1 ./ sqrt(diag(S))) * V';
        K_transformed = M_inv_sqrt' * K(dof, dof) * M_inv_sqrt;
        d = eig(K_transformed);
        L = sort(d);
    end
else
    if nargout == 2
        % Use svds to solve the generalized eigenvalue problem
        [U, S, V] = svds(M, ele_nod);
        M_inv_sqrt = U * diag(1 ./ sqrt(diag(S))) * V';
        K_transformed = M_inv_sqrt' * K * M_inv_sqrt;
        [X1, D] = eig(K_transformed);
        X1 = M_inv_sqrt * X1;
        d = diag(D);
        [L, ii] = sort(d);
        X = X1(:, ii);
        % Normalize eigenvectors
        for jj = 1:ele_nod
            mnorm = sqrt(X(:, jj)' * M * X(:, jj));
            X(:, jj) = X(:, jj) / mnorm;
        end
    else
        [U, S, V] = svds(M, ele_nod);
        M_inv_sqrt = U * diag(1 ./ sqrt(diag(S))) * V';
        K_transformed = M_inv_sqrt' * K * M_inv_sqrt;
        d = eig(K_transformed);
        L = sort(d);
    end
end
end

