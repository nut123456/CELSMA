function [L, X] = eigensGPU(K, M, dof_constrnt)
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
        K_cpu = gather(K(dof, dof));
        M_cpu = gather(M(dof, dof));
        num_eigenvalues = min(length(dof), 6); % Adjust this number as needed
        [X1, D] = eigs(K_cpu, M_cpu, num_eigenvalues, 'smallestabs');
        X1 = gpuArray(X1);
        d = diag(D);
        [L, ii] = sort(d);
        X2 = X1(:, ii);
        X = zeros(ele_nod, num_eigenvalues, 'gpuArray');
        X(dof, :) = X2;
        % Normalize eigenvectors
        for jj = 1:num_eigenvalues
            mnorm = sqrt(gather(X(:, jj)' * M * X(:, jj)));
            X(:, jj) = X(:, jj) / mnorm;
        end
        X = gather(X);
        L = gather(L);
    else
        K_cpu = gather(K(dof, dof));
        M_cpu = gather(M(dof, dof));
        num_eigenvalues = min(length(dof), 6); % Adjust this number as needed
        d = eigs(K_cpu, M_cpu, num_eigenvalues, 'smallestabs');
        L = sort(gather(d));
    end
else
    if nargout == 2
        % Use svds to solve the generalized eigenvalue problem
        K_cpu = gather(K);
        M_cpu = gather(M);
        num_eigenvalues = min(ele_nod, 6); % Adjust this number as needed
        [X1, D] = eigs(K_cpu, M_cpu, num_eigenvalues, 'smallestabs');
        X1 = gpuArray(X1);
        d = diag(D);
        [L, ii] = sort(d);
        X = X1(:, ii);
        % Normalize eigenvectors
        for jj = 1:ele_nod
            mnorm = sqrt(gather(X(:, jj)' * M * X(:, jj)));
            X(:, jj) = X(:, jj) / mnorm;
        end
        X = gather(X);
        L = gather(L);
    else
        K_cpu = gather(K);
        M_cpu = gather(M);
        num_eigenvalues = min(ele_nod, 6); % Adjust this number as needed
        d = eigs(K_cpu, M_cpu, num_eigenvalues, 'smallestabs');
        L = sort(gather(d));
    end
end
end
