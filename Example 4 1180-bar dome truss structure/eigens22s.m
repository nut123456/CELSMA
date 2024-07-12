function [L, X ] = eigens22s(K, M, dof_constrnt)
    % ... (input/output comments unchanged) ...
    
    [ele_nod, ~] = size(K);
    dof = (1:ele_nod)';

    if nargin == 3
        % Handle constrained DOFs
        [K, M_reduced] = applyBoundaryConditions(K, M, dof_constrnt);
        free_dof = 1:size(K, 1);  % Get indices of unconstrained DOFs
    else
        % No boundary conditions applied
        K_reduced = K;
        M_reduced = M;
        free_dof = dof;
    end

    % Set options for eigs
    %h = waitbar(0, 'Calculating eigenvalues...'); % Create progress bar
    options.tol = 1e-8;   
    options.p = 20;       
    options.maxit = 500;  
    options.disp = 2;     
    %options.sigma = 0.1; % Example shift value - adjust based on your problem

    if nargout == 2
        % Compute eigenvalues and eigenvectors with boundary conditions
        num_eigenvalues = min(length(free_dof), 6);
        [X, D, flag] = eigs(K, M_reduced, num_eigenvalues, 'smallestabs', options);
        
        if flag ~= 0 
            warning('Not all eigenvalues converged. Consider adjusting parameters or using a different solver.');
        end
        
        d = diag(D);
        [L, ii] = sort(d);
        X = X(:, ii);
        X = expandEigenvectors(X, free_dof, ele_nod); 
        
        % Normalize eigenvectors (using the original, unreduced M)
        for jj = 1:num_eigenvalues
            mnorm = sqrt(X(:, jj)' * M * X(:, jj));
            X(:, jj) = X(:, jj) / mnorm;
            %waitbar(jj/num_eigenvalues, h); % Update progress bar
        end
    else
        % Compute only eigenvalues with boundary conditions
        d = eigs(K, M_reduced, 'smallestabs', options);
        L = sort(d);
    end
end
%close(h); % Close progress bar
function [K_reduced, M_reduced] = applyBoundaryConditions(K, M, dof_constrnt)
    % Remove rows and columns corresponding to constrained DOFs
    dof_constrnt = unique(dof_constrnt);  % Remove duplicates
    dof_constrnt = dof_constrnt(dof_constrnt > 0 & dof_constrnt <= size(K,1)); % Keep valid indices
    free_dof = setdiff(1:size(K, 1), dof_constrnt);
    K_reduced = K(free_dof, free_dof);
    M_reduced = M(free_dof, free_dof);
end

function X_expanded = expandEigenvectors(X, dof, ele_nod)
    % Expand eigenvectors back to original size
    X_expanded = zeros(ele_nod, size(X, 2));
    X_expanded(dof, :) = X;
end
