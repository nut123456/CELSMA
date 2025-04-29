% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________  
function [c,ceq] = Cstr120(x)
E = 210e9;     % Young's elastic modulus (kg/m^2)
rho = 7971.81;   % density of material (kg/m^3)
A = zeros(1,7);       % area of bar (m^2)
A1 = x(1);
A2 = x(2);
A3 = x(3);
A4 = x(4);
A5 = x(5);
A6 = x(6);
A7 = x(7);
A(1:12)= A1;
A(13:24)= A2;
A(25:36)= A3;
A(37:60)= A4;
A(61:84)= A5;
A(85:96)= A6;
A(97:120)= A7;

% number of elements
num_ele = 120;
% elements nodes
ele_nod = [1 2;1 3;1 4;1 5;1 6;1 7;1 8;1 9;1 10;1 11;1 12;1 13;2 3;3 4;4 5;5 6;6 7;7 8;8 9;9 10;10 11;11 12;12 13;13 2;2 14;3 16;4 18;5 20;6 22;7 24;8 26;9 28;10 30;11 32;12 34;13 36;...
           2 15;15 3;3 17;17 4;4 19;19 5;5 21;21 6;6 23;23 7;7 25;25 8;8 27;27 9;9 29;29 10;10 31;31 11;11 33;33 12;12 35;35 13;13 37;37 2;14 15;15 16;16 17;17 18;18 19;19 20;20 21;21 22;22 23;23 24;24 25;25 26;26 27;27 28;28 29;29 30;30 31;31 32;32 33;33 34;34 35;35 36;36 37;37 14;...
           14 38;16 39;18 40;20 41;22 42;24 43;26 44;28 45;30 46;32 47;34 48;36 49;38 15;15 39;39 17;17 40;40 19;19 41;41 21;21 42;42 23;23 43;43 25;25 44;44 27;27 45;45 29;29 46;46 31;31 47;47 33;33 48;48 35;35 49;49 37;37 38];
% number of nodes
num_nod = 49;
% nodes coordinates
nod_coor = [0 0 7;6.94 0 5.85;6.01 3.47 5.85;3.47 6.01 5.85;0 6.94 5.85;-3.47 6.01 5.85;-6.01 3.47 5.85;-6.94 0 5.85;-6.01 -3.47 5.85;-3.47 -6.01 5.85;...
            0 -6.94 5.85;3.47 -6.01 5.85;6.01 -3.47 5.85;12.04 0 3;11.63 3.116 3;10.427 6.02 3;8.514 8.514 3;6.02 10.417 3;3.116 11.63 3;0 12.04 3;...
            -3.116 11.63 3;-6.02 10.417 3;-8.514 8.514 3;-10.427 6.02 3;-11.63 3.116 3;-12.04 0 3;-11.63 -3.116 3;-10.427 -6.02 3;-8.514 -8.514 3;-6.02 -10.417 3;...
            -3.116 -11.63 3;0 -12.04 3;3.116 -11.63 3;6.02 -10.417 3;8.514 -8.514 3;10.427 -6.02 3;11.63 -3.116 3;15.89 0 0;13.761 7.945 0;7.945 13.761 0;
            0 15.89 0;-7.945 13.761 0;-13.761 7.945 0;-15.89 0 0;-13.761 -7.945 0;-7.945 -13.761 0;0 -15.89 0;7.945 -13.761 0;13.761 -7.945 0];
dofPerNode = 3;
num_dof = num_nod*dofPerNode;
dof = 1:num_dof;
% calculate stiffness matrix 
%[K,M] = Cal_K_and_M(nod_coor,ele_nod,A,rho,E );
% calculate K and M

dofPerNode = 3;  % number of degree of freedom of one node
num_dof    = num_nod*dofPerNode; % total dgree of freedom of system
K = zeros(num_dof);
M = zeros(num_dof);
ele_dof = zeros(num_ele,6);
L = zeros(num_ele,1);
C = zeros(num_ele,3);
for ii=1:num_ele
   index1 = ele_nod(ii,1); 
   index2 = ele_nod(ii,2);
   dx = nod_coor(index2,1)-nod_coor(index1,1);
   dy = nod_coor(index2,2)-nod_coor(index1,2);
   dz = nod_coor(index2,3)-nod_coor(index1,3);
   % compute length of each bar
   L(ii) = sqrt(dx^2+dy^2+dz^2);
   C(ii,1) = dx/L(ii);
   C(ii,2) = dy/L(ii);
   C(ii,3) = dz/L(ii);
   ele_dof(ii,:) = [3*index1-2 3*index1-1 3*index1...
                    3*index2-2 3*index2-1 3*index2];
end
% Construct global stiffness matrix & Mass matrix
for ii = 1:num_ele
    index = ele_dof(ii,:);
    K(index,index) = K(index,index) + A(ii)*E/L(ii)*rotate3(C(ii,1),C(ii,2),C(ii,3)); % rotate() calls rotation matrix;
    M(index,index) = M(index,index) + rho*L(ii)*A(ii)*[2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0;0 1 0 0 2 0;0 0 1 0 0 2]/6;% lumped mass matrix
end

% boundary and loading conditions
dof_constrnt = [112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147]; % fixed constraints
% add non-structural mass 
dof_active = setdiff(dof,dof_constrnt);
%dof_m = [1:111];%[1 2 3 4 5 6 7 8 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60 63 66 69 72 75 78 81 84 87 90 93 96 99 102 105 108 111]; % load at nodes
%add_M1 = 3000;
%add_M2 = 500;
%add_M3 = 100;
add_M = diag([3000;3000;3000;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;500;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]);
%add_M = diag([0;0;3000;0;0;500;0;0;500;0;0;500;0;0;500;0;0;500;0;0;500;0;0;500;0;0;500;0;0;500;0;0;500;0;0;500;0;0;500;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;100;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]);
%M(dof_m) = [3000 3000 3000 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100]; 
%M(dof_m) = [0 0 3000 0 0 500 0 0 500 0 0 500 0 0 500 0 0 500 0 0 500 0 0 500 0 0 500 0 0 500 0 0 500 0 0 500 0 0 500 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100 0 0 100]; 
%dof_1 = 3;
%dof_2 = ([6 9 12 15 18 21 24 27 30 33 36 39]);
%dof_3 = ([42 45 48 51 54 57 60 63 66 69 72 75 78 81 84 87 90 93 96 99 102 105 108 111]);
for dof = 1:num_dof
    %for ii = 1:size(add_M,1)
        %for jj = 1:size(add_M,2)
            %if ii == jj
              %M(dof,dof) = M(dof,dof) + add_M(ii,jj);
            %end
        %end
    %end
   M(dof,dof) = M(dof,dof) + add_M(dof,dof);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%% Version 1.0 
% %%Calculate frequency
% [omega, ~] = eigensss(K, M, dof_constrnt);
% omega = real(omega); % Filter out any imaginary parts
% 
% % Calculate frequencies
% f = sqrt(omega) / (2 * pi);
% 
% % Constraints
% c1 = 7 / f(1) - 1;
% c2 = 9 / f(3) - 1;
% c = [c1, c2];
% ceq = [];
% 
% % Display results
% %disp('Frequencies:');
% %disp(c)
% 
% 
% function [L,X]=eigens(K,M,dof_constrnt)
% % PURPOSE
% %  Solve the generalized eigenvalue problem
% %  [K-LM]X = 0, considering boundary conditions.
% %
% % INPUT:
% %    K : global stiffness matrix, dim(K)= ele_nod x ele_nod
% %    M : global mass matrix, dim(M)= ele_nod x ele_nod
% %    b : boundary condition matrix
% %        dim(dof_active)= dof_active x 1
% % OUTPUT:
% %    L : eigenvalues stored in a vector with length (ele_nod - dof_active) 
% %    X : eigenvectors dim(X)= ele_nod x num_dof, num_dof : number of dof's
% 
%   [ele_nod,ele_nod]=size(K);
%   dof=(1:ele_nod)';
% %
%   if nargin==3
%     dof_constrnt=dof_constrnt(:);
%     dof(dof_constrnt)=[]; 
%     if nargout==2
%       [X1,D]=eig(K(dof,dof),M(dof,dof));
%       [num_dof,num_dof]=size(X1);
%       for jj=1:num_dof
%         mnorm=sqrt(X1(:,jj)'*M(dof,dof)*X1(:,jj));
%         X1(:,jj)=X1(:,jj)/mnorm;
%       end
%       d=diag(D);
%       [L,ii]=sort(d);
%       X2=X1(:,ii);
%       X=zeros(ele_nod,num_dof);
%       X(dof,:)=X2;
%     else
%       d=eig(K(dof,dof),M(dof,dof));
%       L=sort(d);
%     end
%   else
%     if nargout==2
%       [X1,D]=eig(K,M);
%       for jj=1:ele_nod
%         mnorm=sqrt(X1(:,jj)'*M*X1(:,jj));
%         X1(:,jj)=X1(:,jj)/mnorm;
%       end
%       d=diag(D);
%       [L,ii]=sort(d);
%       X=X1(:,ii);
%     else
%       d=eig(K,M);
%       L=sort(d);
%     end
%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %version 2.0
% % Calculate frequency
% [omega, ~] = eigens(K, M, dof_constrnt);
% omega = real(omega); % Filter out any imaginary parts
% 
% % Calculate frequencies
% f = sqrt(omega) / (2 * pi);
% 
% % Constraints
% c1 = 7 / f(1) - 1;
% c2 = 9 / f(3) - 1;
% c = [c1, c2];
% ceq = [];
% 
% % % Display results
% % disp('Frequencies:');
% % disp(f)
% 
% function [L, X] = eigens(K, M, dof_constrnt)
% % PURPOSE
% %  Solve the generalized eigenvalue problem
% %  [K-LM]X = 0, considering boundary conditions.
% %
% % INPUT:
% %    K : global stiffness matrix, dim(K)= ele_nod x ele_nod
% %    M : global mass matrix, dim(M)= ele_nod x ele_nod
% %    dof_constrnt : boundary condition matrix
% %        dim(dof_constrnt)= dof_constrnt x 1
% % OUTPUT:
% %    L : eigenvalues stored in a vector with length (ele_nod - length(dof_constrnt)) 
% %    X : eigenvectors dim(X)= ele_nod x num_dof, num_dof : number of dof's
% 
%   [ele_nod, ~] = size(K);
%   dof = (1:ele_nod)';
% 
%   if nargin == 3
%     dof(dof_constrnt) = []; 
%     reduced_K = K(dof, dof);
%     reduced_M = M(dof, dof);
%   else
%     reduced_K = K;
%     reduced_M = M;
%   end
% 
%   if nargout == 2
%     [X1, D] = eig(reduced_K, reduced_M);
%     d = diag(D);
%     [L, ii] = sort(d);
%     X = zeros(ele_nod, length(d));
%     X(dof, :) = X1(:, ii);
%     % Normalize eigenvectors
%     for jj = 1:length(dof)
%       mnorm = sqrt(X(dof, jj)' * M(dof, dof) * X(dof, jj));
%       X(dof, jj) = X(dof, jj) / mnorm;
%     end
%   else
%     d = eig(reduced_K, reduced_M);
%     L = sort(d);
%   end
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Version 3.0 very fast
% Calculate frequency
[omega, ~] = eigens(K, M, dof_constrnt);
omega = real(omega); % Filter out any imaginary parts

% Calculate frequencies
f = sqrt(omega) / (2 * pi);

% Constraints
c1 = 9 / f(1) - 1;
c2 = 11 / f(3) - 1;
c = [c1, c2];
ceq = [];

% Display results
% disp('Frequencies:');
% disp(f)

function [L, X] = eigens(K, M, dof_constrnt)
% PURPOSE
%  Solve the generalized eigenvalue problem
%  [K-LM]X = 0, considering boundary conditions.
%
% INPUT:
%    K : global stiffness matrix, dim(K)= ele_nod x ele_nod
%    M : global mass matrix, dim(M)= ele_nod x ele_nod
%    b : boundary condition matrix
%        dim(dof_active)= dof_active x 1
% OUTPUT:
%    L : eigenvalues stored in a vector with length (ele_nod - dof_active) 
%    X : eigenvectors dim(X)= ele_nod x num_dof, num_dof : number of dof's

  [ele_nod, ~] = size(K);
  dof = (1:ele_nod)';

  if nargin == 3
    dof_constrnt = dof_constrnt(:);
    dof(dof_constrnt) = []; 
    if nargout == 2
      [X1, D] = eig(K(dof, dof), M(dof, dof));
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
      d = eig(K(dof, dof), M(dof, dof));
      L = sort(d);
    end
  else
    if nargout == 2
      [X1, D] = eig(K, M);
      d = diag(D);
      [L, ii] = sort(d);
      X = X1(:, ii);
      % Normalize eigenvectors
      for jj = 1:ele_nod
        mnorm = sqrt(X(:, jj)' * M * X(:, jj));
        X(:, jj) = X(:, jj) / mnorm;
      end
    else
      d = eig(K, M);
      L = sort(d);
    end
  end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%end

