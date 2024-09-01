%% Benchmark Test functions
function [lb,ub,dim,fobj] = Get_Functions_details(F)
switch F
    case 'F1'
        fobj = @F1;
        % lb=0.0001;
        % ub=0.01293;
        % dim=7;
        dim=7;
        lb = repmat(0.0001, 1, dim); % Lower bound
        ub = repmat(0.01293, 1, dim);   % Upper bound
            
end

% F1

function  W = F1(x)
rho = 7971.81;   % density of material (kg/m^3)
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
% nodes coordinates
nod_coor = [0 0 7;6.94 0 5.85;6.01 3.47 5.85;3.47 6.01 5.85;0 6.94 5.85;-3.47 6.01 5.85;-6.01 3.47 5.85;-6.94 0 5.85;-6.01 -3.47 5.85;-3.47 -6.01 5.85;...
            0 -6.94 5.85;3.47 -6.01 5.85;6.01 -3.47 5.85;12.04 0 3;11.63 3.116 3;10.427 6.02 3;8.514 8.514 3;6.02 10.417 3;3.116 11.63 3;0 12.04 3;...
            -3.116 11.63 3;-6.02 10.417 3;-8.514 8.514 3;-10.427 6.02 3;-11.63 3.116 3;-12.04 0 3;-11.63 -3.116 3;-10.427 -6.02 3;-8.514 -8.514 3;-6.02 -10.417 3;...
            -3.116 -11.63 3;0 -12.04 3;3.116 -11.63 3;6.02 -10.417 3;8.514 -8.514 3;10.427 -6.02 3;11.63 -3.116 3;15.89 0 0;13.761 7.945 0;7.945 13.761 0;
            0 15.89 0;-7.945 13.761 0;-13.761 7.945 0;-15.89 0 0;-13.761 -7.945 0;-7.945 -13.761 0;0 -15.89 0;7.945 -13.761 0;13.761 -7.945 0];
L = zeros(num_ele,1);
W = 0;
for ii = 1:num_ele
    index1 = ele_nod(ii,1); 
    index2 = ele_nod(ii,2);
    dx = nod_coor(index2,1)-nod_coor(index1,1);
    dy = nod_coor(index2,2)-nod_coor(index1,2);
    dz = nod_coor(index2,3)-nod_coor(index1,3);
    % compute length of each bar
    L(ii) = sqrt(dx^2+dy^2+dz^2);
    % compute weight of structure
    W = W + rho*L(ii)*A(ii);
end
% constraint
[ineq,eq] = Cstr120(x);
% Set inequality constraint g(x) <= 0
g    = ineq;
% Set equality constraint h(x) = 0
if ~isempty(eq) 
    h = eq; 
    else h = []; 
end
%compute the degree of constraint violations according to 'normalization'
if ~isempty(h)
    delta = 0; % the error of equal constraints
    f2 = sum(max(0,g),2) + sum(max(0,abs(h)- delta),2);
else
    f2 = sum(max(0,g),2);
end

% Obtain the fitness
maxit = 600;
eps1 = 1;
for it=1:maxit
eps2 = 1.5 + 0.5 *(it/maxit);
fpenalty = (1+eps1*f2).^eps2;
end
W = W.*fpenalty;
end

end