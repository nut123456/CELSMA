% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________  
function Positions = initialization(SearchAgents_no,dim,ub,lb)
% This function initialize the first population of search agents
Boundary_no = size(ub,2); % Numnber of boundaries
% If the boundaries of all variables are equal and user enter a signle
% Number for both ub and lb
if Boundary_no == 1
    Positions = rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no > 1
    for i = 1:dim
        ub_i = ub(i);
        lb_i = lb(i);
        Positions(:,i) = rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end