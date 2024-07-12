% Chaotic enhanced leader slime mold algorithm  (CELSMA) source Code Version 1.0
% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________           


function [Destination_fitness,bestPositions,Convergence_curve]=LSMA(N,Max_iter,lb,ub,dim,fobj)
bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);  
Convergence_curve=zeros(1,Max_iter);

it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

time_per_iteration = zeros(1, Max_iter); % Preallocate array for time per iteration
% Main loop
while  it <= Max_iter
    tic; % Start timer for iteration  
    %sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
    end   
    [SmellOrder,SmellIndex] = sort(AllFitness); 
    worstFitness = SmellOrder(N); 
    bestFitness = SmellOrder(1);  
    bestPositions2=X(SmellIndex(2),:);%Leader 2
    bestPositions3=X(SmellIndex(3),:);%Leader 3
    
    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero

    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(SmellIndex(1),:); %Leader 1
        Destination_fitness = bestFitness;
    end
     
    a = atanh(-(it/Max_iter)+1); 
    b = 1-it/Max_iter;            
    % Update the Position of search agents
    for i=1:N
        if rand<z     %Eq.(24.a)
            X(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness)); 
            vb = unifrnd(-a,a,1,dim);  
            vc = unifrnd(-b,b,1,dim);
            A = randi([1,N]);  % two positions randomly selected from population
            B = randi([1,N]);
            r1 = rand();
            for j=1:dim
                if r1<p    
                    X(i,j) = bestPositions(j)+ vb(j)*(((weight(i,j)*bestPositions2(j)...
                        -X(A,j))+(weight(i,j)*bestPositions3(j)-X(B,j))));
                else      
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    Convergence_curve(it)=Destination_fitness;
    time_per_iteration(it) = toc; % End timer and record elapsed time
   % disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(Convergence_curve(it)) ': Time =' num2str(time_per_iteration(it)) ]);
    it=it+1;
end
disp(['Time(s) = ' num2str(sum(time_per_iteration))]);

% % Plot convergence curve (Best Cost over Iterations)
% figure(1);
% plot(1:Max_iter, Convergence_curve, '-ob', 'LineWidth', 1.4, 'MarkerSize', 3);
% xlabel('\fontsize{12}\bf Iteration');
% ylabel('\fontsize{12}\bf Best Cost');
% title('\fontsize{12}\bf Convergence Curve');
% legend('\fontsize{10}\bf Best Cost');
% 
% % Plot time per iteration
% figure(2);
% plot(1:Max_iter, time_per_iteration, '-or', 'LineWidth', 1.4, 'MarkerSize', 3);
% xlabel('\fontsize{12}\bf Iteration');
% ylabel('\fontsize{12}\bf Time (s)');
% title('\fontsize{12}\bf Computation Time per Iteration');
% legend('\fontsize{10}\bf Time per Iteration');

