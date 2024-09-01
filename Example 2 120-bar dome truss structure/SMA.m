function [Destination_fitness,bestPositions,Convergence_curve]=SMA(N,Max_iter,lb,ub,dim,fobj)
disp('SMA is now tackling your problem')

% Initialize position
bestPositions=zeros(1,dim);
Destination_fitness=inf; % Change this to -inf for maximization problems
AllFitness = inf*ones(N,1); % Record the fitness of all slime mold
weight = ones(N,dim); % Fitness weight of each slime mold
% Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
it=1;  % Number of iterations
lb=ones(1,dim).*lb; % Lower boundary 
ub=ones(1,dim).*ub; % Upper boundary
z=0.03; % Parameter
time_per_iteration = zeros(1, Max_iter); % Preallocate array for time per iteration

% Main loop
while it <= Max_iter
    tic; % Start timer for iteration  
    % Sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
    end
    
    [SmellOrder,SmellIndex] = sort(AllFitness);  % Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);

    S=bestFitness-worstFitness+eps;  % Plus eps to avoid denominator zero

    % Calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i <= (N/2)  % Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    % Update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
    a = atanh(-(it/Max_iter)+1);   % Eq.(2.4)
    b = 1-it/Max_iter;
    % Update the Position of search agents
    for i=1:N
        if rand < z     % Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        else
            p = tanh(abs(AllFitness(i)-Destination_fitness));  % Eq.(2.2)
            vb = unifrnd(-a,a,1,dim);  % Eq.(2.3)
            vc = unifrnd(-b,b,1,dim);
            for j=1:dim
                r = rand();
                A = randi([1,N]);  % Two positions randomly selected from population
                B = randi([1,N]);
                if r < p    % Eq.(2.1)
                    X(i,j) = bestPositions(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    Convergence_curve(it)=Destination_fitness;
    time_per_iteration(it) = toc; % End timer and record elapsed time
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(Convergence_curve(it)) ': Time =' num2str(time_per_iteration(it)) ' seconds']);
    it=it+1;
end

disp(['Total Time(s) = ' num2str(sum(time_per_iteration))]);

% Save time per iteration to a file
fileID = fopen('time_per_iteration.txt','w');
fprintf(fileID, 'Time per iteration (seconds):\n');
fprintf(fileID, '%f\n', time_per_iteration);
fclose(fileID);

% % Plot convergence curve (Best Cost over Iterations)
% figure;
% subplot(2,1,1);
% plot(1:Max_iter, Convergence_curve, '-ob', 'LineWidth', 1.4, 'MarkerSize', 3);
% xlabel('\fontsize{12}\bf Iteration');
% ylabel('\fontsize{12}\bf Best Cost');
% title('\fontsize{12}\bf Convergence Curve');
% legend('\fontsize{10}\bf Best Cost');
% 
% % Plot time per iteration
% subplot(2,1,2);
% plot(1:Max_iter, time_per_iteration, '-or', 'LineWidth', 1.4, 'MarkerSize', 3);
% xlabel('\fontsize{12}\bf Iteration');
% ylabel('\fontsize{12}\bf Time (s)');
% title('\fontsize{12}\bf Computation Time per Iteration');
% legend('\fontsize{10}\bf Time per Iteration');
