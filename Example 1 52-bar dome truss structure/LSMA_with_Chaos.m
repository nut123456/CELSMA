% Leader Slime Mould Algorithm (LSMA) source Code Version 1.0
%
% Developed in MATLAB R2018b
%
% Author and programmer: 
% Dr Manoj Kumar Naik
% Faculty of Engineering and Technology, Siksha O Anusandhan, Bhubaneswar, Odisha – 751030, India 
% e-mail:       naik.manoj.kumar@gmail.com
% ORCID:        https://orcid.org/0000-0002-8077-1811
% SCOPUS:       https://www.scopus.com/authid/detail.uri?authorId=35753522900
% Publons:      https://publons.com/researcher/2057920/manoj-kumar-naik/
% G-Scholar:    https://scholar.google.co.in/citations?user=tX-8Xw0AAAAJ&hl=en 
% Researchgate: https://www.researchgate.net/profile/Manoj_Naik9
% DBLP:         https://dblp.uni-trier.de/pers/k/Kumar:Naik_Manoj
%_____________________________________________________________________________________________________           
% Please cite to the main paper:
% ******************************
% M. K. Naik, R. Panda, and A. Abraham, “Normalized square difference based 
% multilevel thresholding technique for multispectral images using leader slime 
% mould algorithm,” J. King Saud Univ. - Comput. Inf. Sci., Nov. 2020, 
% doi: 10.1016/j.jksuci.2020.10.030.
%
% This program using the framework of SMA by Ali Asghar Heidari
% https://aliasgharheidari.com/SMA.html
%_____________________________________________________________________________________________________

function [Destination_fitness,bestPositions,Convergence_curve]=LSMA_with_Chaos(N,Max_iter,lb,ub,dim,fobj,chaos_index)
bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);  %Eq. (16)
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
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(20)
    worstFitness = SmellOrder(N); %Eq. (23)
    bestFitness = SmellOrder(1);  %Eq. (22)
    bestPositions2=X(SmellIndex(2),:);%Leader 2
    bestPositions3=X(SmellIndex(3),:);%Leader 3
    
    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero

    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(21)
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
     
    a = atanh(-(it/Max_iter)+1);   %Eq.(18)
    b = 1-it/Max_iter;             %Eq.(19) 
    % Update the Position of search agents
    for i=1:N
        if rand<z     %Eq.(24.a)
            chaos_value = chaos(chaos_index, it, Max_iter, 1); % ปรับค่า Value ตามต้องการ
            X(i,:) = (ub-lb)*chaos_value+lb;
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(17)
            vb = unifrnd(-a,a,1,dim);  
            vc = unifrnd(-b,b,1,dim);
            A = randi([1,N]);  % two positions randomly selected from population
            B = randi([1,N]);
            r1 = rand();
            for j=1:dim
                if r1<p    %Eq.(24.b)
                    X(i,j) = bestPositions(j)+ vb(j)*(((weight(i,j)*bestPositions2(j)...
                        -X(A,j))+(weight(i,j)*bestPositions3(j)-X(B,j))));
                else      %Eq.(24.c)
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    Convergence_curve(it)=Destination_fitness;
    time_per_iteration(it) = toc; % End timer and record elapsed time
    %disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(Convergence_curve(it)) ': Time =' num2str(time_per_iteration(it)) ]);
    it=it+1;
end
disp(['Time(s) = ' num2str(sum(time_per_iteration))]);

% Plot convergence curve (Best Cost over Iterations)
figure(1);
plot(1:Max_iter, Convergence_curve, '-ob', 'LineWidth', 1.4, 'MarkerSize', 3);
xlabel('\fontsize{12}\bf Iteration');
ylabel('\fontsize{12}\bf Best Cost');
title('\fontsize{12}\bf Convergence Curve(CLSMA)');
legend('\fontsize{10}\bf Best Cost');

% Plot time per iteration
figure(2);
plot(1:Max_iter, time_per_iteration, '-or', 'LineWidth', 1.4, 'MarkerSize', 3);
xlabel('\fontsize{12}\bf Iteration');
ylabel('\fontsize{12}\bf Time (s)');
title('\fontsize{12}\bf Computation Time per Iteration');
legend('\fontsize{10}\bf Time per Iteration');
