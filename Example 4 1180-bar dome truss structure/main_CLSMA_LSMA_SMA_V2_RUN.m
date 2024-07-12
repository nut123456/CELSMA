clear all 
close all
clc

disp('The SMA is tracking the problem');
N=20; % Number of search agents
Function_name='F1'; % Name of the test function, range from F1-F13
MaxIT=500; % Maximum number of iterations

% Load details of the selected benchmark function
[lb, ub, dim, fobj] = Get_Functions_details(Function_name);

Times = 40; % Number of independent times you want to run the SMA, LSMA, CLSMA
display(['Number of independent runs: ', num2str(Times)]);

Convergence_curve_SMA_all = zeros(Times, MaxIT);
Convergence_curve_LSMA_all = zeros(Times, MaxIT);
Convergence_curve_CLSMA_all = zeros(Times, MaxIT);
best_chaos_index = 1; 
chaos_selection_method = 1; % 1 = random per iteration, 2 = best from previous run
for i = 1:Times
    [Destination_fitness_SMA(i), bestPositions_SMA(i,:), Convergence_curve_SMA(i,:)] = SMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_SMA_all(i,:) = Convergence_curve_SMA(i,:);
    display(['The optimal fitness of SMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_SMA(i))]);
end



disp('The LSMA is tracking the problem');
%N=30; % Number of slime mould
Function_name='F1'; % Name of the test function that can be from F1 to F23 
%MaxIT=250; % Maximum number of iterations

[lb, ub, dim, fobj] = Get_Functions_details(Function_name); % Function details

%Times = 1; % Number of independent times you want to run the LSMA
display(['Number of independent runs: ', num2str(Times)]);

for i = 1:Times
    [Destination_fitness_LSMA(i), bestPositions_LSMA(i,:), Convergence_curve_LSMA(i,:)] = LSMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_LSMA_all(i,:) = Convergence_curve_LSMA(i,:);
    display(['The optimal fitness of LSMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_LSMA(i))]);
end



disp('The CLSMA is tracking the problem');
%N=30; % Number of slime mould
Function_name='F1'; % Name of the test function that can be from F1 to F23 
%MaxIT=500; % Maximum number of iterations

[lb, ub, dim, fobj] = Get_Functions_details(Function_name); % Function details

%Times = 1; % Number of independent times you want to run the LSMA
display(['Number of independent runs: ', num2str(Times)]);

for i = 1:Times
    if chaos_selection_method == 1
        chaos_index = randi([1, 10]);  
    else
        chaos_index = best_chaos_index;
    end
    [Destination_fitness_CLSMA(i), bestPositions_CLSMA(i,:), Convergence_curve_CLSMA(i,:)] = CLSMA2(N, MaxIT, lb, ub, dim, fobj, chaos_index);
    Convergence_curve_CLSMA_all(i,:) = Convergence_curve_CLSMA(i,:);
    display(['The optimal fitness of CLSMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_CLSMA(i))]);
end

[bestfitness_SMA, index_SMA] = min(Destination_fitness_SMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of SMA is: ', num2str(bestfitness_SMA)]);
display(['The average fitness of SMA is: ', num2str(mean(Destination_fitness_SMA))]);
display(['The standard deviation fitness of SMA is: ', num2str(std(Destination_fitness_SMA))]);
display(['The best location of SMA is: ', num2str(bestPositions_SMA(index_SMA,:))]);

[bestfitness_LSMA, index_LSMA] = min(Destination_fitness_LSMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of LSMA is: ', num2str(bestfitness_LSMA)]);
display(['The average fitness of LSMA is: ', num2str(mean(Destination_fitness_LSMA))]);
display(['The standard deviation fitness of LSMA is: ', num2str(std(Destination_fitness_LSMA))]);
display(['The best location of LSMA is: ', num2str(bestPositions_LSMA(index_LSMA,:))]);

[bestfitness_CLSMA, index_CLSMA] = min(Destination_fitness_CLSMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of CLSMA is: ', num2str(bestfitness_CLSMA)]);
display(['The average fitness of CLSMA is: ', num2str(mean(Destination_fitness_CLSMA))]);
display(['The standard deviation fitness of CLSMA is: ', num2str(std(Destination_fitness_CLSMA))]);
display(['The best location of CLSMA is: ', num2str(bestPositions_LSMA(index_CLSMA,:))]);



% Plotting  convergence curves
figure;
semilogy(Convergence_curve_SMA(index_SMA,:), 'LineWidth', 2.5);
hold on;
semilogy(mean(Convergence_curve_SMA_all), '--', 'LineWidth', 2.5);
hold on;

semilogy(Convergence_curve_LSMA(index_LSMA,:), 'LineWidth', 2.5);
hold on;
semilogy(mean(Convergence_curve_LSMA_all), '--', 'LineWidth', 2.5);
hold on;

semilogy(Convergence_curve_CLSMA(index_CLSMA,:), 'LineWidth', 2.5);
hold on;
semilogy(mean(Convergence_curve_CLSMA_all), '--', 'LineWidth', 2.5);
hold on;

xlabel('Iteration');
ylabel('Weight (kg)');


legend('Best, SMA', 'Mean, SMA', 'Best, LSMA', 'Mean, LSMA', 'Best, CLSMA', 'Mean, CLSMA');
title(['Example 5']);
box on; axis tight; grid on;
current_ylim = ylim;
ylim([current_ylim(1)*0.85, current_ylim(2)]);
% Save the figure
saveas(gcf, 'convergence_curves_Example5.png');
