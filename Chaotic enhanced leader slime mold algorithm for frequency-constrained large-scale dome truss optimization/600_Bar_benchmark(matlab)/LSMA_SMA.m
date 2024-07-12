clear all;
close all;
clc;

disp('The SMA is tracking the problem');
N=30; % Number of search agents
Function_name='F1'; % Name of the test function, range from F1-F13
T=500; % Maximum number of iterations

% Load details of the selected benchmark function
[lb, ub, dim, fobj] = Get_Functions_details(Function_name);

Times = 100; % Number of independent times you want to run the ESMA
display(['Number of independent runs: ', num2str(Times)]);

Convergence_curve_SMA_all = zeros(Times, T);
Convergence_curve_LSMA_all = zeros(Times, T);
Destination_fitness_SMA = zeros(1, Times);
bestPositions_SMA = zeros(Times, dim);
Destination_fitness_LSMA = zeros(1, Times);
bestPositions_LSMA = zeros(Times, dim);
time_SMA = zeros(1, Times);
time_LSMA = zeros(1, Times);

total_time_SMA = 0;
total_time_LSMA = 0;

for i = 1:Times
    tic;
    [Destination_fitness_SMA(i), bestPositions_SMA(i,:), Convergence_curve_SMA(i,:)] = SMA(N, T, lb, ub, dim, fobj);
    time_SMA(i) = toc;
    total_time_SMA = total_time_SMA + time_SMA(i);
    Convergence_curve_SMA_all(i,:) = Convergence_curve_SMA(i,:);
    display(['The optimal fitness of SMA is: ', num2str(Destination_fitness_SMA(i))]);
end

[bestfitness_SMA, index_SMA] = min(Destination_fitness_SMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of SMA is: ', num2str(bestfitness_SMA)]);
display(['The average fitness of SMA is: ', num2str(mean(Destination_fitness_SMA))]);
display(['The standard deviation fitness of SMA is: ', num2str(std(Destination_fitness_SMA))]);
display(['The best location of SMA is: ', num2str(bestPositions_SMA(index_SMA,:))]);
display(['Total time for SMA is: ', num2str(total_time_SMA)]);

disp('The LSMA is tracking the problem');
N=30; % Number of slime mould
Function_name='F1'; % Name of the test function that can be from F1 to F23 
MaxIT=500; % Maximum number of iterations

[lb, ub, dim, fobj] = Get_Functions_details(Function_name); % Function details

Times = 100; % Number of independent times you want to run the LSMA
display(['Number of independent runs: ', num2str(Times)]);

for i = 1:Times
    tic;
    [Destination_fitness_LSMA(i), bestPositions_LSMA(i,:), Convergence_curve_LSMA(i,:)] = LSMA(N, MaxIT, lb, ub, dim, fobj);
    time_LSMA(i) = toc;
    total_time_LSMA = total_time_LSMA + time_LSMA(i);
    Convergence_curve_LSMA_all(i,:) = Convergence_curve_LSMA(i,:);
    display(['The optimal fitness of LSMA is: ', num2str(Destination_fitness_LSMA(i))]);
end

[bestfitness_LSMA, index_LSMA] = min(Destination_fitness_LSMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of LSMA is: ', num2str(bestfitness_LSMA)]);
display(['The average fitness of LSMA is: ', num2str(mean(Destination_fitness_LSMA))]);
display(['The standard deviation fitness of LSMA is: ', num2str(std(Destination_fitness_LSMA))]);
display(['The best location of LSMA is: ', num2str(bestPositions_LSMA(index_LSMA,:))]);
display(['Total time for LSMA is: ', num2str(total_time_LSMA)]);

% Save results
save('optimization_results.mat', 'bestfitness_SMA', 'Destination_fitness_SMA', 'bestPositions_SMA', 'Convergence_curve_SMA_all', 'time_SMA', 'total_time_SMA', 'bestfitness_LSMA', 'Destination_fitness_LSMA', 'bestPositions_LSMA', 'Convergence_curve_LSMA_all', 'time_LSMA', 'total_time_LSMA');

% Plotting convergence curves
figure;
semilogy(Convergence_curve_SMA(index_SMA,:), 'LineWidth', 3);
hold on;
semilogy(mean(Convergence_curve_SMA_all), '--', 'LineWidth', 3);
hold on;
semilogy(Convergence_curve_LSMA(index_LSMA,:), 'LineWidth', 3);
xlabel('Iterations');
hold on;
semilogy(mean(Convergence_curve_LSMA_all), '--', 'LineWidth', 3);
xlabel('Iteration');
ylabel('Weight (kg)');
legend('SMA','Average SMA', 'LSMA', 'Average LSMA');
box on;
axis tight;
grid on;
% Adjusting Y-axis limit
current_ylim = ylim;
ylim([current_ylim(1)*0.9, current_ylim(2)]); % Set lower limit to 90% of current lower limit
title('Convergence Curves of SMA and LSMA');

