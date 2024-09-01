% Chaotic enhanced leader slime mold algorithm  (CELSMA) source Code Version 1.0
% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________           

clearvars
close all
clc
disp('The CLSMA is tracking the problem');

N = 35;         % Number of slime molds
Function_name = 'F1';  % Replace 'F1' with your desired test function
MaxIT = 500;     % Maximum number of iterations

% Get function details (assuming you have a function for this)
[lb, ub, dim, fobj] = Get_Functions_details(Function_name);

Times = 5;      % Number of independent runs
display(['Number of independent runs: ', num2str(Times)]);

% Initialize arrays to store results
Destination_fitness = zeros(1, Times);
bestPositions = zeros(Times, dim);
Convergence_curve = zeros(Times, MaxIT);
best_fitness_overall = inf; % Initialize best_fitness_overall here
best_position_overall = []; 
% Main loop
best_chaos_index = 1;  % เริ่มต้นด้วย chaos map index ที่ 1

% เลือกวิธีการเลือก chaos map
chaos_selection_method = 1; % 1 = สุ่มทุก iteration, 2 = ใช้ที่ดีที่สุดจาก run ก่อนหน้า

for i = 1:Times
    if chaos_selection_method == 1
        % สุ่มเลือก chaos map ใหม่ในทุก iteration
        chaos_index = randi([5, 5]);
    else
        % ใช้ chaos map ที่ดีที่สุดจาก run ก่อนหน้า
        chaos_index = best_chaos_index;
    end

    [Destination_fitness(i), bestPositions(i,:), Convergence_curve(i,:)] = CLSMA2(N, MaxIT, lb, ub, dim, fobj, chaos_index);
    display(['The optimal fitness of LSMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness(i)), ' (Chaos Index: ', num2str(chaos_index), ')']);

    if Destination_fitness(i) < best_fitness_overall
        best_fitness_overall = Destination_fitness(i);
        best_position_overall = bestPositions(i,:);
        best_chaos_index = chaos_index; % บันทึก chaos index ที่ดีที่สุด
    end
end

% Post-processing and analysis
[bestfitness, index] = min(Destination_fitness);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of CLSMA is: ', num2str(bestfitness)]);
display(['The average fitness of CLSMA is: ', num2str(mean(Destination_fitness))]);
display(['The standard deviation fitness of CLSMA is: ', num2str(std(Destination_fitness))]);
display(['The best location of CLSMA is: ', num2str(bestPositions(index,:))]);

% Plot convergence curve of the best run
semilogy(Convergence_curve(index,:), 'LineWidth', 3);
xlabel('Iterations');
ylabel('Best fitness obtained so far');
legend('LSMA with Chaos'); 
%title(['Convergence Curve (Function: ' Function_name ', Chaos Index: ' num2str(chaos_index) ')']); % Add chaos index to title
box on;
axis tight;
grid off;
