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

clearvars
close all
clc
disp('The LSMA is tracking the problem');

N = 20;         % Number of slime molds
Function_name = 'F1';  % Replace 'F1' with your desired test function
MaxIT = 500;     % Maximum number of iterations

% Get function details (assuming you have a function for this)
[lb, ub, dim, fobj] = Get_Functions_details(Function_name);

Times = 10;      % Number of independent runs
display(['Number of independent runs: ', num2str(Times)]);

% Initialize arrays to store results
Destination_fitness = zeros(1, Times);
bestPositions = zeros(Times, dim);
Convergence_curve = zeros(Times, MaxIT);

% Main loop
% สร้างชุดของ Chaos Index ที่มีความถี่เท่ากัน
%โค้ดนี้จะทำให้ chaos map แต่ละตัวถูกเลือกอย่างเท่าเทียมกันที่สุดเท่าที่จะเป็นไปได้ หาก Times หารด้วย 10 ไม่ลงตัว อาจมี chaos map บางตัวที่ถูกเลือกมากกว่าตัวอื่น ๆ อยู่ 1 ครั้ง
runs_per_map = ceil(Times / 10);
chaos_indices = repelem(1:10, runs_per_map);
chaos_indices = chaos_indices(randperm(length(chaos_indices)));
chaos_indices = chaos_indices(1:Times);
for i = 1:Times
    % Choose a chaos map randomly for each run (optional)
    chaos_index = chaos_indices(i); 
    %chaos_index = randi([1, 10]);  
    [Destination_fitness(i), bestPositions(i,:), Convergence_curve(i,:)] = LSMA_with_Chaos(N, MaxIT, lb, ub, dim, fobj, chaos_index);
    display(['The optimal fitness of CLSMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness(i))]);
end

% Post-processing and analysis
[bestfitness, index] = min(Destination_fitness);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of LSMA is: ', num2str(bestfitness)]);
display(['The average fitness of LSMA is: ', num2str(mean(Destination_fitness))]);
display(['The standard deviation fitness of LSMA is: ', num2str(std(Destination_fitness))]);
display(['The best location of LSMA is: ', num2str(bestPositions(index,:))]);

% Plot convergence curve of the best run
semilogy(Convergence_curve(index,:), 'LineWidth', 3);
xlabel('Iterations');
ylabel('Best fitness obtained so far');
legend('LSMA with Chaos'); 
title(['Convergence Curve (Function: ' Function_name ', Chaos Index: ' num2str(chaos_index) ')']); % Add chaos index to title
box on;
axis tight;
grid off;