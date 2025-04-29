% Please refer to the main paper:
% Slime Mould Algorithm: A New Method for Stochastic Optimization
% Shimin Li, Huiling Chen, Mingjing Wang, Ali Asghar Heidari, Seyedali Mirjalili
% Future Generation Computer Systems,2020
% DOI: https://doi.org/10.1016/j.future.2020.03.055
% https://www.sciencedirect.com/science/article/pii/S0167739X19320941
% ------------------------------------------------------------------------------------------------------------
% Website of SMA: http://www.alimirjalili.com/SMA.html
% You can find and run the SMA code online at http://www.alimirjalili.com/SMA.html

% You can find the SMA paper at https://doi.org/10.1016/j.future.2020.03.055
% Please follow the paper for related updates in researchgate: https://www.researchgate.net/publication/340431861_Slime_mould_algorithm_A_new_method_for_stochastic_optimization
%---------------------------------------------------------------------------------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all
clc
disp('The SMA is tracking the problem');
N=30; % Number of search agents

Function_name='F1'; % Name of the test function, range from F1-F13

T=500; % Maximum number of iterations



% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

Times=1; %Number of independent times you want to run the ESMA
display(['Number of independent runs: ', num2str(Times)]);

for i=1:Times
[Destination_fitness(i),bestPositions(i,:),Convergence_curve(i,:)]=SMA(N,T,lb,ub,dim,fobj);
display(['The optimal fitness of SMA is: ', num2str(Destination_fitness(i))]);
end

[bestfitness,index]=min(Destination_fitness);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of SMA is: ', num2str(bestfitness)]);
display(['The average fitness of SMA is: ', num2str(mean(Destination_fitness))]);
display(['The standard deviation fitness of SMA is: ', num2str(std(Destination_fitness))]);
display(['The best location of SMA is: ', num2str(bestPositions(index,:))]);

figure(3);
semilogy(Convergence_curve(index,:),'LineWidth',3);
xlabel('Iterations');
ylabel('Best fitness obtained so far');
legend('SMA');
box on;
axis tight;
grid off;

        



