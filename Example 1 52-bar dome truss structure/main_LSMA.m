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

disp('The LSMA is tracking the problem');

N=20; % Number of slime mould
Function_name='F1' % Name of the test function that can be from F1 to F23 
MaxIT=500; % Maximum number of iterations

[lb,ub,dim,fobj]=Get_Functions_details(Function_name); % Function details

Times=1; %Number of independent times you want to run the AOSMA
display(['Number of independent runs: ', num2str(Times)]);

for i=1:Times
[Destination_fitness(i),bestPositions(i,:),Convergence_curve(i,:)]=LSMA(N,MaxIT,lb,ub,dim,fobj);
display(['The optimal fitness of LSMA is: ', num2str(Destination_fitness(i))]);
end

[bestfitness,index]=min(Destination_fitness);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of LSMA is: ', num2str(bestfitness)]);
display(['The average fitness of LSMA is: ', num2str(mean(Destination_fitness))]);
display(['The standard deviation fitness of LSMA is: ', num2str(std(Destination_fitness))]);
display(['The best location of LSMA is: ', num2str(bestPositions(index,:))]);

semilogy(Convergence_curve(index,:),'LineWidth',3);
xlabel('Iterations');
ylabel('Best fitness obtained so far');
legend('LSMA');
box on;
axis tight;
grid off;









