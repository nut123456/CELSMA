clearvars
close all
clc

disp('CLPPSO');

nvar = 25;
xmin = 0.0001;
xmax = 0.01;
dx = xmax-xmin;
vmax = 0.5*dx;
npop = 30;
maxit = 10;

ai = zeros(npop, nvar);
f_pbest = 1 : npop;
f_pbest = repmat(f_pbest', 1, nvar);
t = 0 : 1 / (npop - 1) : 1;
t = 5 .* t;
Pc = 0.05 + (0.5 - 0.05) .* (exp(t) - exp(t(1))) ./ (exp(t(npop)) - exp(t(1)));
m = 0 .* ones(npop, 1);
cost_history = zeros(maxit, 1); % To store cost in each iteration
time_history = zeros(maxit, 1); % To store time in each iteration

for it = 1:maxit
    tic; % Start timing for the current iteration
    if it == 1
       gbestcost(1) = inf;    
       for ii = 1:npop
           velocity(ii, :) = zeros(1, nvar);
           theta(ii) = unifrnd(0, 2*pi);
           position(ii, :) = xmin + (xmax - xmin) * rand(1, nvar);
           cost(ii) = Bar600_truss(position(ii, :)); % Replace with your cost function
           MM(ii) = 7;
           FF(ii) = 0;
           pbest(ii, :) = position(ii, :);
           pbestcost(ii) = cost(ii);
            
           if pbestcost(ii) < gbestcost(it)
              gbest = pbest(ii, 1:nvar);
              gbestcost(it) = pbestcost(ii);
           end
        end
        
        for kk = 1:npop
            ar = randperm(nvar);
            ai(kk, ar(1 : m(kk))) = 1;
            fi1 = ceil(npop * rand(1, nvar));
            fi2 = ceil(npop * rand(1, nvar));
            fi = (pbestcost(fi1) < pbestcost(fi2)) .* fi1 + (pbestcost(fi1) >= pbestcost(fi2)) .* fi2;
            bi = ceil(rand(1, nvar) - 1 + Pc(kk));
            
            if bi == zeros(1, nvar)
                rc = randperm(nvar);
                bi(rc(1)) = 1;
            end
            
            f_pbest(kk, :) = bi .* fi + (1 - bi) .* f_pbest(kk, :);
        end
    else
        gbestcost(it) = gbestcost(it - 1);
        
        for ii = 1:npop
            if FF(ii) >= 2
               FF(ii) = 0;
               ai(ii, :) = zeros(1, nvar);
               f_pbest(ii, :) = ii .* ones(1, nvar);
               ar = randperm(nvar);
               ai(ii, ar(1 : m(ii))) = 1;
               fi1 = ceil(npop * rand(1, nvar));
               fi2 = ceil(npop * rand(1, nvar));
                
               fi = (pbestcost(fi1) < pbestcost(fi2)) .* fi1 + (pbestcost(fi1) >= pbestcost(fi2)) .* fi2;
               bi = ceil(rand(1, nvar) - 1 + Pc(ii));
                
               if bi == zeros(1, nvar)
                  rc = randperm(nvar);
                  bi(rc(1)) = 1;
               end
                
               f_pbest(ii, :) = bi .* fi + (1 - bi) .* f_pbest(ii, :);
            end
            
            aa = 2 * (sin(theta(ii)));
            bb = 2 * (cos(theta(ii)));
            ee = abs(cos(theta(ii)))^aa;
            tt = abs(sin(theta(ii)))^bb;
            
            for jj = 1:nvar
                velocity(ii, jj) = ((ee / (npop)) * velocity(ii, jj)) + (ee) * (pbest(f_pbest(ii, jj), jj) - position(ii, jj));
            end
            
            velocity(ii, :) = min(max(velocity(ii, :), -vmax), vmax);
            
            position(ii, :) = position(ii, :) + velocity(ii, :);
            
            position(ii, :) = min(max(position(ii, :), xmin), xmax);
            
            cost(ii) = Bar600_truss(position(ii, :)); % Replace with your cost function
            
            theta(ii) = theta(ii) + (abs(aa + bb) * (pi));
            
            vmax = (abs(cos(theta(ii)))^2) * dx;
            
            if cost(ii) < pbestcost(ii)
               FF(ii) = 0;
            else
               FF(ii) = FF(ii) + 1;
            end
            if cost(ii) < pbestcost(ii)
               pbest(ii, :) = position(ii, :);
               pbestcost(ii) = cost(ii);
            else
               if pbestcost(ii) < gbestcost(it)
                  gbest = pbest(ii, :);
                  gbestcost(it) = pbestcost(ii);
               end
            end
        end
    end
    time_history(it) = toc; % End timing for the current iteration
    cost_history(it) = gbestcost(it); % Store the best cost
    disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(gbestcost(it)) '   Time = ' num2str(time_history(it)) ' seconds']);
end
disp(['Total Time = ' num2str(sum(time_history)) ' seconds']);
hold on;
figure(1)
it = 1:1:maxit;

semilogy((it), (gbestcost(it)), '-or', 'linewidth', 1.4, 'Markersize', 3);

xlabel('\fontsize{12}\bf Iteration');
ylabel('\fontsize{12}\bf Best value');
legend('\fontsize{10}\bf CLPPSO');

Cost_Rsult = gbestcost(end);

% Display or save the time and cost history
figure(2)
subplot(2,1,1);
plot(1:maxit, cost_history, '-o');
xlabel('Iteration');
ylabel('Best Cost');
title('Cost History');

subplot(2,1,2);
plot(1:maxit, time_history, '-o');
xlabel('Iteration');
ylabel('Time (seconds)');
title('Time History');
