clearvars
close all
clc

disp('PPSO');

nvar = 47;
xmin = 0.0001;
xmax = 0.01;
dx = xmax-xmin;
vmax = 0.5*dx;
npop = 30;
maxit = 1;

for it=1:maxit
    if it==1
        gbestcost(1)=inf;
        for ii=1:npop
            velocity(ii,:)=zeros(1,nvar);
            theta(ii)=unifrnd(0,2*pi);
            position(ii,:)=xmin+(xmax-xmin)*rand(1,nvar);
            cost(ii)=Bar1410_truss(position(ii,:));
            
            pbest(ii,:)=position(ii,:);
            pbestcost(ii)=cost(ii);
            
            if pbestcost(ii)<gbestcost(it)
                gbest=pbest(ii,1:nvar);
                gbestcost(it)=pbestcost(ii);
            end
        end
    else
        
        gbestcost(it)=gbestcost(it-1);
        
        for ii=1:npop
            
            aa = 2*(sin(theta(ii)));
            bb = 2*(cos(theta(ii)));
            ee = abs(cos(theta(ii)))^aa;
            tt = abs(sin(theta(ii)))^bb;
            
            KK=(ee)*(pbest(ii,:)-position(ii,:)) +(tt)*(gbest-position(ii,:));
            
            velocity(ii,:)=((ee/(ii))*velocity(ii,:))+KK;
            
            velocity(ii,:)=min(max(velocity(ii,:),-vmax),vmax);
            
            position(ii,:)=position(ii,:)+velocity(ii,:);
            
            position(ii,:)=min(max(position(ii,:),xmin),xmax);
            
            cost(ii)=Bar1410_truss(position(ii,:));
            
            theta(ii)=theta(ii)+(abs(aa+bb)*(2*pi));
            
            vmax=(abs(cos(theta(ii)))^2)*dx;
            if cost(ii)<pbestcost(ii)
                pbest(ii,:)=position(ii,:);
                pbestcost(ii)=cost(ii);
                
                if pbestcost(ii)<gbestcost(it)
                    gbest=pbest(ii,:);
                    gbestcost(it)=pbestcost(ii);
                end
                
            end
        end
    end
  
    disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(gbestcost(it))]);
    
end
hold on;
figure(1)
it=1:10:maxit;

semilogy((it),(gbestcost(it)),'-ob','linewidth',1.4,'Markersize',3);

xlabel('\fontsize{12}\bf Iteration');
ylabel('\fontsize{12}\bf Best value');
legend('\fontsize{10}\bf PPSO');
Cost_Rsult=gbestcost(end);




