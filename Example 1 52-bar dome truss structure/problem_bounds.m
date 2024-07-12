function [LB UB]=problem_bounds(problem_type)
% LB : Lower Bounds
% UB : Upper Bounds


switch problem_type
            case 1
                LB=ones(1,17)*0.1;
                UB=ones(1,17)*50;
            case 2
                LB=ones(1,12);
                LB(1,1:4)=LB(1,1:4)*3.5;   %in2
                LB(1,5)=LB(1,5)*775;  %in (x3)
                LB(1,6)=LB(1,6)*-225;  %in (y3)
                LB(1,7)=LB(1,7)*525;  %in (x5)
                LB(1,8)=LB(1,8)*-225;  %in (y5)
                LB(1,9)=LB(1,9)*275;  %in (x7)
                LB(1,10)=LB(1,10)*-225;  %in (y7)
                LB(1,11)=LB(1,11)*25;  %in (x9)
                LB(1,12)=LB(1,12)*-225;  %in (y9)

                UB=ones(1,12);
                UB(1,1:4)=UB(1,1:4)*18;   %in2  
                UB(1,5)=UB(1,5)*1225;  %in (x3)
                UB(1,6)=UB(1,6)*245;  %in (y3)
                UB(1,7)=UB(1,7)*975;  %in (x5)
                UB(1,8)=UB(1,8)*245;  %in (y5)
                UB(1,9)=UB(1,9)*725;  %in (x7)
                UB(1,10)=UB(1,10)*245;  %in (y7)
                UB(1,11)=UB(1,11)*475;  %in (x9)
                UB(1,12)=UB(1,12)*245;  %in (y9)
            case 3
                LB=ones(1,29)*0.1;
                UB=ones(1,29)*20;
			case 4
                LB=ones(1,8)*0.01;
                UB=ones(1,8)*3.4;
			case 5
                LB=ones(1,16)*0.1;
                UB=ones(1,16)*20;
			case 6
                LB=ones(1,7)*0.775;
                UB=ones(1,7)*20;
            case 7
                LB=ones(1,10)*0.645e-4;
                UB=ones(1,10)*50e-4;
            case 8          			
                LB=ones(1,19);
                LB(1,1:5)=LB(1,1:5)*0.3;   %m
                LB(1,6:19)=LB(1,6:19)*1e-4;  %m2

                UB=ones(1,19);
                UB(1,1:5)=UB(1,1:5)*3;   %m
                UB(1,6:19)=UB(1,6:19)*10e-4;  %m2 
            case 9
                LB=ones(1,13);
                LB(1,1)=LB(1,1)*4;   %m (ZA)
                LB(1,2)=LB(1,2)*0.01;   %m  (XB)
                LB(1,3)=LB(1,3)*3.7;   %m (ZB)
                LB(1,4)=LB(1,4)*2.01;   %m (XF)
                LB(1,5)=LB(1,5)*2.5;   %m (ZF)
                LB(1,6:13)=LB(1,6:13)*0.0001;  %m2

                UB=ones(1,13);
                UB(1,1)=UB(1,1)*8;   %m (ZA)
                UB(1,2)=UB(1,2)*4;   %m  (XB)
                UB(1,3)=UB(1,3)*7.7;   %m (ZB)
                UB(1,4)=UB(1,4)*6;   %m (XF)
                UB(1,5)=UB(1,5)*6.5;   %m (ZF)
                UB(1,6:13)=UB(1,6:13)*0.001;  %m2
end