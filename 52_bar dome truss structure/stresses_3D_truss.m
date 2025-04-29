% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________   
function [stresses]=stresses_3D_truss(elsay,eldn,dnkoor,E,yer,elboy)

%% Stresses
for e=1:elsay
      
        indis=eldn(e,:);
        elsd=[indis(1)*3-2 indis(1)*3-1 indis(1)*3 indis(2)*3-2 indis(2)*3-1 indis(2)*3];
        xa=dnkoor(indis(2),1)-dnkoor(indis(1),1);
        ya=dnkoor(indis(2),2)-dnkoor(indis(1),2);
        za=dnkoor(indis(2),3)-dnkoor(indis(1),3);
        CX=xa/elboy(e);
        CY=ya/elboy(e);
        CZ=za/elboy(e);
        memFORCE(e)=E/elboy(e)*[-CX -CY -CZ CX CY CZ]*yer(elsd);
end

stresses=memFORCE';
