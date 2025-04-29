% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________   
function [rijitlik,elboy]=stiffness_3D_truss(topser,elsay,eldn,dnkoor,E,A)
rijitlik=zeros(topser);

for e=1:elsay
      
        indis=eldn(e,:);
        elsd=[indis(1)*3-2 indis(1)*3-1 indis(1)*3 indis(2)*3-2 indis(2)*3-1 indis(2)*3];
        xa=dnkoor(indis(2),1)-dnkoor(indis(1),1);
        ya=dnkoor(indis(2),2)-dnkoor(indis(1),2);
        za=dnkoor(indis(2),3)-dnkoor(indis(1),3);
        elboy(e)=sqrt(xa^2+ya^2+za^2);
        CX=xa/elboy(e);
        CY=ya/elboy(e);
        CZ=za/elboy(e);
        DD=sqrt(CX^2+CY^2);
        if CZ==1
        t=[ 0   0   1   0   0   0;
            0   1   0   0   0   0;
            -1  0   0   0   0   0;
            0   0   0   0   0   1;
            0   0   0   0   1   0;
            0   0   0   -1  0   0];
        elseif CZ==-1
        t=[ 0   0   -1  0   0   0;
            0   1   0   0   0   0;
            1   0   0   0   0   0;
            0   0   0   0   0   -1;
            0   0   0   0   1   0;
            0   0   0   1  0   0];
        else
        t=[ CX              CY          CZ  0           0           0;
            -(CY)/DD        (CX)/DD     0   0           0           0;
            -(CX*CZ)/DD     -(CY*CZ)/DD DD  0           0           0;
            0               0           0   CX          CY          CZ;
            0               0           0   -(CY)/DD    (CX)/DD     0;
            0               0           0   -(CX*CZ)/DD -(CY*CZ)/DD DD];
        end
        stiff_el=E*A(e)/elboy(e)*[1  0 0 -1 0 0;
                                  0  0 0  0 0 0;
                                  0  0 0  0 0 0;
                                  -1 0 0  1 0 0;
                                  0  0 0  0 0 0;
                                  0  0 0  0 0 0];
        stiff_glob=(t'*stiff_el)*t;
        rijitlik(elsd,elsd)=rijitlik(elsd,elsd)+stiff_glob;
end