function [mass_mat]=mass_truss(topser,elsay,eldn,dnkoor,GS,A)
mass_mat=zeros(topser);

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
            0   0   0   1   0   0];
        else
        t=[ CX              CY          CZ  0           0           0;
            -(CY)/DD        (CX)/DD     0   0           0           0;
            -(CX*CZ)/DD     -(CY*CZ)/DD DD  0           0           0;
            0               0           0   CX          CY          CZ;
            0               0           0   -(CY)/DD    (CX)/DD     0;
            0               0           0   -(CX*CZ)/DD -(CY*CZ)/DD DD];
        end


%         T=[CX^2 CX*CY CX*CZ;CY*CX CY^2 CY*CZ;CZ*CX CZ*CY CZ^2];

%         rijitlik(elsd,elsd)=rijitlik(elsd,elsd)+1/2*GS*A(e)/elboy(e)*[T -T; -T T];

        mass=1/2*GS*A(e)*elboy(e)*[1 0 0 0 0 0;
                                   0 1 0 0 0 0;
                                   0 0 1 0 0 0;
                                   0 0 0 1 0 0;
                                   0 0 0 0 1 0;
                                   0 0 0 0 0 1];

        mass_glo=(t'*mass)*t;

mass_mat(elsd,elsd)=mass_mat(elsd,elsd)+mass_glo;
end