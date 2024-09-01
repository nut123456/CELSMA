function [TRUSS]=Truss_modal_analysis_52bar(Individual)
%SMA% Individual=[5.7036,	2.1294,	3.7268,	3.9222,	2.5001,	0.0001,
%0.00013971,	0.00012738,	0.00015981,	0.00014042,	0.0001,	0.00016354,	0.00015141]

%LSMA% Individual=[5.8221,	2.1102,3.7604	,3.9616	,2.5,	0.00010211,0.0001573,	0.00014216,	0.00013456,	0.00014141,	0.00010029,	0.00015474,	0.00015685]

%CELSMA Individual=[4.0063,	2.213,	3.7001,	3.9448,	2.5,	0.00010016,0.00010009,	0.00014864,	0.00019687,	0.00012841,	0.00010011,	0.00015517,	0.00018154]


Truss_Data=TrussData_52bar_modal;

% ZA=Individual(5)+Individual(3)+Individual(1);
% XB=Individual(2);
% ZB=Individual(5)+Individual(3);
% XF=Individual(2)+Individual(4);
% ZF=Individual(5);

ZA=Individual(1);
XB=Individual(2);
ZB=Individual(3);
XF=Individual(4);
ZF=Individual(5);

Truss_Data.dnkoor=[0	0	ZA;
XB	0	ZB;
0	XB	ZB;
-XB	0	ZB;
0	-XB	ZB;
XF	0	ZF;
XF*cos(45/180*pi)	XF*cos(45/180*pi)	ZF;
0	XF	ZF;
-XF*cos(45/180*pi)	XF*cos(45/180*pi)	ZF;
-XF	0	ZF;
-XF*cos(45/180*pi)	-XF*cos(45/180*pi)	ZF;
0	-XF	ZF;
XF*cos(45/180*pi)	-XF*cos(45/180*pi)	ZF;
6	0	0;
6*cos(45/180*pi)	6*cos(45/180*pi)	0;
0	6	0;
-6*cos(45/180*pi)	6*cos(45/180*pi)	0;
-6	0	0;
-6*cos(45/180*pi)	-6*cos(45/180*pi)	0;
0	-6	0;
6*cos(45/180*pi) -6*cos(45/180*pi) 0];

Truss_Data.A(1:4)=Individual(6);
Truss_Data.A(5:8)=Individual(7);
Truss_Data.A(9:16)=Individual(8);
Truss_Data.A(17:20)=Individual(9);
Truss_Data.A(21:28)=Individual(10);
Truss_Data.A(29:36)=Individual(11);
Truss_Data.A(37:44)=Individual(12);
Truss_Data.A(45:52)=Individual(13);

if ZB==ZF && XB==XF
    ObjVal=10^8;
    gg(1:2)=10^8;
    
else 
    
dnsay=size(Truss_Data.dnkoor,1);
elsay=size(Truss_Data.eldn,1);

topser=dnsay*3;
yer=zeros(topser,1);

[rijitlik,elboy]=stiffness_3D_truss(topser,elsay,Truss_Data.eldn,Truss_Data.dnkoor,Truss_Data.E,Truss_Data.A);
[mass_mat]=mass_truss(topser,elsay,Truss_Data.eldn,Truss_Data.dnkoor,Truss_Data.GS,Truss_Data.A);

for i=1:dnsay
mass_mat(3*i-2,3*i-2)=mass_mat(3*i-2,3*i-2)+Truss_Data.add_mass(i);
mass_mat(3*i-1,3*i-1)=mass_mat(3*i-1,3*i-1)+Truss_Data.add_mass(i);
mass_mat(3*i,3*i)=mass_mat(3*i,3*i)+Truss_Data.add_mass(i);
end

tutser=find(Truss_Data.MesKos'==1);
ser=setdiff([1:topser]',[tutser]);
[FI,LAMDA]=eig(inv(mass_mat(ser,ser))*rijitlik(ser,ser));
[LAMDA_sorted, ind]=sort(diag(LAMDA),'ascend');


Frekans=sqrt(LAMDA_sorted)/(2*pi);

%% Constraints

g(1)=(Frekans(1)/Truss_Data.limfre(1))-1;
gg(1)=g(1)*(g(1)>0);

g(2)=(Truss_Data.limfre(2)/Frekans(2))-1;
gg(2)=g(2)*(g(2)>0);

%% Object Function
ObjVal=0;
for i=1:size(Truss_Data.eldn,1)
ObjVal=ObjVal+elboy(i)*Truss_Data.A(i)*Truss_Data.GS;
end

end
%% Penalized Obj. Func.
PEN=10^8;
Z=ObjVal;
for k=1:length(gg)
    out = imag(gg(k)) ~= 0;
    if out == true
        gg(k)=abs(gg(k));
    end
     Z=Z+ PEN*gg(k);
end
TRUSS.PENALIZED=Z;