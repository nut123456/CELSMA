function Truss_Data=TrussData_52bar_modal

% Modulus of Elasticity
E=2.10e11;

% Bar connections between nodes 
eldn=[1	4	;
1	5	;
1	2	;
1	3	;
4	10	;
5	12	;
2	6	;
3	8	;
4	11	;
11	5	;
5	13	;
2	13	;
2	7	;
3	7	;
3	9	;
4	9	;
4	5	;
5	2	;
2	3	;
3	4	;
10	11	;
11	12	;
12	13	;
13	6	;
6	7	;
7	8	;
8	9	;
9	10	;
10	18	;
11	19	;
12	20	;
13	21	;
6	14	;
7	15	;
8	16	;
9	17	;
10	19	;
12	19	;
12	21	;
6	21	;
6	15	;
8	15	;
8	17	;
10	17	;
11	18	;
11	20	;
13	20	;
13	14	;
7	14	;
7	16	;
9	16	;
9	18	];

% Support conditions
MesKos=[0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
0	0	0	;
1	1	1	;
1	1	1	;
1	1	1	;
1	1	1	;
1	1	1	;
1	1	1	;
1	1	1	;
1	1	1	];

% Frequency limits (Hz)
limfre(1)=15.916;
limfre(2)=28.648;

% Unit weight
GS=7800;%kg/m^3

% added masses
add_mass=zeros(size(MesKos,1));
add_mass(1:13)=50; %kg

Truss_Data=struct('MesKos',MesKos,'eldn',eldn,'E',E,'GS',GS,'add_mass',add_mass,'limfre',limfre);

