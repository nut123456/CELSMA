% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________   
function s_rotate = rot_stress3D(Cx,Cy,Cz)
s_rotate = [-Cx -Cy -Cz Cx Cy Cz];
end