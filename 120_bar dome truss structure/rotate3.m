% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________  
function rot_mat = rotate3(Cx,Cy,Cz)
rot_mat = zeros(6,6);
rot_mat(1,1) = Cx^2;
rot_mat(1,2) = Cx*Cy;
rot_mat(2,1) = rot_mat(1,2);
rot_mat(1,3) = Cx*Cz;
rot_mat(3,1) = rot_mat(1,3);
rot_mat(1,4) = -Cx^2;
rot_mat(4,1) = rot_mat(1,4);
rot_mat(1,5) = -Cx*Cy;
rot_mat(5,1) = rot_mat(1,5);
rot_mat(1,6) = -Cx*Cz;
rot_mat(6,1) = rot_mat(1,6);
rot_mat(2,2) = Cy^2;
rot_mat(2,3) = Cy*Cz;
rot_mat(3,2) = rot_mat(2,3);
rot_mat(2,4) = -Cx*Cy;
rot_mat(4,2) = rot_mat(2,4);
rot_mat(2,5) = -Cy^2;
rot_mat(5,2) = rot_mat(2,5);
rot_mat(2,6) = -Cy*Cz;
rot_mat(6,2) = rot_mat(2,6);
rot_mat(3,3) = Cz^2;
rot_mat(3,4) = -Cx*Cz;
rot_mat(4,3) = rot_mat(3,4);
rot_mat(3,5) = -Cy*Cz;
rot_mat(5,3) = rot_mat(3,5);
rot_mat(3,6) = -Cz*Cz;
rot_mat(6,3) = rot_mat(3,6);
rot_mat(4,4) = Cx^2;
rot_mat(4,5) = Cx*Cy;
rot_mat(5,4) = rot_mat(4,5);
rot_mat(4,6) = Cx*Cz;
rot_mat(6,4) = rot_mat(4,6);
rot_mat(5,5) = Cy^2;
rot_mat(5,6) = Cy*Cz;
rot_mat(6,5) = rot_mat(5,6);
rot_mat(6,6) = Cz^2;
end