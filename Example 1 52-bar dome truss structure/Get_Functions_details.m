%% Benchmark Test functions
function [lb, ub, dim, fobj] = Get_Functions_details(F)
    switch F
        case 'F1'
            fobj = @F1;
            lb = [4, 0.01, 3.7, 2.01, 2.5, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001];
            ub = [8, 4,    7.7,    6, 6.5, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
            dim = length(lb);
    end
end

function W = F1(Individual)
    Truss_Data = TrussData_52bar_modal;

    % Extract individual components
    ZA = Individual(1);
    XB = Individual(2);
    ZB = Individual(3);
    XF = Individual(4);
    ZF = Individual(5);

    % Define truss node coordinates
    Truss_Data.dnkoor = [0, 0, ZA;
                         XB, 0, ZB;
                         0, XB, ZB;
                         -XB, 0, ZB;
                         0, -XB, ZB;
                         XF, 0, ZF;
                         XF * cos(45 / 180 * pi), XF * cos(45 / 180 * pi), ZF;
                         0, XF, ZF;
                         -XF * cos(45 / 180 * pi), XF * cos(45 / 180 * pi), ZF;
                         -XF, 0, ZF;
                         -XF * cos(45 / 180 * pi), -XF * cos(45 / 180 * pi), ZF;
                         0, -XF, ZF;
                         XF * cos(45 / 180 * pi), -XF * cos(45 / 180 * pi), ZF;
                         6, 0, 0;
                         6 * cos(45 / 180 * pi), 6 * cos(45 / 180 * pi), 0;
                         0, 6, 0;
                         -6 * cos(45 / 180 * pi), 6 * cos(45 / 180 * pi), 0;
                         -6, 0, 0;
                         -6 * cos(45 / 180 * pi), -6 * cos(45 / 180 * pi), 0;
                         0, -6, 0;
                         6 * cos(45 / 180 * pi), -6 * cos(45 / 180 * pi), 0];

    % Assign individual components to truss areas
    Truss_Data.A(1:4) = Individual(6);
    Truss_Data.A(5:8) = Individual(7);
    Truss_Data.A(9:16) = Individual(8);
    Truss_Data.A(17:20) = Individual(9);
    Truss_Data.A(21:28) = Individual(10);
    Truss_Data.A(29:36) = Individual(11);
    Truss_Data.A(37:44) = Individual(12);
    Truss_Data.A(45:52) = Individual(13);

    % Check constraints
    if ZB == ZF && XB == XF
        ObjVal = 10^8;
        gg(1:2) = 10^8;
    else
        % Initialize variables
        dnsay = size(Truss_Data.dnkoor, 1);
        elsay = size(Truss_Data.eldn, 1);
        topser = dnsay * 3;
        yer = zeros(topser, 1);

        % Calculate stiffness and mass matrices
        [rijitlik, elboy] = stiffness_3D_truss(topser, elsay, Truss_Data.eldn, Truss_Data.dnkoor, Truss_Data.E, Truss_Data.A);
        [mass_mat] = mass_truss(topser, elsay, Truss_Data.eldn, Truss_Data.dnkoor, Truss_Data.GS, Truss_Data.A);

        % Add additional mass
        for i = 1:dnsay
            mass_mat(3 * i - 2, 3 * i - 2) = mass_mat(3 * i - 2, 3 * i - 2) + Truss_Data.add_mass(i);
            mass_mat(3 * i - 1, 3 * i - 1) = mass_mat(3 * i - 1, 3 * i - 1) + Truss_Data.add_mass(i);
            mass_mat(3 * i, 3 * i) = mass_mat(3 * i, 3 * i) + Truss_Data.add_mass(i);
        end

        % Constraint calculations
        tutser = find(Truss_Data.MesKos' == 1);
        ser = setdiff([1:topser]', [tutser]);
        [FI, LAMDA] = eig(inv(mass_mat(ser, ser)) * rijitlik(ser, ser));
        [LAMDA_sorted, ind] = sort(diag(LAMDA), 'ascend');

        % Frequency calculation
        Frekans = sqrt(LAMDA_sorted) / (2 * pi);

        % Constraints
        g(1) = (Frekans(1) / Truss_Data.limfre(1)) - 1;
        gg(1) = g(1) * (g(1) > 0);

        g(2) = (Truss_Data.limfre(2) / Frekans(2)) - 1;
        gg(2) = g(2) * (g(2) > 0);

        % Objective Function
        ObjVal = 0;
        for i = 1:size(Truss_Data.eldn, 1)
            ObjVal = ObjVal + elboy(i) * Truss_Data.A(i) * Truss_Data.GS;
        end
    end

    % Penalized Objective Function
    PEN = 10^8;
    Z = ObjVal;
    for k = 1:length(gg)
        out = imag(gg(k)) ~= 0;
        if out
            gg(k) = abs(gg(k));
        end
        Z = Z + PEN * gg(k);
    end
    
    W = Z; % Assign value to W to ensure output is returned correctly
end
%end