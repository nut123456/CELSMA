%Individual=[5.9098    ,  2.5359    ,  3.7213   ,   4.1294   ,   2.5001   ,   0.0001 , 0.00010485 , 0.00011884,  0.00014388 , 0.00014152  ,0.00010001 , 0.00015372 , 0.00016568];
function [TRUSS] = Truss_modal_analysis_52bar(Individual)
    Truss_Data = TrussData_52bar_modal();

    ZA = Individual(1);
    XB = Individual(2);
    ZB = Individual(3);
    XF = Individual(4);
    ZF = Individual(5);

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

    Truss_Data.A(1:4) = Individual(6);
    Truss_Data.A(5:8) = Individual(7);
    Truss_Data.A(9:16) = Individual(8);
    Truss_Data.A(17:20) = Individual(9);
    Truss_Data.A(21:28) = Individual(10);
    Truss_Data.A(29:36) = Individual(11);
    Truss_Data.A(37:44) = Individual(12);
    Truss_Data.A(45:52) = Individual(13);

    if ZB == ZF && XB == XF
        ObjVal = 10^8;
        gg(1:2) = 10^8;
    else
        dnsay = size(Truss_Data.dnkoor, 1);
        elsay = size(Truss_Data.eldn, 1);
        topser = dnsay * 3;
        yer = zeros(topser, 1);

        [rijitlik, elboy] = stiffness_3D_truss(topser, elsay, Truss_Data.eldn, Truss_Data.dnkoor, Truss_Data.E, Truss_Data.A);
        [mass_mat] = mass_truss(topser, elsay, Truss_Data.eldn, Truss_Data.dnkoor, Truss_Data.GS, Truss_Data.A);

        for i = 1:dnsay
            mass_mat(3 * i - 2, 3 * i - 2) = mass_mat(3 * i - 2, 3 * i - 2) + Truss_Data.add_mass(i);
            mass_mat(3 * i - 1, 3 * i - 1) = mass_mat(3 * i - 1, 3 * i - 1) + Truss_Data.add_mass(i);
            mass_mat(3 * i, 3 * i) = mass_mat(3 * i, 3 * i) + Truss_Data.add_mass(i);
        end

        tutser = find(Truss_Data.MesKos' == 1);
        ser = setdiff([1:topser]', [tutser]);
        [FI, LAMDA] = eig(inv(mass_mat(ser, ser)) * rijitlik(ser, ser));
        [LAMDA_sorted, ind] = sort(diag(LAMDA), 'ascend');

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
        if out == true
            gg(k) = abs(gg(k));
        end
        Z = Z + PEN * gg(k);
    end
    TRUSS.PENALIZED = Z;

    % Plot the undeformed and deformed shapes
    plot_dome_structure(Truss_Data.dnkoor, FI(:, 1));
end

function plot_dome_structure(nodes, displacement)
    figure;
    
    % Plot undeformed shape
    subplot(1, 2, 1);
    plot3(nodes(:, 1), nodes(:, 2), nodes(:, 3), 'bo-');
    title('Undeformed Shape');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    grid on;
    
    % Plot deformed shape
    subplot(1, 2, 2);
    scale = 10; % Scaling factor for visualization
    % Ensure displacement is reshaped correctly
    num_nodes = size(nodes, 1);
    % Verify displacement size
    if length(displacement) < num_nodes * 3
        displacement = [displacement; zeros(num_nodes * 3 - length(displacement), 1)];
    end
    displacement = displacement(1:num_nodes*3); % Ensure correct size
    deformed_nodes = nodes + scale * reshape(displacement, num_nodes, 3);
    plot3(deformed_nodes(:, 1), deformed_nodes(:, 2), deformed_nodes(:, 3), 'ro-');
    title('Deformed Shape');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    grid on;
end

function Truss_Data = TrussData_52bar_modal()
    % Modulus of Elasticity
    E = 2.10e11;

    % Bar connections between nodes 
    eldn = [1, 4;
            1, 5;
            1, 2;
            1, 3;
            4, 10;
            5, 12;
            2, 6;
            3, 8;
            4, 11;
            11, 5;
            5, 13;
            2, 13;
            2, 7;
            3, 7;
            3, 9;
            4, 9;
            4, 5;
            5, 2;
            2, 3;
            3, 4;
            10, 11;
            11, 12;
            12, 13;
            13, 6;
            6, 7;
            7, 8;
            8, 9;
            9, 10;
            10, 18;
            11, 19;
            12, 20;
            13, 21;
            6, 14;
            7, 15;
            8, 16;
            9, 17;
            10, 19;
            12, 19;
            12, 21;
            6, 21;
            6, 15;
            8, 15;
            8, 17;
            10, 17;
            11, 18;
            11, 20;
            13, 20;
            13, 14;
            7, 14;
            7, 16;
            9, 16;
            9, 18];

    % Support conditions
    MesKos = [0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              0, 0, 0;
              1, 1, 1;
              1, 1, 1;
              1, 1, 1;
              1, 1, 1;
              1, 1, 1;
              1, 1, 1;
              1, 1, 1;
              1, 1, 1];

    % Frequency limits (Hz)
    limfre(1) = 15.916;
    limfre(2) = 28.648;

    % Unit weight
    GS = 7800; % kg/m^3

    % added masses
    add_mass = zeros(size(MesKos, 1));
    add_mass(1:13) = 50; % kg

    Truss_Data = struct('MesKos', MesKos, 'eldn', eldn, 'E', E, 'GS', GS, 'add_mass', add_mass, 'limfre', limfre);
end

function [rijitlik, elboy] = stiffness_3D_truss(topser, elsay, eldn, dnkoor, E, A)
    rijitlik = zeros(topser);
    elboy = zeros(elsay, 1);

    for e = 1:elsay
        indis = eldn(e, :);
        elsd = [indis(1) * 3 - 2, indis(1) * 3 - 1, indis(1) * 3, indis(2) * 3 - 2, indis(2) * 3 - 1, indis(2) * 3];
        xa = dnkoor(indis(2), 1) - dnkoor(indis(1), 1);
        ya = dnkoor(indis(2), 2) - dnkoor(indis(1), 2);
        za = dnkoor(indis(2), 3) - dnkoor(indis(1), 3);
        elboy(e) = sqrt(xa^2 + ya^2 + za^2);
        CX = xa / elboy(e);
        CY = ya / elboy(e);
        CZ = za / elboy(e);
        DD = sqrt(CX^2 + CY^2);
        if CZ == 1
            t = [ 0, 0, 1, 0, 0, 0;
                  0, 1, 0, 0, 0, 0;
                 -1, 0, 0, 0, 0, 0;
                  0, 0, 0, 0, 0, 1;
                  0, 0, 0, 0, 1, 0;
                  0, 0, 0, -1, 0, 0];
        elseif CZ == -1
            t = [ 0, 0, -1, 0, 0, 0;
                  0, 1, 0, 0, 0, 0;
                  1, 0, 0, 0, 0, 0;
                  0, 0, 0, 0, 0, -1;
                  0, 0, 0, 0, 1, 0;
                  0, 0, 0, 1, 0, 0];
        else
            t = [ CX, CY, CZ, 0, 0, 0;
                 -CY / DD, CX / DD, 0, 0, 0, 0;
                 -CX * CZ / DD, -CY * CZ / DD, DD, 0, 0, 0;
                  0, 0, 0, CX, CY, CZ;
                  0, 0, 0, -CY / DD, CX / DD, 0;
                  0, 0, 0, -CX * CZ / DD, -CY * CZ / DD, DD];
        end
        stiff_el = E * A(e) / elboy(e) * [1, 0, 0, -1, 0, 0;
                                          0, 0, 0, 0, 0, 0;
                                          0, 0, 0, 0, 0, 0;
                                         -1, 0, 0, 1, 0, 0;
                                          0, 0, 0, 0, 0, 0;
                                          0, 0, 0, 0, 0, 0];
        stiff_glob = (t' * stiff_el) * t;
        rijitlik(elsd, elsd) = rijitlik(elsd, elsd) + stiff_glob;
    end
end

function [mass_mat] = mass_truss(topser, elsay, eldn, dnkoor, GS, A)
    mass_mat = zeros(topser);
    for e = 1:elsay
        indis = eldn(e, :);
        elsd = [indis(1) * 3 - 2, indis(1) * 3 - 1, indis(1) * 3, indis(2) * 3 - 2, indis(2) * 3 - 1, indis(2) * 3];
        xa = dnkoor(indis(2), 1) - dnkoor(indis(1), 1);
        ya = dnkoor(indis(2), 2) - dnkoor(indis(1), 2);
        za = dnkoor(indis(2), 3) - dnkoor(indis(1), 3);
        elboy = sqrt(xa^2 + ya^2 + za^2);
        CX = xa / elboy;
        CY = ya / elboy;
        CZ = za / elboy;
        DD = sqrt(CX^2 + CY^2);
        if CZ == 1
            t = [ 0, 0, 1, 0, 0, 0;
                  0, 1, 0, 0, 0, 0;
                 -1, 0, 0, 0, 0, 0;
                  0, 0, 0, 0, 0, 1;
                  0, 0, 0, 0, 1, 0;
                  0, 0, 0, -1, 0, 0];
        elseif CZ == -1
            t = [ 0, 0, -1, 0, 0, 0;
                  0, 1, 0, 0, 0, 0;
                  1, 0, 0, 0, 0, 0;
                  0, 0, 0, 0, 0, -1;
                  0, 0, 0, 0, 1, 0;
                  0, 0, 0, 1, 0, 0];
        else
            t = [ CX, CY, CZ, 0, 0, 0;
                 -CY / DD, CX / DD, 0, 0, 0, 0;
                 -CX * CZ / DD, -CY * CZ / DD, DD, 0, 0, 0;
                  0, 0, 0, CX, CY, CZ;
                  0, 0, 0, -CY / DD, CX / DD, 0;
                  0, 0, 0, -CX * CZ / DD, -CY * CZ / DD, DD];
        end

        mass = 1 / 2 * GS * A(e) * elboy * [1, 0, 0, 0, 0, 0;
                                            0, 1, 0, 0, 0, 0;
                                            0, 0, 1, 0, 0, 0;
                                            0, 0, 0, 1, 0, 0;
                                            0, 0, 0, 0, 1, 0;
                                            0, 0, 0, 0, 0, 1];

        mass_glo = (t' * mass) * t;

        mass_mat(elsd, elsd) = mass_mat(elsd, elsd) + mass_glo;
    end
end

% Example call
Individual = [5.9098, 2.5359, 3.7213, 4.1294, 2.5001, 0.0001, 0.00010485, 0.00011884, 0.00014388, 0.00014152, 0.00010001, 0.00015372, 0.00016568];
Truss_modal_analysis_52bar(Individual);