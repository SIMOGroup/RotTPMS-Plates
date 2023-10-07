function [D, I] = cal_Material_Matrices_2D_6dof_RotTPMS(Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate material matrices for RotTPMS with Quasi-3D 6-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from Plate
E_s = Plate.mat.E; nu_s = Plate.mat.nu; rho_s = Plate.mat.rho;
h = Plate.geo.h;
shear_func = Plate.theory.shear_func; stretch_func = Plate.theory.stretch_func;
porous_type = Plate.por_mat.type; RD  = Plate.por_mat.RD; alpha = Plate.por_mat.alpha;

%% ===== Initial material matrices =====
for i = 6:-1:1
    for j = 6:-1:1
        D{i,j} = zeros(6,6);
    end
end
for i = 4:-1:1
    for j = 4:-1:1
        I{i,j} = zeros(3,3);
    end
end

%% ===== Material properties =====
G_s = E_s/(2*(1+nu_s));

%% ===== Gauss integration =====
[Pg, Wg] = gauleg(1,-1,20) ;
z_upper = h/2; % FG upper
z_lower = -h/2; % FG lower

%% ===== Define material properties =====
% --- TPMS parameters ---
switch porous_type
    case {1, 2, 3}
        if porous_type == 1 % Primitive
            k_e = 0.25; C_1e = 0.317; n_1e = 1.264; n_2e = 2.006;
            k_g = 0.25; C_1g = 0.705; n_1g = 1.189; n_2g = 1.715;
            k_nu = 0.55; a_1 = 0.314; b_1 = -1.004; a_2 = 0.152;
        elseif porous_type == 2 % Gyroid
            k_e = 0.45; C_1e = 0.596; n_1e = 1.467; n_2e = 2.351;
            k_g = 0.45; C_1g = 0.777; n_1g = 1.544; n_2g = 1.982;
            k_nu = 0.50; a_1 = 0.192; b_1 = -1.349; a_2 = 0.402;
        elseif porous_type == 3 % IWP
            k_e = 0.35; C_1e = 0.597; n_1e = 1.225; n_2e = 1.782;
            k_g = 0.35; C_1g = 0.529; n_1g = 1.287; n_2g = 2.188;
            k_nu = 0.13; a_1 = 2.597; b_1 = -0.157; a_2 = 0.201;
        end
        C_2e = (C_1e*k_e^(n_1e) - 1)/(k_e^(n_2e) - 1); C_3e = 1 - C_2e;
        C_2g = (C_1g*k_g^(n_1g) - 1)/(k_g^(n_2g) - 1); C_3g = 1 - C_2g;
        d_1 = 0.3 - a_1*exp(b_1*k_nu); b_2 = - a_2*(k_nu + 1); d_2 = 0.3 - a_2*(1)^2 - b_2(1);
end

% --- Porous material properties ---
switch porous_type
    case {1, 2, 3}
        e = (RD <= k_e) * (C_1e*RD^(n_1e)) + ...
            (RD >  k_e) * (C_2e*RD.^(n_2e) + C_3e);
        g = (RD <= k_g) * (C_1g*RD^(n_1g)) + ...
            (RD >  k_g) * (C_2g*RD.^(n_2g) + C_3g);
        nu =(RD <= k_nu) * (a_1.*exp(b_1*RD) + d_1) + ...
            (RD >  k_nu) * (a_2*RD^(2) + b_2*RD + d_2);
        E = E_s*e;
        G = G_s*g;
    case 4
        E = E_s*((RD + 0.121)/1.121)^2.3;
        nu = 0.221*(1-RD) + nu_s*(0.342*(1-RD)^(2) - 1.21*(1-RD) + 1);
        G = E / (2*(1+nu));
    case 5
        E = E_s*(RD)^2;
        nu = 1/3;
        G = E /(2*(1+nu));
    case 6
        E = E_s*(RD)^2;
        nu = 0.3;
        G = E /(2*(1+nu));
end
rho = rho_s * RD;
C_11 = (E*(1-nu)) / ((1+nu)*(1-2*nu));
C_12 = (E*nu) / ((1+nu)*(1-2*nu));
C_44 = G;
C = [C_11, C_12, C_12,    0,    0,    0; ...
     C_12, C_11, C_12,    0,    0,    0; ...
     C_12, C_12, C_11,    0,    0,    0; ...    
        0,    0,    0, C_44,    0,    0; ...
        0,    0,    0,    0, C_44,    0; ...
        0,    0,    0,    0,    0, C_44];

% --- Rotating material constitutive matrix ---
alpha = deg2rad(alpha);
cx = cos(alpha(1)); sx = sin(alpha(1));
cy = cos(alpha(2)); sy = sin(alpha(2));
cz = cos(alpha(3)); sz = sin(alpha(3));

Rx = [1, 0, 0; 0, cx, -sx; 0, sx, cx]; % Counter-clockwise rotation
Ry = [cy, 0, sy; 0, 1, 0; -sy, 0, cy]; % Counter-clockwise rotation
Rz = [cz, -sz, 0; sz, cz, 0; 0, 0, 1]; % Counter-clockwise rotation
R = Rz*Ry*Rx;  % Rxyz order with Extrinsic rotation

T = [     R(1,1)^2,      R(1,2)^2,      R(1,3)^2,             2*R(1,2)*R(1,3),             2*R(1,3)*R(1,1),             2*R(1,1)*R(1,2); ...
          R(2,1)^2,      R(2,2)^2,      R(2,3)^2,             2*R(2,2)*R(2,3),             2*R(2,3)*R(2,1),             2*R(2,1)*R(2,2); ...
          R(3,1)^2,      R(3,2)^2,      R(3,3)^2,             2*R(3,2)*R(3,3),             2*R(3,3)*R(3,1),             2*R(3,1)*R(3,2); ...
     R(2,1)*R(3,1), R(2,2)*R(3,2), R(2,3)*R(3,3), R(2,2)*R(3,3)+R(2,3)*R(3,2), R(2,1)*R(3,3)+R(2,3)*R(3,1), R(2,2)*R(3,1)+R(2,1)*R(3,2); ...
     R(3,1)*R(1,1), R(3,2)*R(1,2), R(3,3)*R(1,3), R(1,2)*R(3,3)+R(1,3)*R(3,2), R(1,3)*R(3,1)+R(1,1)*R(3,3), R(1,1)*R(3,2)+R(1,2)*R(3,1); ...
     R(1,1)*R(2,1), R(1,2)*R(2,2), R(1,3)*R(2,3), R(1,2)*R(2,3)+R(1,3)*R(2,2), R(1,3)*R(2,1)+R(1,1)*R(2,3), R(1,1)*R(2,2)+R(1,2)*R(2,1)];
Q = T*C*T';

%% ===== Calculate material matrices =====
for iGauss = 1:size(Wg,1)  % Loop over the integration points
    % --- Gauss points & weights ---
    gpt = Pg(iGauss);
    gwt = Wg(iGauss);
    
    % --- Map the point to global ---
    gpt = (z_upper-z_lower)/2*gpt + (z_upper+z_lower)/2;
    gwt = (z_upper-z_lower)/2*gwt;
             
    % --- Shear deformation function ---
    [fz, dfz, ~] = compute_shear_deformation_function(gpt,h,shear_func);
    
    % --- Stretching effect deformation function ---
    [gz, dgz] = compute_stretching_deformation_function(gpt,h,shear_func,stretch_func);
    
    % --- General stiffness matrix ---
    h_e = [1 gpt fz dfz gz dgz];
    for i = 1:6
        for j = 1:6
            D{i,j} = D{i,j} + h_e(i)*h_e(j)*Q*gwt;
        end
    end
    
    % --- Inertia matrix ---
    h_u = [1 gpt fz gz];
    for i = 1:4
        for j = 1:4
            I{i,j} = I{i,j} + h_u(i)*h_u(j)*rho*eye(3)*gwt;
        end
    end
end
end
