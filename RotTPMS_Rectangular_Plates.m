%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Title: Isogeometric analysis for RotTPMS plates with Quasi-3D six-variable plate model %%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% ! Please reference to paper: ............................................
% ! This work can be used, modified, and shared under the MIT License
% ! This work can be found in https://github.com/SIMOGroup/RotTPMS-Plates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========================== Initialization =============================
tic
addpath(genpath('./'));

% clc
clear all
% close all
format long

% global Plate IGA

%% ============================ Plate geometry ============================
% === Physical geometric properties ===
Plate.geo.L = 1; 
Plate.geo.W = 1;
Plate.geo.h = Plate.geo.L/(100);
% Plate.geo.h = Plate.geo.W/20;

% === Plate theory ===
% [1]: fz=0, [2]: Reddy, [3]: Shimpi, [4]: H. Nguyen-Xuan, [5]: Hoang X. Nguyen, [6]: Tuan N. Nguyen, [7]: Chien H. Thai
Plate.theory.shear_func = 5;
% [1]: f'z, [2]: 3/20*f'z, [3]: 1/8*f'z, [4]: 1/12*f'z, [5]: f'z + 1
Plate.theory.stretch_func = 2;

% === Define NURBS functions ===
IGA.NURBS.deg = 3; % Degree of basis functions
IGA.NURBS.ref = 21; % Number of mesh refinement

%% ============================ Material ==================================
% === Base material ===
% [1]: Steel, [2]: Alumium, [3]: Titanium, [4]: Copper, [5]: Brass, [6]: Steel 2
Plate.mat.type = 1;

% === Porous material ===
% [1]: Primitive, [2]: Gyroid, [3]: IWP, [4]: Closed-cell, [5]: Open-cell (\nu = 0.33), [6]: Mod Open-cell (\nu = 0.3)
Plate.por_mat.type = 1;
Plate.por_mat.RD = 1;
% Rxyz order with Extrinsic rotation, Counter-clockwise angles
Plate.por_mat.alpha = [0, 0, 0];  % [0, 0, 1] and Rz = 0
% Plate.por_mat.alpha = [0, 0, 45];  % [0, 0, 1] and Rz = 45
% Plate.por_mat.alpha = [90, 45, 0];  % [0, 1, 1] and Rz = 0
% Plate.por_mat.alpha = [45, -35.264, 25];  % [1, 1, 1] and Rz = 0
% Plate.por_mat.alpha = [90, -45, 0];  % [0, 1, 1] and Rz = 0

%% ========================== Problem type ================================
% [1]: Static, [2]: Vibration
Plate.prob.type = 1;
Plate.prob.q_uniform = -1;  % Static load

%% ========================= Boundary condition ===========================
% [1]: Fully simply supported (SSSS), [2]: Fully clamped (CCCC)
Plate.bc.bc_case = 2;

%% ======================== Material matrices =============================
[Plate.mat.E, Plate.mat.nu, Plate.mat.rho] = compute_material(Plate.mat.type);
Plate.mat.D = Plate.mat.E*Plate.geo.h^(3)/(12*(1-Plate.mat.nu^2));
[Plate.mat_mat.D, Plate.mat_mat.I] = cal_Material_Matrices_2D_6dof_RotTPMS(Plate);

%% =============================== IGA mesh ===============================
% === Generate NURBS mesh ===    
IGA.NURBS = Mesh_2D(Plate, IGA.NURBS);
IGA.NURBS = Gen_Ien_Inn_2D(IGA.NURBS);

% === NURBS properties ===
IGA.NURBS.nsd   = 2;                                                             % Number of spatial dimension
IGA.NURBS.nnode = IGA.NURBS.mcp * IGA.NURBS.ncp;                                 % Number of control point
IGA.NURBS.nshl  = (IGA.NURBS.p + 1) * (IGA.NURBS.q + 1);                         % Number of local shape functions (= degree + 1 per element due to k refinement)
IGA.NURBS.nel   = (IGA.NURBS.mcp - IGA.NURBS.p) * (IGA.NURBS.ncp - IGA.NURBS.q); % Number of element

% === IGA properties ===
IGA.params.ndof   = 6;                                                           % Number of dofs of a control point
IGA.params.sdof   = IGA.NURBS.nnode * IGA.params.ndof;                           % Total number of dofs of the structure
IGA.params.nGauss = IGA.NURBS.p + 1;                                             % Number of gauss point in integration

%% ========================= IGA for linear geometry ======================
% === Building global matrices ===
IGA.result.K = cal_Stiffness_Matrices_2D_6dof(IGA,Plate);    % Stiffness
IGA.result.M = cal_Mass_Matrices_2D_6dof(IGA,Plate);         % Mass
IGA.result.F = cal_Load_Vector_Uniform_2D_6dof(IGA,Plate);   % Extenal load

% === Imposing boundary conditions ===
[IGA.params.bcdof, IGA.params.bcval] = cal_bcdof_2D_6dof(IGA,Plate);
IGA.params.fdof = setdiff((1:IGA.params.sdof)', IGA.params.bcdof');  % Free dofs

% === Solving weak form ===
switch Plate.prob.type
    case 1  % Static
        bcdof = IGA.params.bcdof; bcval = IGA.params.bcval;
        sdof = IGA.params.sdof; fdof = IGA.params.fdof;
        IGA.result.U = zeros(sdof, 1); 
        IGA.result.U(bcdof') = bcval';
        IGA.result.F(fdof) = IGA.result.F(fdof) - IGA.result.K(fdof, bcdof')*bcval';
        IGA.result.U(fdof) = IGA.result.K(fdof, fdof) \ IGA.result.F(fdof);
        
        norm_method = 1;
        cen_def_norm = cal_central_deflection_2D_6dof(IGA,Plate,norm_method);
        sprintf("Normalized central deflection = " + sprintf('%.4f', cen_def_norm))
%         CP_deform = plot_deflection_2D_6dof(IGA,Plate,norm_method);
        clear sdof bcdof bcval fdof 
     case 2  % Free vibration
        bcdof = IGA.params.bcdof;
        KK = full(IGA.result.K); MM = full(IGA.result.M);
        [Lambda, ModeShape] = Eigen(KK,MM,bcdof);
        [IGA.result.Lambda, sort_index] = sort(Lambda,'ascend');
        IGA.result.ModeShape = ModeShape(:, sort_index);
        
        % --- Normalization ---
        norm_method = 1; nmode = 3;
        E_mat = Plate.mat.E; nu_mat = Plate.mat.nu; rho_mat = Plate.mat.rho; D_mat = Plate.mat.D;
        L = Plate.geo.L; h = Plate.geo.h;
        switch norm_method
            case 1
                Lambda_norm = (Lambda(1:nmode)*rho_mat*L^4*h/D_mat).^0.25;
            case 2
                norm = h*sqrt(rho_mat/E_mat);
                Lambda_norm = sqrt(Lambda(1:nmode))*norm;
            case 3
                norm = L*sqrt(rho_mat*(1-nu_mat^(2))/E_mat);
                Lambda_norm = sqrt(Lamda0(1:nmode))*norm;
            case 4
                norm = L^2/h*sqrt(rho_mat/E_mat);
                Lambda_norm = sqrt(Lamda0(1:nmode))*norm;
        end
        sprintf("Normalized frequencies = " + sprintf('%.4f, ', Lambda_norm))
        clear E_mat nu_mat rho_mat D_mat L h 
        clear bcdof KK MM Lambda ModeShape sort_index nmode
end
toc
