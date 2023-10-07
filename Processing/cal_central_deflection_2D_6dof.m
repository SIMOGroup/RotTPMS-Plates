function [cen_def_norm] = cal_central_deflection_2D_6dof(IGA,Plate,norm_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate central deflection of rectangular plate %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
uKnot = IGA.NURBS.uKnot; vKnot = IGA.NURBS.vKnot;
Inn = IGA.NURBS.Inn; Ien = IGA.NURBS.Ien; gcoord = IGA.NURBS.gcoord; nel = IGA.NURBS.nel;
ndof = IGA.params.ndof;
U = IGA.result.U;

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
shear_func = Plate.theory.shear_func; stretch_func = Plate.theory.stretch_func;
q_uniform = Plate.prob.q_uniform;
D_mat = Plate.mat.D;

%% ===== Central deflection =====
% --- Initial central element and point ---
c_ele = floor((nel+1)/2); % Central element
c_point = [L/2,W/2]; % Physical coordinate of the center point

% --- Element parameters ---
sctr = Ien(c_ele,:);  % Control points indexes
ni = Inn(Ien(c_ele,1),1);  % Index of the element in parametric domain
nj = Inn(Ien(c_ele,1),2);
nn = length(sctr);  % Number of control points in the element
for idof = 1:ndof  % Dofs of control points
    sctrF(idof:ndof:ndof*(nn-1) + idof) = ndof.*(sctr-1) + idof;
end
nodes = gcoord(sctr,:);

% --- Point parameters ---
nat_coord = cal_xi_eta(IGA,ni,nj,nodes,c_point);  % Natural coordinate
nat_coord(1) = (uKnot(ni+1)-uKnot(ni))/2*nat_coord(1) + (uKnot(ni+1)+uKnot(ni))/2;  % Map the point to parametric domain
nat_coord(2) = (vKnot(nj+1)-vKnot(nj))/2*nat_coord(2) + (vKnot(nj+1)+vKnot(nj))/2;

% --- Kinematic matrices of NURBS shape function ---
[N, dNdxi, dNdxi2, dNdxy, dN2dxy, detJ1] = Kine_Shape_2nd_2D(IGA,ni,nj,nat_coord(1),nat_coord(2));

% --- Stretching effect deformation function ---
[gz, ~] = compute_stretching_deformation_function(0,h,shear_func,stretch_func);
                
% --- Deflection ---
N0 = zeros(1,ndof*nn);
N0(1,3:ndof:ndof*nn) = N';
N0(1,6:ndof:ndof*nn) = gz*N';
cen_def = N0*U(sctrF);

% --- Normalization ---
switch norm_method
    case 1
        norm = 100*D_mat/(q_uniform*L^4);
        cen_def_norm = cen_def*norm;
    case 2
        norm = 1/h;
        cen_def_norm = cen_def*norm;
end
end
