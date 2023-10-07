function [CP_deform] = plot_deflection_2D_6dof(IGA,Plate,norm_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot deflection of plate %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
p = IGA.NURBS.p; uKnot = IGA.NURBS.uKnot; 
q = IGA.NURBS.q; vKnot = IGA.NURBS.vKnot;
CP = IGA.NURBS.CP;
ndof = IGA.params.ndof; 
U = IGA.result.U;

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
shear_func = Plate.theory.shear_func; stretch_func = Plate.theory.stretch_func;
q_uniform = Plate.prob.q_uniform;
D_mat = Plate.mat.D;

%% ===== Initial deformation plot =====
CP_deform = CP;

%% ====== Plot deflection ======
size_factor = 1;
count = 0;
for i = 1:size(CP,1)
    for j = 1:size(CP,2) 
        count = count + 1;
        
        % --- Stretching effect deformation function ---
        [gz, ~] = compute_stretching_deformation_function(0,h,shear_func,stretch_func);
        def = U(ndof*(count-1)+3) + gz*U(ndof*(count-1)+6);
        
        % --- Normalization ---
        switch norm_method
            case 1
                norm = 100*D_mat/(q_uniform*L^4);
                def_norm = def*norm;
            case 2
                norm = 1/h;
                def_norm = def*norm;
        end
        CP_deform(i,j,3) = CP(i,j,3) + size_factor*sign(q_uniform)*def_norm;
    end
end
plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP_deform);
end
