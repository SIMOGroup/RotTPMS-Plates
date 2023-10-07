function [alpha_z] = compute_Rotation_Angle_z(alpha_x, alpha_y, rot_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute rotation angle alpha_z from Rotation Rxy mapping oirigin vector to [0;0;1] %%%
% Author: Kim Q. Tran
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reversed angles to find origin vector
alpha_rev = - deg2rad([alpha_x, alpha_y]);
cx_rev = cos(alpha_rev(1)); sx_rev = sin(alpha_rev(1));
cy_rev = cos(alpha_rev(2)); sy_rev = sin(alpha_rev(2));

Rx = [1, 0, 0; 0, cx_rev, -sx_rev; 0, sx_rev, cx_rev]; % Counter-clockwise rotation
Ry = [cy_rev, 0, sy_rev; 0, 1, 0; -sy_rev, 0, cy_rev]; % Counter-clockwise rotation
R_xy = Rx*Ry;  % Extrinsic % Reversed order
Ori_vec = R_xy*[0; 0; 1]; 

%% 
alpha = deg2rad([alpha_x, alpha_y]);
cx = cos(alpha(1)); sx = sin(alpha(1));
cy = cos(alpha(2)); sy = sin(alpha(2));

%% Calculate rotation angle alpha_z
% [1]: X' _|_ OYZ, [2]: [1;1;1] _|_ OYZ, [3]: Depends on subspaces

switch rot_method
    case 1
        alpha_z = 0;
    case 2
        alpha_z = rad2deg(atan2(cy + cx*sy + sx*sy, cx - sx));
    case 3
        if prod(Ori_vec) == 0
            if Ori_vec(1) == 0
                subspace = 1;  % Subspace YZ
            elseif Ori_vec(2) == 0
                subspace = 2;  % Subspace XZ
            elseif Ori_vec(3) == 0
                subspace = 3;  % Subspace XY
            end
        else
            yx = Ori_vec(2)/Ori_vec(1); xz = Ori_vec(1)/Ori_vec(3); zy = Ori_vec(3)/Ori_vec(2);
            if yx >= 1 && xz < 1
                subspace = 1;  % Subspace YZ
            elseif yx < 1 && zy >= 1
                subspace = 2;  % Subspace XZ
            elseif xz >= 1 && zy < 1
                subspace = 3;  % Subspace XY
            end
        end

        if subspace == 1  % Subspace YZ: X' _|_ OYZ
            alpha_z = rad2deg(-atan2(0,cy));
        elseif subspace == 2  % Subspace XZ: Y' _|_ OYZ
            alpha_z = rad2deg(-atan2(cx,sx*sy));
        elseif subspace == 3  % Subspace XY: Z' _|_ OYZ
            alpha_z = rad2deg(-atan2(-sx,cx*sy));
        end
end
end
