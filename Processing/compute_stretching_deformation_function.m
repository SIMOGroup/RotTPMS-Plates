function [gz, dgz] = compute_stretching_deformation_function(z,h,shear_func,stretch_func)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute stretching effect deformation function at point z %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%% === Shear deformation function ===
[fz, dfz, ddfz] = compute_shear_deformation_function(z,h,shear_func);

%% === Stretching deformation function ===
% [1]: f'z, [2]: 3/20*f'z, [3]: 1/8*f'z, [4]: 1/12*f'z, [5]: f'z + 1
switch stretch_func
    case 1
        gz = dfz;
        dgz = ddfz;
    case 2
        gz = 3/20*dfz;
        dgz = 3/20*ddfz;
    case 3
        gz = 1/8*dfz;
        dgz = 1/8*ddfz;
    case 4
        gz = 1/12*dfz;
        dgz = 1/12*ddfz;
    case 5
        gz = dfz + 1;
        dgz = ddfz;
    otherwise
        disp('Error')
        pause
end
end
