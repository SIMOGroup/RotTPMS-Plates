function [alpha_x, alpha_y] = compute_Rotation_Angles_xy(dir_vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute rotation angles alpha_x and alpha_y from Rotation Rxy mapping oirigin vector to [0;0;1] %%%
% Author: Kim Q. Tran
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
dir_vec = dir_vec / norm(dir_vec);
if ~ isempty(find(isnan(dir_vec), 1))
    dir_vec = double(isnan(dir_vec));
end
alpha_x = rad2deg(atan2(dir_vec(2), dir_vec(3)));
alpha_y = rad2deg(- atan2(dir_vec(1), sqrt(dir_vec(2)^2 + dir_vec(3)^2)));

end
