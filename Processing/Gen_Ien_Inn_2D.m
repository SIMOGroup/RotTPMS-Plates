function NURBS = Gen_Ien_Inn_2D(NURBS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate the IEN matrix, which relates element numbers and local node numbers to the appropriate global node numbers for 1D element; %%%
%%% Generates the INN matrix, which relates global node number to the "NURBS coordinates" of the node for 1D element. %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
p = NURBS.p; mcp = NURBS.mcp; 
q = NURBS.q; ncp = NURBS.ncp; 

%% ===== Initial index matrices =====
g = 0; e = 0;

%% ===== Compute Ien and Inn =====
for j = 1:ncp  % Loop through control points in V direction
    for i = 1:mcp  % Loop through control points in U direction
        g = g+1;
        Inn(g,1) = i;
        Inn(g,2) = j;
        if i >= p+1 && j >= q+1
            e = e + 1;
            for loop1 = 0:q
                for loop2 = 0:p
                    gtemp = g - mcp*loop1 - loop2;
                    ln = (p+1)*loop1 + loop2 + 1;
                    Ien(e,ln) = gtemp;
                end
            end
        end
    end
end

%% ===== Save NURBS =====
NURBS.Inn = Inn; NURBS.Ien = Ien;
end
