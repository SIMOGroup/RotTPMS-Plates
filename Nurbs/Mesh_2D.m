function NURBS = Mesh_2D(Plate, NURBS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Isogeometric mesh for 2D element %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
deg = NURBS.deg; ref = NURBS.ref; 

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W;

%% ===== Square mesh generation ======
[CP, uKnot, vKnot, p, q] = Square_Coarse_Mesh(L,W,deg);
R1 = Refinement_vec(uKnot,ref); R2 = Refinement_vec(vKnot,ref); % Refinement
[CP, uKnot, vKnot] = Knot_Refine_Surf(p,q,uKnot,vKnot,CP,R1,R2); % K Refinement

%% ====== Reduced to 2D mesh (xy plane) ======
mcp = length(CP(:,1,1)); ncp = length(CP(1,:,1));
B_net(:,:,1)= CP(:,:,1); B_net(:,:,2)= CP(:,:,2); B_net(:,:,3)= CP(:,:,4);

%% ===== Global Coordinate of Control Points ======
gcoord(:,1) = reshape(B_net(:,:,1),mcp*ncp,1);
gcoord(:,2) = reshape(B_net(:,:,2),mcp*ncp,1);
gcoord(:,3) = reshape(B_net(:,:,3),mcp*ncp,1);

%% ===== Save NURBS =====
NURBS.p = p; NURBS.uKnot = uKnot; NURBS.mcp = mcp;
NURBS.q = q; NURBS.vKnot = vKnot; NURBS.ncp = ncp;
NURBS.CP = CP; NURBS.B_net = B_net; NURBS.gcoord = gcoord;
NURBS.deg = deg; NURBS.ref = ref; 

%% ====== Plot Mesh ======
plot_option = "No";
if plot_option == "Yes"
%     set(gcf,'color','white')
    plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
    plot_ctrlnet(CP,'bo'); hold on
    % Show node index
    count = 0;
    for i = 1:size(CP,1)
       for j = 1:size(CP,2) 
           count = count +1;
           X = CP(i,j,1); Y = CP(i,j,2);
           text(X,Y,num2str(count));
       end
    end
end
end
