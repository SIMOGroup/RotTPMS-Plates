function [nat_coord] = cal_xi_eta(IGA,ni,nj,nodes,c_point)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate natural coordinate (xi, eta) of a point %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
uKnot = IGA.NURBS.uKnot; vKnot = IGA.NURBS.vKnot;

%% ===== Initial parametric coordinate =====
nat_coord = zeros(1,2);
xi = 0; eta = 0;

%% ===== Find coordinate =====
% --- Loop to find the parametric parameter ---
EPS = 1e-10; err = inf;
while (err > EPS)
    % Map the point to parametric domain
    nat_coord(1) = (uKnot(ni+1)-uKnot(ni))/2*nat_coord(1) + (uKnot(ni+1)+uKnot(ni))/2;
    nat_coord(2) = (vKnot(nj+1)-vKnot(nj))/2*nat_coord(2) + (vKnot(nj+1)+vKnot(nj))/2;
    
    % Kinematic matrices of NURBS shape function
    [N, dNdxi, dNdxi2, dNdxy, dN2dxy, detJ1] = Kine_Shape_2nd_2D(IGA,ni,nj,nat_coord(1),nat_coord(2));
    
    % Error in physical domain
    x = N'*nodes(:,1); y = N'*nodes(:,2);
    deltaX = x - c_point(1); deltaY = y - c_point(2);
    delta = [deltaX; deltaY];
    
    % Derivatives
    dxdxi = dNdxi(:,1)' * nodes(:,1);
    dxdeta = dNdxi(:,2)' * nodes(:,1);
    dydxi = dNdxi(:,1)' * nodes(:,2);
    dydeta = dNdxi(:,2)' * nodes(:,2);
    F = [dxdxi, dxdeta;...
         dydxi, dydeta];
    invF = inv(F);

    % Compute new natural coordinate for next iteration
    xi = xi - invF(1,:)*delta; eta = eta - invF(2,:)*delta;
    nat_coord(1) = xi; nat_coord(2) = eta;
    
    err = sqrt((deltaX)^2+(deltaY)^2);
end
