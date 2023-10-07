function [R, dRdxi, dRdxi2, dRdx, dRdx2, detJ] = Kine_Shape_2nd_2D(IGA,ni,nj,u,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate non-zero basis functions, the matrix of gradients of all non-zero basis functions wrt parameter variable xi and physical variable %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
p = IGA.NURBS.p; mcp = IGA.NURBS.mcp; uKnot = IGA.NURBS.uKnot; 
q = IGA.NURBS.q; ncp = IGA.NURBS.ncp; vKnot = IGA.NURBS.vKnot;
B_net = IGA.NURBS.B_net; 
nsd = IGA.NURBS.nsd; nshl = IGA.NURBS.nshl;

%% ===== Initial Kinematic matrices =====
R = zeros(nshl,1);
dRdxi = zeros(nshl,nsd);
dRdxi2 = zeros(nshl,nsd+1);
denom_sum = 0; derv_sum_u = 0; derv_sum_v = 0;
derv_sum_uu = 0; derv_sum_vv = 0; derv_sum_uv = 0;
dxdxi = zeros(nsd,nsd);
dxdxi2 = zeros(nsd,nsd);

%% =====  Basis functions and their derivatives =====
M = dersbasisfuns(ni,p,mcp,u,uKnot) ;
N = dersbasisfuns(nj,q,ncp,v,vKnot) ;
% evaluate 1d size functions and derivatives each direction size and M and N = number of derviatives+shape functions, degree of poynomial.
% row 1 of M and N => shape functions i^{th} row (i > 1) => i^{th} derivative of the shape function calculate in u direction

%% ===== NURBS functions and their derivatives on parametric domain (R and dR/dxi) =====
icount = 0;
for j = 0:q
    for i = 0:p
        icount = icount+1;
        
        % Basis functions
        R(icount,1) = M(1,p+1-i)*N(1,q+1-j)*B_net(ni-i,nj-j,nsd+1);
        denom_sum = denom_sum + R(icount);
        
        % First derivatives
        dRdxi(icount,1) = M(2,p+1-i)*N(1,q+1-j)*B_net(ni-i,nj-j,nsd+1);  % du
        derv_sum_u = derv_sum_u + dRdxi(icount,1);
        dRdxi(icount,2) = M(1,p+1-i)*N(2,q+1-j)*B_net(ni-i,nj-j,nsd+1);  % dv
        derv_sum_v = derv_sum_v + dRdxi(icount,2);
        
        % Second derivatives
        dRdxi2(icount,1) = M(3,p+1-i)*N(1,q+1-j)*B_net(ni-i,nj-j,nsd+1);  % du2
        derv_sum_uu = derv_sum_uu + dRdxi2(icount,1) ;
        dRdxi2(icount,2) = M(1,p+1-i)*N(3,q+1-j)*B_net(ni-i,nj-j,nsd+1) ;  % dv2
        derv_sum_vv = derv_sum_vv + dRdxi2(icount,2) ;
        dRdxi2(icount,3) = M(2,p+1-i)*N(2,q+1-j)*B_net(ni-i,nj-j,nsd+1);  %duv
        derv_sum_uv = derv_sum_uv + dRdxi2(icount,3) ;
    end
end

% Basis functions
R = R/denom_sum;

% First derivatives ... divide through by denominator
dRdxi(:,1) = dRdxi(:,1)/denom_sum -(R(:)*derv_sum_u)/(denom_sum^2);
dRdxi(:,2) = dRdxi(:,2)/denom_sum -(R(:)*derv_sum_v)/(denom_sum^2);

% Second derivatives ... divide through by denominator
dRdxi2(:,1) = dRdxi2(:,1)/denom_sum - 2*derv_sum_u*dRdxi(:,1)/denom_sum^2 + ...
              - derv_sum_uu*R(:)/denom_sum^2+ 2*derv_sum_u^2*R(:)/denom_sum^3 ;
dRdxi2(:,2) = dRdxi2(:,2)/denom_sum - 2*derv_sum_v*dRdxi(:,2)/denom_sum^2 + ...
              - derv_sum_vv*R(:)/denom_sum^2+ 2*derv_sum_v^2*R(:)/denom_sum^3 ;
dRdxi2(:,3) = dRdxi2(:,3)/denom_sum - dRdxi(:,1)*derv_sum_v/denom_sum^2 + ...
              - dRdxi(:,2)*derv_sum_u/denom_sum^2 + ... 
              - R(:)*derv_sum_uv/denom_sum^2 + 2*R(:)*derv_sum_u*derv_sum_v/denom_sum^3;

%% ===== Derivative of physical domain wrt parametric domain (dxdxi, dxdxi2) =====
icount = 0;
for j = 0:q
    for i = 0:p
        icount = icount + 1;
        dxdxi(1,1) = dxdxi(1,1) + B_net(ni-i,nj-j,1)*dRdxi(icount,1);
        dxdxi(1,2) = dxdxi(1,2) + B_net(ni-i,nj-j,1)*dRdxi(icount,2);
        dxdxi(2,1) = dxdxi(2,1) + B_net(ni-i,nj-j,2)*dRdxi(icount,1);
        dxdxi(2,2) = dxdxi(2,2) + B_net(ni-i,nj-j,2)*dRdxi(icount,2);
    end
end

%% ===== NURBS functions and their derivatives on physical domain (R and dR/dx) =====
dxidx = inv(dxdxi); detJ = det(dxdxi);
dRdx = dRdxi*dxidx ;

% Note that DetJ resides in common
%  if(detj < 0) % chu y doan code doi dau detj nay khi can thiet.
%  detj = -detj;
%  end

% for higher order derivatives
icount = 0 ;
jac1 = zeros(3,2);
for j = 0:q
    for i = 0:p
        icount = icount + 1 ;
        jac1(1,1) = jac1(1,1) + B_net(ni-i,nj-j,1)*dRdxi2(icount,1);
        jac1(1,2) = jac1(1,2) + B_net(ni-i,nj-j,2)*dRdxi2(icount,1);
        jac1(2,1) = jac1(2,1) + B_net(ni-i,nj-j,1)*dRdxi2(icount,2);
        jac1(2,2) = jac1(2,2) + B_net(ni-i,nj-j,2)*dRdxi2(icount,2);
        jac1(3,1) = jac1(3,1) + B_net(ni-i,nj-j,1)*dRdxi2(icount,3);
        jac1(3,2) = jac1(3,2) + B_net(ni-i,nj-j,2)*dRdxi2(icount,3);
    end
end
jac2 = [         dxdxi(1,1)^2,          dxdxi(2,1)^2,                     2*dxdxi(1,1)*dxdxi(2,1); ...
                 dxdxi(1,2)^2,          dxdxi(2,2)^2,                     2*dxdxi(1,2)*dxdxi(2,2); ...
        dxdxi(1,1)*dxdxi(1,2), dxdxi(2,1)*dxdxi(2,2), dxdxi(2,1)*dxdxi(1,2)+dxdxi(1,1)*dxdxi(2,2)];
  
term1 = dRdx*jac1';
term2 = dRdxi2 - term1;
dxidx2 = inv(jac2)*term2';
dRdx2 = dxidx2';

end