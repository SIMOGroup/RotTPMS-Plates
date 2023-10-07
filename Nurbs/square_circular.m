% -------------------------------------------------------------------------
%                                                                        
%  NURBS based, single patch, 3D linear elasticity code               
%  Main routine to call all the subroutines                        
%  June 17, 2003
%  J. Austin Cottrell
%  CES Graduate Student
%  Texas Institute for Computational Engineering Science
%  University of Texas at Austin
%
%  Modified to 2D linear elasticity using Matlab by :
%  Hung Nguyen Xuan 
%  Date: Oct 2009
%  Faculty of Mathematics & Informatics, University of Natural Sciences
%  Vietnam   National University?HCM
%--------------------------------------------------------------------------
% Chu y: do cac cach goi trong set path ve data la tuong tu, nen khi chay
% bai toan nao thi remove phan data trong set path cua cac bai toan khac
clear all;
format long;
addpath Nurbs
addpath Processing
addpath Boundary

deg=2; ref=1;  %%% refinement


p=2;
q=1;
U = [0 0 0 1/2 1/2 1 1 1]; %m=5
V = [0 0 1 1]; %n=2

R=4;
CP(:,:,1)=[0 R;0 2*R;0 2*R;0 2*R;0 R];
CP(:,:,2)=[R R;R/2 R;0 0;-R/2 -R;-R -R];
CP(:,:,3)=[1 1;1 1;1 1;1 1;1 1]*5;
CP(:,:,4)=[1 1;1 1/sqrt(2); 1 1;1 1/sqrt(2); 1 1];


%=====================================================================
% REFINE
[CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);
R1 = refinement_vec_repeated_p2(U,ref);
R2 = refinement_vec_repeated_p2(V,ref);

[CP,U,V] = Knot_Refine_Surf(p,q,U,V,CP,R1,R2);
plotNURBS_surf_El_CP(p,q,U,V,CP); hold on

%view(2)
%break
% section computational IgA
p=deg;
q=p;
ngauss=p+1;                 % number of gauss point in integratio

%CP=importdata('P2_15.mat');
b_net(:,:,1)=CP(:,:,1);
b_net(:,:,2)=CP(:,:,2);
b_net(:,:,3)=CP(:,:,4);

mcp= length(b_net(:,1,1));
ncp= length(b_net(1,:,1));
gcoord(:,1)=reshape(b_net(:,:,1),mcp*ncp,1);
gcoord(:,2)=reshape(b_net(:,:,2),mcp*ncp,1);
gcoord(:,3)=reshape(b_net(:,:,3),mcp*ncp,1);
%clear CP U V R1 R2

% generate connectivities     
[ien,inn]=genIEN_INN_2D_repeated_Ck(p,q,mcp,ncp);


nodes=connecElementQ4(mcp-1,ncp-1);

for iel=1:size(nodes,1)
    nd=nodes(iel,:);
    x=gcoord(nd,1);
    y=gcoord(nd,2);    
%    xx=mean(x); 
%    yy=mean(y); 
     patch(x,y,'w')
     text(x(1),y(1),num2str(nd(1)));
     text(x(2),y(2),num2str(nd(2)));
     text(x(3),y(3),num2str(nd(3)));
     text(x(4),y(4),num2str(nd(4)));
end


