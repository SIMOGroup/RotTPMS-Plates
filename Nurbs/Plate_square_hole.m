clear all
%clc
clf
format long

deg=2; ref=11;  %%% refinement

%plate circle 
p=1;
q=1;
mcp=8;
U = [zeros(1,p) linspace(0,1,mcp-p+1) ones(1,p)];

%U = [0 0 1/4 1/4 2/4 2/4 3/4 3/4 1 1]; %m=9
V = [0 0 1 1]; %n=2
%weights= [ 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1]';
R=2;
r=R*0.6;
CP(:,:,1)=[R r;-R -r;-R -r;-R -r;-R -r;R r;R r;R r];
CP(:,:,2)=[R r;R r;R r;-R -r;-R -r;-R -r;-R -r;R r];
CP(:,:,3)=[1 1; 1 1; 1 1; 1 1 ;1 1;1 1; 1 1; 1 1]*5;
CP(:,:,4)=[ 1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1];

%CP(:,:,4)=[1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2)];


% Gauss Points
ngauss = [deg+1 deg+1];

%=====================================================================
% REFINE
[CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);
R1 = refinement_vec_repeated_p2(U,ref);
R2 = refinement_vec_repeated_p2(V,ref);

[CP,U,V] = Knot_Refine_Surf(p,q,U,V,CP,R1,R2);
plotNURBS_surf_El_CP(p,q,U,V,CP); hold on
size(CP)
view(2)
% %break
% save CP2 CP
% save U2 U
% save V2 V
