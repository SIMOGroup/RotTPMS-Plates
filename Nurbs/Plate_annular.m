clear all
%clc
clf
format long

deg=3; ref=1;  %%% refinement

%plate circle 
p=2;
q=2;
U = [0 0 0 1/4 1/4 2/4 2/4 3/4 3/4 1 1 1]; %m=9
V = [0 0 0 1 1 1]; %n=2
%weights= [ 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1]';
R=2;
r=1;
CP(:,:,1)=[R (R+r)/2 r;R (R+r)/2 r;0 0 0;-R -(R+r)/2 -r;-R -(R+r)/2 -r;-R -(R+r)/2 -r;0 0 0;R (R+r)/2 r;R (R+r)/2 r];
CP(:,:,2)=[0 0 0;R (R+r)/2 r;R (R+r)/2 r;R (R+r)/2 r;0 0 0;-R -(R+r)/2 -r;-R -(R+r)/2 -r;-R -(R+r)/2 -r;0 0 0];
CP(:,:,3)=[1 1 1;1 1 1;1 1 1;1 1 1 ;1 1 1;1 1 1;1 1 1;1 1 1 ;1 1 1]*5;
%CP(:,:,4)=[ 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2);1 1];
CP(:,:,4)=[ 1 1 1; 1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1];


% U = [0 0 0 1/4 1/4 2/4 2/4 3/4 3/4 1 1 1]; %m=9
% V = [0 0 1 1]; %n=2
% %weights= [ 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1]';
% R=2;
% r=1;
% CP(:,:,1)=[R r;R r;0 0;-R -r;-R -r;-R -r;0 0;R r;R r];
% CP(:,:,2)=[0 0;R r;R r;R r;0 0;-R -r;-R -r;-R -r;0 0];
% CP(:,:,3)=[1 1; 1 1; 1 1; 1 1 ;1 1; 1 1; 1 1; 1 1 ;1 1]*5;
% %CP(:,:,4)=[ 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2);1 1];
% CP(:,:,4)=[ 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2);1 1];


% Gauss Points
ngauss = [deg+1 deg+1];

%=====================================================================
% REFINE
[CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);

R1 = refinement_vec_repeated_p3(U,ref);

R2 = refinement_vec_repeated_p3(V,ref);

[CP,U,V] = knot_refine_surf_repeated(p,q,U,V,CP,R1,R2);

plotNURBS_surf_El_CP(p,q,U,V,CP); hold on
% 
% save U3_1 U 
% save V3_1 V
% save CP3_1 CP

view(2)
