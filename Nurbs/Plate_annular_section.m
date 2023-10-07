clear all
%clc
clf
format long

deg=2; ref=8;  %%% refinement

%%%%%%%% section annular plate

theta=90*pi/180;
L=4;
R=1;
p=2;
q=2;
U=[0 0 0 0.5 1 1 1];%m=3
V=[0 0 0 1 1 1];%n=2
%Control Point coordinates
CP(:,:,1)=-[-R -L/2 -L; -R -L/2 -L;1-sqrt(2) -0.75 -L;0 0 0];
CP(:,:,2)=[0 0 0;sqrt(2)-1 0.75 L;R L/2 L;R L/2 L ];
CP(:,:,3)=[1 1 1;1 1 1 ;1 1 1; 1 1 1];
CP(:,:,4)=[1 1 1; 1/sqrt(2) 1 1;1/sqrt(2) 1 1; 1 1 1];

% %%%%%%%% section annular plate
% 
% theta=90*pi/180;
% R=2;
% r=1;
% p=2;
% q=1;
% U=[0 0 0 1 1 1];
% V=[0 0 1 1];
% %Control Point coordinates
% CP(:,:,1)=[R r; R r; R*cos(theta) r*cos(theta)];
% CP(:,:,2)=[0 0; R*tan(theta/2) r*tan(theta/2); R*sin(theta) r*sin(theta)];
% CP(:,:,3)=[0 R; 0 R ; 0 R];
% CP(:,:,4)=[1 1;  cos(theta/2)  cos(theta/2); 1 1];

%=====================================================================
% REFINE
[CP,U,V,p,q] = degree_elevate_surf(p,q,U,V,CP,deg-p,deg-q);
R1 = Refinement_vec(U,ref);
R2 = Refinement_vec(V,ref);

[CP,u_knot,v_knot] = Knot_Refine_Surf(p,q,U,V,CP,R1,R2);

plotNURBS_surf_El_CP(p,q,u_knot,v_knot,CP); hold on
%view(2)

% section computational IgA

p=deg;
q=p;
ngauss=p+1;                 % number of gauss point in integration

% input control net
b_net(:,:,1)=CP(:,:,1);
b_net(:,:,2)=CP(:,:,2);
b_net(:,:,3)=CP(:,:,4);

mcp = length(b_net(:,1,1));
ncp = length(b_net(1,:,1));

clear CP U V R1 R2
gcoord(:,1)=reshape(b_net(:,:,1),mcp*ncp,1);
gcoord(:,2)=reshape(b_net(:,:,2),mcp*ncp,1);
gcoord(:,3)=reshape(b_net(:,:,3),mcp*ncp,1);


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

R=2;    %radius of the circular
h= R*0.1; %thickness
f0=1;
l=-4/(3*h^2);
n=0; %volume fraction exponent range

% ======================= material ===============================
%{
%  isotropic
 material=[ 1      1;    % Young modulus
           0.3      0.3;    % poisson
           1e-6    1e-6;    % thermal expansion
           204     10.4;    % thermal conductivity
            1        1];    % density
%}
%{
 % Al(metal)/Al2O3(ceramic)
 material=[70e9     380e9;    % Young modulus
           0.3      0.3;    % poisson
           23e-6    7.4e-6;    % thermal expansion
           204      10.4;    % thermal conductivity
           2707     3800];    % density
%}       
%{
% Al/ZrO2-1
 material=[70e9     200e9;    % Young modulus
           0.3      0.3;    % poisson
           23e-6    10e-6;    % thermal expansion
           204      2.09;    % thermal conductivity
           2707     5700];    % density       
 %}
 %{
 % Al/ZrO2-2
 material=[70e9     151e9;    % Young modulus
           0.3      0.3;    % poisson
           23e-6    10e-6;    % thermal expansion
           204      2.09;    % thermal conductivity
           2707     3000];    % density  
%}       
%{
 % Ti-Al-4V/Al2O3
 material=[107.5e9    320.2e9;    % Young modulus
           0.298      0.26;    % poisson
           6.9e-6     7.2e-6;    % thermal expansion
           18.1       10.4;    % thermal conductivity
           4429       3750];    % density       
%}
 
 %
 % Al(metal)/Silicon carbide
 
  material=[70e9  70e9   ;    % Young modulus
            0.3      0.3;       % poisson
            23e-6    10e-6;     % thermal expansion
            204      2.09;      % thermal conductivity
            2707     5700];     % density       
     
%}


 % khai bao lai vat lieu
 E_m=material(1,1);     E_c=material(1,2);
 nu_m=material(2,1);    nu_c=material(2,2);
 an_m=material(3,1);    an_c=material(3,2);
 k_m=material(4,1);     k_c=material(4,2);
 rho_m=material(5,1);   rho_c=material(5,2);
% ====== caculate membrane,bending, curvature, shear, mass matrix========
[A1,B,D1,E,F1,H,A2,D2,F2]=strain_matrix(material,h,n);


% *************************************************************************
%                     P R E - P R O C E S S I N G                         %                      
% *************************************************************************
%number of elements, degrees of freedom
nnode=mcp*ncp;         % number of control point
nshl = (p+1)*(q+1);    % number of local shape functions
nel=(mcp-p)*(ncp-q); % number of element
nsd=2;                       % number of spatial dimension
ndof=5;
sdof=nnode*ndof;
nnel=nshl;
%%%%% boundary conditions
%[bcdof, bcval]=bcdof_CCCC_composite_5dof(ndof,mcp-1,mcp-1);
[bcdof, bcval]=bcdof_SSSS_composite_5dof(ndof,mcp-1,mcp-1);


%----------------------------------------------
%  initialization of matrices and vectors
%----------------------------------------------
K=sparse(sdof,sdof); % stiffness matrix
F=zeros(sdof,1); % mass matrix
%Read in mesh information from file mesh.dat to load b_net
% icount=0;
% for j = 1:ncp
%       for i = 1:mcp
%           icount = icount + 1;
%             b_net(i,j,:) = gcoord(icount,:);
%       end
% end  
% generate connectivities     
[ien,inn]=genIEN_INN_2D(p,q,mcp,ncp);
  
%*************************************************************************
%Loop over interior elements to build lhsk and force vector contributions 
%of elements (body loads).
%*************************************************************************
[K]=KmatNurbs_plate_composite_5dof_C1(ngauss,nel,inn,ien,b_net,A1,B,D1,E,F1,H,A2,D2,F2,l,K);
[F]=ForceNurbs_plate_composite_uniform_5dof(ngauss,nel,inn,ien,b_net,f0,F);
%-----------------------------
%   apply boundary conditions
%-----------------------------
 
[KK,FF]=feaplyc2(K,F,bcdof,bcval);
% SOLVE SYSTEM
[LL UU]=lu(KK);
utemp=LL\FF;
disp=UU\utemp;
dis=disp(3:5:sdof-2);

M=(numx+1)*numy/2+(numx+2)/2;
% dis=disp(5*M-2,1);

Dc=(E_m*h^3)/(12*(1-nu_m^2));
w_exact=Dc*dis(M)/(f0*(R-r)^4)


