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
clear all; close all
addpath Nurbs
addpath Processing
addpath Boundary
global  nsd nshl nnode nel 
global u_knot v_knot p q mcp ncp 
global u_knot1 v_knot1 p1 q1 mcp1 ncp1 
global gcoord order
format long

% ===================== MATERIAL PROPERTY =================================
% % geometry
R= 0.5;     % length of the plate
L=1; D=1;
h= R/100;
% h= input('thickness of plate = ');
c_cof=-4/(3*h^2);
f0=1;  % distributed force

% % isotropic
e2=1e7;e1=e2;
miu12=0.3;miu21=miu12*e2/e1;
g23=0.5*e2/(1+miu12);g13=g23;g12=g13;
an1=1e-6;   an2=an1;                      % thermal expansion cofficent

% %  orthotropic                          ?/????? property of materials????
% e2=1;     e1=15*e2;
% miu12=0.3; miu21=miu12;%*e2/e1;
% g23=0.3356*e2;    g13=0.5*e2;    g12=g13;
% an2=1e-6;   an1=0.015*an2;                  % thermal expansion cofficent

% % % material
% e2=1;e1=40*e2;
% g23=0.5*e2;g13=0.6*e2;g12=g13;
% miu12=0.25;miu21=miu12*e2/e1;
% an1=1e-6;   an2=2*an1;                  % thermal expansion cofficent
% %!!!!!!!!!!!!!!!!!!!!!!!!1 check

% input number layer of laminae
nplayer=1;  %[0/90/90/0]
E_module(1,1:nplayer)=e1;
E_module(2,1:nplayer)=e2;
nuy(1,1:nplayer)=miu12;
nuy(2,1:nplayer)=miu21;
alpha(1,1:nplayer)=an1;
alpha(2,1:nplayer)=an2;
G(1,1:nplayer)=g12;
G(2,1:nplayer)=g23;
G(3,1:nplayer)=g13;

for i=1:nplayer+1
    z(i)=-h/2+(i-1)*h/nplayer;      
end

% input orientation of fiber
ale=0;
theta=[0]*pi/180;
A1=0;B=0;D1=0;E=0;F1=0;H=0;A2=0;D2=0;F2=0;
for k=1:nplayer
    [Q,Qs]=CALCQ(k,theta,E_module,nuy,G);
    [Q_Bar]=CALCQBAR(k,theta,Q);
    [Q_Bars]=CALCQsBAR(k,theta,Qs);
    A1=A1+Q_Bar*(z(k+1)-z(k));
    B=B+0.5*Q_Bar*(z(k+1)^2-z(k)^2);
    D1=D1+(1/3)*Q_Bar*(z(k+1)^3-z(k)^3);
    E=E+(1/4)*Q_Bar*(z(k+1)^4-z(k)^4);
    F1=F1+(1/5)*Q_Bar*(z(k+1)^5-z(k)^5);
    H=H+(1/7)*Q_Bar*(z(k+1)^7-z(k)^7);
    A2=A2+Q_Bars*(z(k+1)-z(k));
    D2=D2+(1/3)*Q_Bars*(z(k+1)^3-z(k)^3);
    F2=F2+(1/5)*Q_Bars*(z(k+1)^5-z(k)^5);
end

% ========================= DATA OF NURBS ================================
% Read in mesh information from file mesh_knot.dat

% p = input('degree of NURBS basis function = ');
% q = p;

p=2;                      % degree of curve in u
q=p;                      % degree of curve in v
%Input num element along x, y direction
% nx = input('Number element in x-dir = ');
% ny = input('Number element in y-dir = ');
nx=3;ny=nx;
order=p-1;                      % C^k continuity element
[mcp,ncp,u_knot,v_knot] = repeated_knot(p,q,nx,ny,order);

p1=p%-1;
q1=q%-1;
[mcp1,ncp1,u_knot1,v_knot1] = repeated_knot(p1,q1,nx,ny,order);

% ========================== GEOMETRY AND MESHING =========================
numx     = mcp-1;% the number of elements in the x-direction (beam length)
numy     = ncp-1;% the number of elements in the y-direction (width ength)

%  input data for nodal coordinate values

CP=importdata('P2_3.mat');
b_net(:,:,1)=CP(:,:,1)';
b_net(:,:,2)=CP(:,:,2)';
b_net(:,:,3)=CP(:,:,4)';

gcoord(:,1)=reshape(b_net(:,:,1),mcp*ncp,1);
gcoord(:,2)=reshape(b_net(:,:,2),mcp*ncp,1);
gcoord(:,3)=reshape(b_net(:,:,3),mcp*ncp,1);

% icount=0;
% for j = 1:ncp
%       for i = 1:mcp
%           icount = icount + 1;
%             b_net(i,j,:) = gcoord(icount,:);
%       end
% end  
%----------------------------------------------
% gcoord=meshRectangularCoord(L,W,numx,numy);
nodes=connecElementQ4(numx,numy);
%
figure('Name','mesh of deflection')

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
 axis square
 axis([-0.6 0.6 -0.6 0.6])
 %}

% connectivity for pressure
gcoordP=gcoord;
b_net1=b_net;
%{
figure('Name','mesh of displacement')
nodesP=connecElementQ4(mcp1-1,ncp1-1);   % CP for p-1 order
 for iel=1:size(nodesP,1)
    nd=nodesP(iel,:);
    x=gcoordP(nd,1);
    y=gcoordP(nd,2);    
%    xx=mean(x); 
%    yy=mean(y); 
     patch(x,y,'w')
     text(x(1),y(1),num2str(nd(1)));
     text(x(2),y(2),num2str(nd(2)));
     text(x(3),y(3),num2str(nd(3)));
     text(x(4),y(4),num2str(nd(4)));
 end
% break
 %}

%%% --------------------------------------------------------------------%%%
%number of elements, degrees of freedom
nel=nx*ny;                  % number of element
nsd=2;                      % number of spatial dimension
ngauss=max(p,q)+1;
nnode=mcp*ncp;              % number of control point for displacement
nshl = (p+1)*(q+1);         % number of local shape functions
ndofU=1;                    % w
nunk = ndofU*nnode;

ndofP=4;                    % u, v, theta_x, theta_y
nunkP = ndofP*mcp1*ncp1;    % global dof of pressure 

sdof  = nunk+nunkP;         % global dof

nnelU = (p+1)*(q+1);        % number of node per element for displacement
nnelP = (p1+1)*(q1+1);      % number of node per element for pressure

%%% =====================nurbs connectivity===================
% generate connectivities     
[ien,inn]=genIEN_INN_2D_repeated_Ck(p,q,mcp,ncp);
% generate connectivities     
[ienP,innP]=genIEN_INN_2D_repeated_Ck(p1,q1,mcp1,ncp1);
 
% % ========================= THERMAL LOAD ==================================
% delta_T=1;
% N_th=0; M_th=0; P_th=0;
% for k=1:nplayer
%     [Q,Qs]=CALCQ(k,theta,E_module,nuy,G);
%     [Q_Bar]=CALCQBAR(k,theta,Q);
%     glo_al=CALC_ALPHA(k, theta, alpha);         % thermal expansion 3x1 
%     N_th=N_th+Q_Bar*glo_al*delta_T*(z(k+1)-z(k));           % 3x1
%     M_th=M_th+1/2*Q_Bar*glo_al*delta_T*(z(k+1)^2-z(k)^2);   % 3x1
%     P_th=P_th+1/4*Q_Bar*glo_al*delta_T*(z(k+1)^4-z(k)^4);   % 3x1
% end

% ========================= BUILD STIFFNESS MATRIX ========================
K =sparse(sdof,sdof);     % stiffness matrix
% Kg=zeros(sdof,sdof);    % mass matrix
F =zeros(sdof,1);

% % Pre-bucklling stress
%  prestress=[N_th(1), N_th(3);                
%             N_th(3), N_th(2)];

% stiffness matrix
% get gaussian points and weights;
[gp,gw]=genGP_GW(ngauss);
tol=1e-8;% in order to check u_knot(ni) matches u_knot(ni+1) or not ...
% loop over elements;
for iel = 1:nel  
  edofU=ien(iel,:);           % element scatter vector for w
  nn=length(edofU);
  sctr=ienP(iel,:);           % element scatter vector for u,v, tx,ty
  np=length(sctr);
  edofP(1:4:4*np) = 4.*sctr-3 ;
  edofP(2:4:4*np) = 4.*sctr-2 ;
  edofP(3:4:4*np) = 4.*sctr-1 ;
  edofP(4:4:4*np) = 4.*sctr ;
  edofP=edofP+nunk;
  sctrB=[edofU, edofP];
%  check to see ifmlv current element has nonzero area;
  ni = inn(ien(iel,1),1);% get NURBS coordinates
  nj = inn(ien(iel,1),2);
% element has positive area in the parametric domain
  if(abs(u_knot(ni)-u_knot(ni+1))>tol)&&(abs(v_knot(nj)-v_knot(nj+1))>tol)
     da =(u_knot(ni+1) - u_knot(ni))*(v_knot(nj+1) - v_knot(nj))/4;
    
     %  set up element parameters;
     for igauss = 1: ngauss
        for jgauss = 1: ngauss
            gwt = gw(igauss)*gw(jgauss)*da; 
            [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd_repeated_Ck(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net,p,q,nnelU,mcp,ncp);
            [NP,dNdxiP,dNdxyP,dN2dxyP,detjP]=Kine_SHAPE_2nd_repeated_Ck(iel,gp(igauss),gp(jgauss),u_knot1,v_knot1,b_net1,p1,q1,nnelP,mcp1,ncp1);
            
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [B] strain matrix of membrane and bending
        B_m=B_menbrance(dNdxyP,nnelU,nnelP);
        B_b0=B_bending(dNdxyP,nnelU,nnelP);
        B_b1=B_bending_1(dNdxyP,dN2dxy,nnelU,nnelP,c_cof);
        
 % [B] strain matrix of shear
        B_s0=B_shear(NP,dNdxy,nnelU,nnelP);
        B_s1=3*c_cof*B_s0;    
        
  % [B] strain matrix of geometry
        B_g=B_geometry(dNdxy,nnelU,nnelP);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        k=(B_m'*A1*B_m  + B_b0'*B*B_m  + B_m'*B*B_b0 + B_b0'*D1*B_b0 +...      %)*gwt*detjP + ...
           B_b1'*E*B_m  + B_b1'*F1*B_b0+ B_m'*E*B_b1 + B_b0'*F1*B_b1 +...
           B_b1'*H*B_b1 + B_s0'*A2*B_s0+ B_s1'*D2*B_s0 + ...
           B_s0'*D2*B_s1 +B_s1'*F2*B_s1)*gwt*detj;
     
% Assemply stiffness matrix       
         K(sctrB,sctrB)=K(sctrB,sctrB)+k; 

         end
     end % end integration loop
  end %end if
end % end loop over elements

% % thermal load
% [F]=thermal_load(ngauss,nel,inn,ien,b_net,N_th, M_th, P_th,l,F);
[F,nel_nza]=ForceNurbs_uniform(ngauss,nel,inn,ien,b_net,f0,F);

% boundary conditions
%bc for disp
bcdof=[];   bcval=[];
[bcdof, bcval]=bcdof_CCCC(4,mcp1,ncp1);
bcdof=bcdof+nunk;

% bc for deflection
lowerEdge=[1:mcp];
upperEdge=[mcp*(ncp-1)+1:ncp*mcp];
leftEdge=[mcp+1:mcp:mcp*(ncp-2)+1];
rightEdge=[2*mcp:mcp: mcp*(ncp-1)];
bcdof=[bcdof lowerEdge, upperEdge, leftEdge, rightEdge];
bcval=[bcval zeros(1,length([lowerEdge, upperEdge, leftEdge, rightEdge]))];

%
% apply constrain
[Stiff,force]=feaplyc2(K,F,bcdof,bcval);

% ===========================STATIC ANALYSIS ==============================
% % displacement
U=Stiff\force;
% strain energy
EU=1/2*U'*K*U

Dm=(e2*h^3)/(12*(1-miu12^2));
M=mcp*(ncp-1)/2+(mcp+1)/2;

W_center=U(M);
Dis_Nurbs=W_center/(f0*R^4/(64*Dm))

%%% central moment
coord = zeros(1,2);
ksi   = 0;
eta   = 0;
iter  = 10;
c_ele=(nel+1)/2;   % for p=2 or p=4
% c_ele=nx*((ncp+1)/2-2)+(mcp+1)/2-2;   % for p=3
sctr = ien(c_ele,:);
nodes = gcoord(sctr,:);
% xm=mean(xx);
% ym=mean(yy);
inc = 1;
while (inc < iter)
   [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd_repeated_Ck(c_ele,coord(1,1),coord(1,2),u_knot,v_knot,b_net,p,q,nnelU,mcp,ncp);
    x = N'*nodes(:,1);
    y = N'*nodes(:,2);
    dx_dxi = dNdxi(:,1)' * nodes(:,1);
    dy_deta = dNdxi(:,2)' * nodes(:,1);
    dy_dxi = dNdxi(:,1)' * nodes(:,2);
    dx_deta = dNdxi(:,2)' * nodes(:,2);
 
    deltaX = x - gcoord(M,1);
    deltaY = y - gcoord(M,2);
    delta=[deltaX;deltaY];
    F=[dx_dxi dy_deta; dy_dxi dx_deta];
    invF=inv(F);
    
    ksi = ksi - invF(1,:)*delta;
    eta = eta - invF(2,:)*delta;
    
    coord(1) = ksi;
    coord(2) = eta;
    inc = inc + 1;
end

  np=length(sctr);
  edofP(1:4:4*np) = 4.*sctr-3 ;
  edofP(2:4:4*np) = 4.*sctr-2 ;
  edofP(3:4:4*np) = 4.*sctr-1 ;
  edofP(4:4:4*np) = 4.*sctr ;
  edofP=edofP+nunk;
  sctrB=[sctr, edofP];
  
  [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd_repeated_Ck(c_ele,coord(1,1),coord(1,2),u_knot,v_knot,b_net,p,q,nnelU,mcp,ncp);
   B_m=B_menbrance(dNdxy,nnelU,nnelP);
   B_b0=B_bending(dNdxy,nnelU,nnelP);
   B_b1=B_bending_1(dNdxy,dN2dxy,nnelU,nnelP,c_cof); 
  
  Moment=B*B_m*U(sctrB)+D1*B_b0*U(sctrB)+F1*B_b1*U(sctrB);
  M_x=Moment(1)/((1+miu12)*f0*R^2/16)
  
  cen_def=N'*U(ien(c_ele))/(f0*R^4/(64*Dm))
