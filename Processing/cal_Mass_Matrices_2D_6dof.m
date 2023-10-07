function [M] = cal_Mass_Matrices_2D_6dof(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate mass matrices with Quasi-3D 6-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
uKnot = IGA.NURBS.uKnot; vKnot = IGA.NURBS.vKnot;
Inn = IGA.NURBS.Inn; Ien = IGA.NURBS.Ien; nel = IGA.NURBS.nel;
sdof = IGA.params.sdof; ndof = IGA.params.ndof; nGauss = IGA.params.nGauss; 

%% Used parameters from Plate
I = Plate.mat_mat.I;

%% ===== Initial stiffness matrices =====
M = sparse(sdof,sdof);

%% ===== Initial Gauss integration =====
[Pg, Wg] = gauleg(1,-1,nGauss);
nel_nza = 0; % Number of non-zero-area elements
tol = 1e-8;

%% ===== Mass matrices =====
for iel = 1:nel  % Loop over the elements
    % --- Element parameters ---
    sctr = Ien(iel,:);  % Control points indexes
    nn = length(sctr);  % Number of control points in the element
    for idof = 1:ndof  % Dofs of control points
        sctrM(idof:ndof:ndof*(nn-1) + idof) = ndof.*(sctr-1) + idof;
    end
    ni = Inn(Ien(iel,1),1);  % Index of the element in parametric domain
    nj = Inn(Ien(iel,1),2);

    % --- Gauss integration ---
    if abs(uKnot(ni)-uKnot(ni+1)) > tol && abs(vKnot(nj)-vKnot(nj+1)) > tol  % Check if the current element has nonzero area in the parametric domain
        nel_nza = nel_nza + 1;
        detJ2_xi = (uKnot(ni+1) - uKnot(ni))/2;  % Mapping parametric domain into natural domain of [[-1, 1]; [-1, 1]]
        detJ2_eta = (vKnot(nj+1) - vKnot(nj))/2;
  
        for iGauss = 1: nGauss  % Loop over the integration points
            for jGauss = 1: nGauss
                % Gauss points & weights
                gpt_xi = Pg(iGauss); gwt_xi = Wg(iGauss);
                gpt_eta = Pg(jGauss); gwt_eta = Wg(jGauss);
                
                % Map the point to parametric domain
                gpt_xi = (uKnot(ni+1)-uKnot(ni))/2*gpt_xi + (uKnot(ni+1)+uKnot(ni))/2; gwt_xi = gwt_xi*detJ2_xi;
                gpt_eta = (vKnot(nj+1)-vKnot(nj))/2*gpt_eta + (vKnot(nj+1)+vKnot(nj))/2; gwt_eta = gwt_eta*detJ2_eta;
                gwt = gwt_xi*gwt_eta;
                
                % Kinematic matrices of NURBS shape function
                [N, dNdxi, dNdxi2, dNdxy, dN2dxy, detJ1] = Kine_Shape_2nd_2D(IGA,ni,nj,gpt_xi,gpt_eta);
                
                % Mass matrix
                N0 = zeros(3,ndof*nn);  % (u_0, v_0, w_0)
                N0(1,1:ndof:ndof*nn) = N';
                N0(2,2:ndof:ndof*nn) = N';
                N0(3,3:ndof:ndof*nn) = N';
                
                N1 = zeros(3,ndof*nn);  % (-w_{0,x}, -w_{0,y}, 0)
                N1(1,3:ndof:ndof*nn) = -dNdxy(:,1)';
                N1(2,3:ndof:ndof*nn) = -dNdxy(:,2)';
                
                N2 = zeros(3,ndof*nn);  % (\beta_{x}, \beta_{y}, 0)
                N2(1,4:ndof:ndof*nn) = N';
                N2(1,5:ndof:ndof*nn) = N';
                
                N3 = zeros(3,ndof*nn);  % (0, 0, \beta_{z})
                N3(3,6:ndof:ndof*nn) = N';
                
                N_shape{1} = N0; N_shape{2} = N1; N_shape{3} = N2; N_shape{4} = N3;
                m = 0;
                for i = 1:4
                    for j = 1:4
                        m = m + N_shape{i}'*I{i,j}*N_shape{j}*gwt*detJ1;
%                         k = k + B{i}'*D{i,j}*B{j}*gwt*detJ1;
                    end
                end
                
                M(sctrM,sctrM) = M(sctrM,sctrM) + m;
            end
        end
    end
end
end
