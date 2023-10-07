function [bcdof, bcval] = cal_bcdof_2D_6dof(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate dofs of boundary conditions with 5-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
mcp = IGA.NURBS.mcp; ncp = IGA.NURBS.ncp;
ndof = IGA.params.ndof;

%% Used parameters from Plate
bc_case = Plate.bc.bc_case;

%% ===== Initial boundary condition arrays =====
bcdof = []; bcval = [];

%% ===== Control points on boundaries =====
left1 = [1: mcp: mcp*(ncp-1)+1];           % on perimeter
left2 = [2: mcp: mcp*(ncp-1)+2];           % 01 control point inside
right1 = [mcp: mcp: mcp*ncp];              % on perimeter
right2 = [mcp-1: mcp: mcp*ncp-1];          % 01 control point inside
upper1 = [mcp*(ncp-1)+1: mcp*ncp];         % on perimeter
upper2 = [mcp*(ncp-2)+1: mcp*(ncp-1)];     % 01 control point inside
lower1 = [1: mcp];                         % on perimeter
lower2 = [mcp+1: 2*mcp];                   % 01 control point inside

%% ===== Dofs of control points on boundaries (dofs = u_{0}, v_{0}, w_{0}, \beta_{x}, \beta_{y}, \beta_{z}) =====
switch bc_case
    case 1 % Fully simply supported (SSSS)
        for idof = [3, 6]  % Constraint w_{0}, \beta_{z}) at perimeter
            bcdof = [bcdof ndof*(left1-1)+idof];    % Left
            bcdof = [bcdof ndof*(right1-1)+idof];   % Right
            bcdof = [bcdof ndof*(upper1-1)+idof];   % Upper
            bcdof = [bcdof ndof*(lower1-1)+idof];   % Lower
        end
        for idof = [1, 4]  % Constraint u_{0}, \beta_{x} at perimeter on Upper and Lower
            bcdof = [bcdof ndof*(upper1-1)+idof];   % Upper
            bcdof = [bcdof ndof*(lower1-1)+idof];   % Lower
        end
        for idof = [2, 5]  % Constraint v_{0}, \beta_{y} at perimeter on Left and Right
            bcdof = [bcdof ndof*(left1-1)+idof];    % Left
            bcdof = [bcdof ndof*(right1-1)+idof];   % Right
        end

%         % If using immovable edge SSSS2
%         for idof = [1, 2, 4, 5]  % Constraint u_{0}, v_{0}, \beta_{x}, \beta_{y} at perimeter on Left and Right
%             bcdof = [bcdof ndof*(left1-1)+idof];    % Left
%             bcdof = [bcdof ndof*(right1-1)+idof];   % Right
%             bcdof = [bcdof ndof*(upper1-1)+idof];   % Upper
%             bcdof = [bcdof ndof*(lower1-1)+idof];   % Lower
%         end
    case 2 % Fully clamped (CCCC)
        for idof = [1, 2, 3, 4, 5, 6]  % Constraint u_{0}, v_{0}, w_{0}, \beta_{x}, \beta_{y} at perimeter
            bcdof = [bcdof ndof*(left1-1)+idof];    % Left
            bcdof = [bcdof ndof*(right1-1)+idof];   % Right
            bcdof = [bcdof ndof*(upper1-1)+idof];   % Upper
            bcdof = [bcdof ndof*(lower1-1)+idof];   % Lower
        end
        
        for idof = [3, 6]  % Constraint dw_{0}dx, dw_{0}dy, d\beta_{z}dx, d\beta_{z}dy at perimeter (with a fine mesh)
            bcdof = [bcdof ndof*(left2-1)+idof];    % Left
            bcdof = [bcdof ndof*(right2-1)+idof];   % Right
            bcdof = [bcdof ndof*(upper2-1)+idof];   % Upper
            bcdof = [bcdof ndof*(lower2-1)+idof];   % Lower
        end
end
bcval = [bcval, zeros(1,length(bcdof))];
end
