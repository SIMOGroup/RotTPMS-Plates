function [fz, dfz, ddfz] = compute_shear_deformation_function(z,h,shear_func)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute shear deformation function at point z %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%% === Shear deformation function ===
switch shear_func
    case 1
        %E. Reissner, The effect of transverse shear deformation on the bending of elastic plates
        fz = 0; % 3DOF with no higher-order term
        dfz = 0;
        ddfz = 0;
    case 2
        %J. N. Reddy, A Simple Higher-Order Theory for Laminated Composite Plates,
        fz = z-4*z^3/(3*h^2); % Reddy
        dfz = 1-4*z^2/h^2;
        ddfz = -4*2*z/h^2;
    case 3
        %R. P. Shimpi, Refined Plate Theory and Its Variants
        fz = z-5*z^3/(3*h^2)+z/4; % Shimpi
        dfz = 1-5*z^2/h^2+1/4;
        ddfz = -5*2*z/h^2;
    case 4
        %H. Nguyen-Xuan, Isogeometric finite element analysis of composite sandwich plates using a higher order shear deformation theory
        fz = z-1/8*z-2/h^2*z^3+2/h^4*z^5;% H. Nguyen-Xuan
        dfz = 1-1/8-6/h^2*z^2+10/h^4*z^4;
        ddfz = -6/h^2*2*z+10/h^4*4*z^3;
    case 5
        % H. X. Nguyen, T. N. Nguyen, A refined quasi-3D isogeometric analysis for functionally graded microplates based on the modified couple stress theory
        fz = z-9*z+(10/h^2)*z^3+6/(5*h^4)*z^5+8/(7*h^6)*z^7;%Hoang X. Nguyen
        dfz = 1-9+(10/h^2)*3*z^2+6/(5*h^4)*5*z^4+8/(7*h^6)*7*z^6;
        ddfz = (10/h^2)*3*2*z+6/(5*h^4)*5*4*z^3+8/(7*h^6)*7*6*z^5;
    case 6
        %T. N. Nguyen, On the general framework of high order shear deformation theories for laminated composite plate structures: A novel unified approach
        fz = z - (17*z^3)/(10*h^2) + (22*z^5)/(25*h^4);%Tuan N. Nguyen
        dfz = 1 - (17*3*z^2)/(10*h^2) + (22*5*z^4)/(25*h^4);
        ddfz = (17*3*2*z)/(10*h^2) + (22*5*4*z^3)/(25*h^4);
    case 7
        %C. H. Thai, Generalized shear deformation theory for functionally graded isotropic and sandwich plates based on isogeometric approach
        fz = atan(sin(pi/h*z));%Chien H. Thai (Computers & Structures)
        dfz = pi/h*cos(pi/h*z)/(1+(sin(pi/h*z))^2);
        ddfz = - (pi^2*sin((pi*z)/h))/(h^2*(sin((pi*z)/h)^2 + 1)) - (2*pi^2*cos((pi*z)/h)^2*sin((pi*z)/h))/(h^2*(sin((z*pi)/h)^2 + 1)^2);
    otherwise
        disp('Error')
        pause
end

% --- Stretching effect deformation function ---
gz = 3/20 * dfz;
dgz = 3/20 * ddfz;
end
