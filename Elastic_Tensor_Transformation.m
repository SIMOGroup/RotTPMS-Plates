%% Initialization
clear
clc

% Voight notation 6x6 tensor
C = sym('C_%d%d', [6,6]);
assume(C(2,1) == C(1,2)); assume(C(3,1) == C(1,3)); assume(C(4,1) == C(1,4)); assume(C(5,1) == C(1,5)); assume(C(6,1) == C(1,6))
assume(C(3,2) == C(2,3)); assume(C(4,2) == C(2,4)); assume(C(5,2) == C(2,5)); assume(C(6,2) == C(2,6));
assume(C(4,3) == C(3,4)); assume(C(5,3) == C(3,5)); assume(C(6,3) == C(3,6)); 
assume(C(5,4) == C(4,5)); assume(C(6,4) == C(4,6)); 
assume(C(6,5) == C(5,6));

% Rotation matrix
syms cx sx cy sy cz sz real
assumeAlso(cx^2 + sx^2 == 1); assumeAlso(cy^2 + sy^2 == 1); assumeAlso(cz^2 + sz^2 == 1);
% disp(assumptions)

% Rx = [1, 0, 0; 0, cx, sx; 0, -sx, cx]; % Clockwise rotation
% Ry = [cy, 0, -sy; 0, 1, 0; sy, 0, cy]; % Clockwise rotation
% Rz = [cz, sz, 0; -sz, cz, 0; 0, 0, 1]; % Clockwise rotation

Rx = [1, 0, 0; 0, cx, -sx; 0, sx, cx]; % Counter-clockwise rotation
Ry = [cy, 0, sy; 0, 1, 0; -sy, 0, cy]; % Counter-clockwise rotation
Rz = [cz, -sz, 0; sz, cz, 0; 0, 0, 1]; % Counter-clockwise rotation

% Test
% alpha = deg2rad([45, 30, 60]); 
alpha = deg2rad([90, -45, 90]);  % [1; 1; 0]
% alpha = deg2rad([45, -35.264, 0]);  % [1; 1; 1]

%% Selection condition
% Tensor case
i_tensor = 3;
% [0] ~ Anisotropic:       C = [ C_11, C_12, C_13, C_14, C_15, C_16]
%                              [ C_21, C_22, C_23, C_24, C_25, C_26]
%                              [ C_31, C_32, C_33, C_34, C_35, C_36]
%                              [ C_41, C_42, C_43, C_44, C_45, C_46]
%                              [ C_51, C_52, C_53, C_54, C_55, C_56]
%                              [ C_61, C_62, C_63, C_64, C_65, C_66] %
% [1] ~ Monoclinic:        C = [ C_11, C_12, C_13,    0,    0, C_16]
%                              [ C_21, C_22, C_23,    0,    0, C_26]
%                              [ C_31, C_32, C_33,    0,    0, C_36]
%                              [    0,    0,    0, C_44, C_45,    0]
%                              [    0,    0,    0, C_54, C_55,    0]
%                              [ C_61, C_62, C_63,    0,    0, C_66] %
% [2] ~ Orthotropic:       C = [ C_11, C_12, C_13,    0,    0,    0]
%                              [ C_21, C_22, C_23,    0,    0,    0]
%                              [ C_31, C_32, C_33,    0,    0,    0]
%                              [    0,    0,    0, C_44,    0,    0]
%                              [    0,    0,    0,    0, C_55,    0]
%                              [    0,    0,    0,    0,    0, C_66] %
% [3] ~ Cubic-symmetric:   C = [ C_11, C_12, C_12,    0,    0,    0]
%                              [ C_12, C_11, C_12,    0,    0,    0]
%                              [ C_12, C_12, C_11,    0,    0,    0]
%                              [    0,    0,    0, C_44,    0,    0]
%                              [    0,    0,    0,    0, C_44,    0]
%                              [    0,    0,    0,    0,    0, C_44] %
% [4] ~ Isotropic:         C = [ C_11, C_12, C_12,             0,             0,             0]
%                              [ C_12, C_11, C_12,             0,             0,             0]
%                              [ C_12, C_11, C_11,             0,             0,             0]
%                              [    0,    0,    0, (C_11-C_12)/2,             0,             0]
%                              [    0,    0,    0,             0, (C_11-C_12)/2,             0]
%                              [    0,    0,    0,             0,             0, (C_11-C_12)/2] %

switch i_tensor
    case 1 
        C(1, [4,5]) = 0; C([4,5], 1) = 0;
        C(2, [4,5]) = 0; C([4,5], 2) = 0;
        C(3, [4,5]) = 0; C([4,5], 3) = 0;
        C(4, [6]) = 0; C([6], 4) = 0;
        C(5, [6]) = 0; C([6], 5) = 0;
    case 2
        C(1, [4,5,6]) = 0; C([4,5,6], 1) = 0;
        C(2, [4,5,6]) = 0; C([4,5,6], 2) = 0;
        C(3, [4,5,6]) = 0; C([4,5,6], 3) = 0;
        C(4, [5,6]) = 0; C([5,6], 4) = 0;
        C(5, [6]) = 0; C([6], 5) = 0;
    case 3
        C(1, [4,5,6]) = 0; C([4,5,6], 1) = 0;
        C(2, [4,5,6]) = 0; C([4,5,6], 2) = 0;
        C(3, [4,5,6]) = 0; C([4,5,6], 3) = 0;
        C(4, [5,6]) = 0; C([5,6], 4) = 0;
        C(5, [6]) = 0; C([6], 5) = 0;
        C(2, 2) = C(1, 1); C(3, 3) = C(1, 1);
        C(5, 5) = C(4, 4); C(6, 6) = C(4, 4); 
        C(1, 3) = C(1, 2); C(3, 1) = C(2, 1);
        C(2, 3) = C(1, 2); C(3, 2) = C(2, 1);
    case 4
        C(1, [4,5,6]) = 0; C([4,5,6], 1) = 0;
        C(2, [4,5,6]) = 0; C([4,5,6], 2) = 0;
        C(3, [4,5,6]) = 0; C([4,5,6], 3) = 0;
        C(4, [5,6]) = 0; C([5,6], 4) = 0;
        C(5, [6]) = 0; C([6], 5) = 0;
        C(2, 2) = C(1, 1); C(3, 3) = C(1, 1);
        C(4,4) = (C(1,1) - C(1,2)) / 2; C(5, 5) = C(4, 4); C(6, 6) = C(4, 4); 
        C(1, 3) = C(1, 2); C(3, 1) = C(2, 1);
        C(2, 3) = C(1, 2); C(3, 2) = C(2, 1);
end

% Rotation case
i_rotation = 0;
% [0] ~ Full rotation: R = Rx*Ry*Rz
% [1] ~ Rotation x: R = Rx
% [2] ~ Rotation y: R = Ry
% [3] ~ Rotation z: R = Rz

switch i_rotation
    case 0
%         R = Rx*Ry*Rz;        
        R = Rz*Ry*Rx;
    case 1
        R = Rx;
    case 2
        R = Ry;
    case 3
        R = Rz;
end

%% Method 1: Apply transformation on 3x3x3x3 tensor
C_tensor = generate(C);

ne = numel(C_tensor);                % number of tensor elements
nd = ndims(C_tensor);                % number of tensor dimensions, i.e. order of tensor
if (ne==3), nd = 1; end         % order of tensor is 1 in case of a 3x1 or 1x3 vector

tensor_new = C_tensor;                      % create output tensor
tensor_new(:) = 0;                     % fill output tensor with zeros; this way a symbolic tensor remains symbolic

iie = zeros(nd,1);              % initialise vector with indices of input tensor element
ioe = zeros(nd,1);              % initialise vector with indices of output tensor element
cne = cumprod(3*ones(nd,1))/3;  % vector with cumulative number of elements for each dimension (divided by three)

for oe = 1:ne                                 % loop over all output elements
   ioe = mod(floor((oe-1)./cne),3)+1;          % calculate indices of current output tensor element
   for ie = 1:ne                              % loop over all input elements
      pmx = 1;                                 % initialise product of transformation matrices
      iie = mod(floor((ie-1)./cne),3)+1;       % calculate indices of current input tensor element
      for id = 1:nd                           % loop over all dimensions
         pmx = pmx * R(ioe(id), iie(id));  % create product of transformation matrices
      end
      tensor_new(oe) = tensor_new(oe) + pmx * C_tensor(ie);       % add product of transformation matrices and input tensor element to output tensor element
   end
end
C_new_1 = ToMatrix(tensor_new);
% cs = floor([cos(alpha(1)), sin(alpha(1)), cos(alpha(2)), sin(alpha(2)), cos(alpha(3)), sin(alpha(3))]/1e-10)*1e-10;
% simplify(subs(C_new_1, {cx, sx, cy, sy, cz, sz}, {cs}))
clear C_tensor tensor_new

%% Method 2: Apply transformation on 6x6 tensor with paper: https://doi.org/10.1002/j.1538-7305.1943.tb01304.x
T = [     R(1,1)^2,      R(1,2)^2,      R(1,3)^2,             2*R(1,2)*R(1,3),             2*R(1,3)*R(1,1),             2*R(1,1)*R(1,2); ...
          R(2,1)^2,      R(2,2)^2,      R(2,3)^2,             2*R(2,2)*R(2,3),             2*R(2,3)*R(2,1),             2*R(2,1)*R(2,2); ...
          R(3,1)^2,      R(3,2)^2,      R(3,3)^2,             2*R(3,2)*R(3,3),             2*R(3,3)*R(3,1),             2*R(3,1)*R(3,2); ...
     R(2,1)*R(3,1), R(2,2)*R(3,2), R(2,3)*R(3,3), R(2,2)*R(3,3)+R(2,3)*R(3,2), R(2,1)*R(3,3)+R(2,3)*R(3,1), R(2,2)*R(3,1)+R(2,1)*R(3,2); ...
     R(3,1)*R(1,1), R(3,2)*R(1,2), R(3,3)*R(1,3), R(1,2)*R(3,3)+R(1,3)*R(3,2), R(1,3)*R(3,1)+R(1,1)*R(3,3), R(1,1)*R(3,2)+R(1,2)*R(3,1); ...
     R(1,1)*R(2,1), R(1,2)*R(2,2), R(1,3)*R(2,3), R(1,2)*R(2,3)+R(1,3)*R(2,2), R(1,3)*R(2,1)+R(1,1)*R(2,3), R(1,1)*R(2,2)+R(1,2)*R(2,1)];

% Tx = [1,      0,     0,         0,  0,   0;...
%       0,   cx^2,  sx^2,   2*sx*cx,  0,   0; ...
%       0,   sx^2,  cx^2,  -2*sx*cx,  0,   0; ...
%       0, -sx*cx, sx*cx, cx^2-sx^2,  0,   0; ...
%       0,      0,     0,         0, cx, -sx; ...
%       0,      0,     0,         0, sx,  cx];   
% Ty = [ cy^2, 0,   sy^2,   0,  -2*sy*cy, 0; ...
%           0, 1,      0,   0,         0, 0; ...
%        sy^2, 0,   cy^2,   0,   2*sy*cy, 0; ...
%           0, 0,      0,  cy,         0, sy; ...
%       sy*cy, 0, -sy*cy,   0, cy^2-sy^2, 0; ...
%           0, 0,      0, -sy,         0, cy];
% Tz = [ cz^2,  sz^2, 0,  0,   0,   2*cz*sz; ...
%        sz^2,  cz^2, 0,  0,   0,  -2*cz*sz; ...
%           0,     0, 1,  0,   0,         0; ...
%           0,     0, 0, cz, -sz,         0; ...
%           0,     0, 0, sz,  cz,         0; ...
%      -cz*sz, cz*sz, 0,  0,   0, cz^2-sz^2];
% T = Tx*Ty*Tz;
% T = Ty;
C_new_2 = T*C*T';
disp('--- Verify method 2 ---')
% Err_2 = simplify(C_new_1 - C_new_2);
% disp('Err_2 =')
% disp(Err_2);
simplify(C_new_2)
cs = floor([cos(alpha(1)), sin(alpha(1)), cos(alpha(2)), sin(alpha(2)), cos(alpha(3)), sin(alpha(3))]/1e-10)*1e-10;
simplify(subs(C_new_2, {cx, sx, cy, sy, cz, sz}, {cs}))

%% Method 3: Apply transformation on specific case only for Rz and 
% C = [ C_11, C_12, C_13,    0,    0,    0]
%     [ C_21, C_22, C_23,    0,    0,    0]
%     [ C_31, C_32, C_33,    0,    0,    0]
%     [    0,    0,    0, C_44,    0,    0]
%     [    0,    0,    0,    0, C_55,    0]
%     [    0,    0,    0,    0,    0, C_66]
% Note that: This method has an error in C36 of the result tensor
if i_tensor == 2 && i_rotation == 3
    c = cz; s = sz;
    H_1 = [    c^4,     2*c^2*s^2,     s^4,       4*c^2*s^2; ...
           c^2*s^2,     c^4 + s^4, c^2*s^2,      -4*c^2*s^2; ...
               s^4,     2*c^2*s^2,     c^4,       4*c^2*s^2; ...
            -c^3*s, c*s*(c^2-s^2),   c*s^3, 2*c*s*(c^2-s^2); ...
            -c*s^3, c*s*(s^2-c^2),   c^3*s, 2*c*s*(s^2-c^2); ...
           c^2*s^2,    -2*c^2*s^2, c^2*s^2,    (c^2-s^2)^2];
    H_2 = [c^2, s^2; c*s, -c*s; s^2, c^2];
    temp_1 = H_1*[C(1,1); C(1,2); C(2,2); C(6,6)];
    temp_2 = H_2*[C(4,4); C(5,5)];
    temp_3 = H_2*[C(1,3); C(2,3)];

    C_new_3 = [temp_1(1), temp_1(2), temp_3(1),         0,         0, temp_1(4); ...
               temp_1(2), temp_1(3), temp_3(3),         0,         0, temp_1(5); ...
               temp_3(1), temp_3(3),    C(3,3),         0,         0, temp_3(2); ...
                       0,         0,         0, temp_2(1), temp_2(2),         0; ...
                       0,         0,         0, temp_2(2), temp_2(3),         0; ...
               temp_1(4), temp_1(5), temp_3(2),         0,         0, temp_1(6)];
    disp('--- Verify method 3 ---')
    Err_3 = simplify(C_new_1 - C_new_3);
    disp('Err_3 =')
    disp(Err_3);

    % Fix the error by fixing C36
    H_2_mod = [c^2, s^2; -c*s, c*s; s^2, c^2];
    temp_3_mod = H_2_mod*[C(1,3); C(2,3)];

    C_new_3_mod = [    temp_1(1),     temp_1(2), temp_3_mod(1),         0,         0,     temp_1(4); ...
                       temp_1(2),     temp_1(3), temp_3_mod(3),         0,         0,     temp_1(5); ...
                   temp_3_mod(1), temp_3_mod(3),        C(3,3),         0,         0, temp_3_mod(2); ...
                               0,             0,             0, temp_2(1), temp_2(2),             0; ...
                               0,             0,             0, temp_2(2), temp_2(3),             0; ...
                       temp_1(4),     temp_1(5), temp_3_mod(2),         0,         0,     temp_1(6)];
    disp('--- Verify method 3 modified ---')
    Err_3_mod = simplify(C_new_1 - C_new_3_mod);
    disp('Err_3_mod =')
    disp(Err_3_mod);
end
% subs(C_new_3, {cx, sx, cy, sy, cz, sz}, {cos(alpha(1)), sin(alpha(1)), cos(alpha(2)), sin(alpha(2)), cos(alpha(3)), sin(alpha(3))})

%% Functions
function C = generate(CH)
C = sym('ts_%d%d', [3,3,3,3]);
% C = zeros(3,3,3,3);
for i = 1:6
    for j = 1:6
        [a,b] = change(i);
        [c,d] = change(j);
        C(a,b,c,d) = CH(i,j);
    end
end
for i = 1:3
    if (i == 3)
        j = 1;
    else
        j = i+1;
    end
    for m = 1:3
        if (m == 3)
            n = 1;
        else
            n = m+1;
        end
        C(j,i,n,m) = C(i,j,m,n);
        C(j,i,m,n) = C(i,j,m,n);
        C(i,j,n,m) = C(i,j,m,n);
        C(j,i,m,m) = C(i,j,m,m);
        C(m,m,j,i) = C(m,m,i,j);
    end
end
end

% change the index 4 5 6 to 23 31 12
function [a,b] = change(w)

if (w < 4)
    a = w;
    b = w;
else
    if (w == 4)
        a = 2;
        b = 3;
    else
        if (w == 5)
            a = 3;
            b = 1;
        else
            if (w==6)
                a = 1;
                b = 2;
            end
        end
    end
end

end

function CH = ToMatrix(C)
% CH = zeros(6,6);
CH = sym('tsn_%d%d', [6,6]);
for i = 1:6
    for j = 1:6
        [a,b] = change(i);
        [c,d] = change(j);
        CH(i,j) = C(a,b,c,d);     
    end
end
end