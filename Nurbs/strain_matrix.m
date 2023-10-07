function [A,B,D,E,F,H,As,Ds,Fs]=strain_matrix(mat,t,n)
% material property
E_m=mat(1,1);           E_c=mat(1,2);            % unit  N/m2
nu_m=mat(2,1);           nu_c=mat(2,2);              % 
% t : thickness
% n : order

% membrane matrix
A11=@(z)((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
A12=@(z)((nu_c-nu_m).*(1/2+z./t).^n+nu_m).*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
A66=@(z)((E_c-E_m).*(1/2+z./t).^n+E_m)./(2.*(1+((nu_c-nu_m).*(1/2+z./t).^n+nu_m)));
a11=quad(A11,-t/2,t/2);
a12=quad(A12,-t/2,t/2);
a66=quad(A66,-t/2,t/2);
A=[a11 a12 0;
   a12 a11 0;
   0   0   a66];

%bending matrix
B11=@(z)z.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
B12=@(z)z.*((nu_c-nu_m).*(1/2+z./t).^n+nu_m).*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
B66=@(z)z.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(2.*(1+((nu_c-nu_m).*(1/2+z./t).^n+nu_m)));
b11=quad(B11,-t/2,t/2);
b12=quad(B12,-t/2,t/2);
b66=quad(B66,-t/2,t/2);
B=[b11 b12 0;
   b12 b11 0;
   0   0   b66];

%curvature matrix
D11=@(z)z.^2.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
D12=@(z)z.^2.*((nu_c-nu_m).*(1/2+z./t).^n+nu_m).*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
D66=@(z)z.^2.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(2.*(1+((nu_c-nu_m).*(1/2+z./t).^n+nu_m)));
d11=quad(D11,-t/2,t/2);
d12=quad(D12,-t/2,t/2);
d66=quad(D66,-t/2,t/2);
D=[d11 d12 0;
   d12 d11 0;
   0   0   d66];

E11=@(z)z.^3.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
E12=@(z)z.^3.*((nu_c-nu_m).*(1/2+z./t).^n+nu_m).*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
E66=@(z)z.^3.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(2.*(1+((nu_c-nu_m).*(1/2+z./t).^n+nu_m)));
e11=quad(E11,-t/2,t/2);
e12=quad(E12,-t/2,t/2);
e66=quad(E66,-t/2,t/2);
E=[e11 e12 0;
   e12 e11 0;
   0   0  e66];

F11=@(z)z.^4.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
F12=@(z)z.^4.*((nu_c-nu_m).*(1/2+z./t).^n+nu_m).*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
F66=@(z)z.^4.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(2.*(1+((nu_c-nu_m).*(1/2+z./t).^n+nu_m)));
f11=quad(F11,-t/2,t/2);
f12=quad(F12,-t/2,t/2);
f66=quad(F66,-t/2,t/2);
F=[f11 f12 0;
   f12 f11 0;
   0   0  f66];

H11=@(z)z.^6.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
H12=@(z)z.^6.*((nu_c-nu_m).*(1/2+z./t).^n+nu_m).*((E_c-E_m).*(1/2+z./t).^n+E_m)./(1-((nu_c-nu_m).*(1/2+z./t).^n+nu_m).^2);
H66=@(z)z.^6.*((E_c-E_m).*(1/2+z./t).^n+E_m)./(2.*(1+((nu_c-nu_m).*(1/2+z./t).^n+nu_m)));
h11=quad(H11,-t/2,t/2);
h12=quad(H12,-t/2,t/2);
h66=quad(H66,-t/2,t/2);
H=[h11 h12 0;
   h12 h11 0;
   0   0  h66];

% shear
As=[a66  0; 0 a66];
Ds=[d66  0; 0 d66];
Fs=[f66  0; 0 f66];



