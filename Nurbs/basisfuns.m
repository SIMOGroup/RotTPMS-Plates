%   This Function evaluates the basis functions at a given parameter value u.
% 
%  Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag: 
%     Berlin 1995; pp. 72-73.
% 
%  June 17, 2003
%  J. Austin Cottrell
%  CES Graduate Student
%  Texas Institute for Computational Engineering Science
%  University of Texas at Austin
%
%  Modified to codes Matlab by :
%  Chien Thai Hoang & Hung Nguyen Xuan
%
%   Faculty of Mathematics & Informatics, University of Natural Sciences
%   Vietnam   National University–HCM

function [shape]=basisfuns(i,pl,ml,u,u_knotl)

%     --------------variable declarations--------------------------------
% Input: i,p,m   %knot span, degree of curve, number of control points
%        u, u_knot  %parameter value, vector of knots
% Out:   ders(1,:)  % shape matrix
left=zeros(1,pl+1);
right=zeros(1,pl+1);
ndu=zeros(pl+1,pl+1);
%     -------------------------------------------------------------------

ndu(1,1) = 1;
for j = 1:pl;
left(j+1) = u - u_knotl(i+1-j);
right(j+1) = u_knotl(i+j) - u;
saved = 0;
   for r = 0:j-1;
       ndu(j+1,r+1) = right(r+2) + left(j-r+1);
       temp = ndu(r+1,j)./ndu(j+1,r+1);
       ndu(r+1,j+1) = saved + right(r+2).*temp;
       saved = left(j-r+1).*temp;
   end 
ndu(j+1,j+1) = saved;
end 

% load basis functions
for j = 0:pl;
shape(1,j+1) = ndu(j+1,pl+1);
end
shape;
return