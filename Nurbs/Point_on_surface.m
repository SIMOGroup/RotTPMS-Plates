% Subroutine eval_SHAPE.f consumes an element number and the coordinates in the
% parent element of an integration point and returns the vector of all local
% basis functions evaluated at the point and the matrix of gradients for all
% nonzero bais functions with respect to parameters u and v and with respect
% to x and y.
%
%  June 17, 2003
%  J. Austin Cottrell
%  CES Graduate Student
%  Texas Institute for Computational Engineering Science
%  University of Texas at Austin
%
%  Modify to codes Matlab by :
%  Hung Nguyen Xuan
%
%   Faculty of Mathematics & Informatics, University of Natural Sciences
%   Vietnam   National University HCM

function [S]=Point_on_surface(e,u_hat,v_hat,u_knot,v_knot,b_net)

global nsd nshl p q mcp ncp

shl=zeros(nshl,1);
denom_sum=0;
% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;
[ien,inn]=genIEN_INN_2D(p,q,mcp,ncp);
ni = inn(ien(e,1),1);
nj = inn(ien(e,1),2);
%     get u and v coordinates of integration point;
u =((u_knot(ni+1)-u_knot(ni))*u_hat +u_knot(ni+1) + u_knot(ni))/2;
v =((v_knot(nj+1)-v_knot(nj))*v_hat +v_knot(nj+1) + v_knot(nj))/2;

% evaluate 1d size functions and derivatives each direction
% size and M and N = number of derviatives+shape functions, degree of
% poynomial.
% row 1 of M and N => shape functions
% i^{th} row (i > 1) => i^{th} derivative of the shape function
% calculate in u direction
M = dersbasisfuns(ni,p,mcp,u,u_knot) ;

% calculate in v direction
N = dersbasisfuns(nj,q,ncp,v,v_knot) ;

% form basis functions and derivatives dr./du and dr./dv;
icount = 0;

for j = 0:q
    for i = 0:p
        icount = icount+1;
        
        % basis functions
        shl(icount,1) = M(1,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1);
        denom_sum = denom_sum + shl(icount);
    end
end

% basis functions
shl = shl/denom_sum;

S=zeros(nsd,1);
icount = 0;
for j = 0: q
    for i = 0: p
        icount = icount + 1;
        S(1,1) = S(1,1) + b_net(ni-i,nj-j,1)*shl(icount,1);
        S(2,1) = S(2,1) + b_net(ni-i,nj-j,2)*shl(icount,1);

    end
end



end