function [dR,ddR] = deriv_Nurbsbasisfunc2(p,i,u,U,q,j,v,V,CP)
 

if (i==0); i = findspan(u,U,length(CP(:,1,1)));  end
if (j==0); j = findspan(v,V,length(CP(1,:,1)));  end

ne = (p+1)*(q+1);       
N = deriv2(i,p,u,U);    
M = deriv2(j,q,v,V);    

R = zeros(ne,1);
dR = zeros(ne,2);
ddR = zeros(ne,3);    
sum = 0;
dsum = zeros(2,1);
ddsum = zeros(3,1);

k = 0;
for c = 0:q
  for b = 0:p
    k = k+1;
    R(k) = N(1,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    sum = sum + R(k);
            
    dR(k,1) = N(2,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    dsum(1)  = dsum(1) + dR(k,1);
    dR(k,2) = N(1,b+1)*M(2,c+1)*CP(i-p+b,j-q+c,4);
    dsum(2)  = dsum(2) + dR(k,2);
            
    ddR(k,1) = N(3,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    ddsum(1)  = ddsum(1) + ddR(k,1);
    ddR(k,2) = N(1,b+1)*M(3,c+1)*CP(i-p+b,j-q+c,4);
    ddsum(2)  = ddsum(2) + ddR(k,2);
    ddR(k,3) = N(2,b+1)*M(2,c+1)*CP(i-p+b,j-q+c,4);
    ddsum(3)  = ddsum(3) + ddR(k,3);
  end
end

for k = 1:ne
  ddR(k,1) = ddR(k,1)/sum - 2*dR(k,1)*dsum(1)/sum^2 ...
             -R(k)*ddsum(1)/sum^2 + 2*R(k)*dsum(1)^2/sum^3;
  ddR(k,2) = ddR(k,2)/sum - 2*dR(k,2)*dsum(2)/sum^2 ...
             -R(k)*ddsum(2)/sum^2 + 2*R(k)*dsum(2)^2/sum^3;
  ddR(k,3) = ddR(k,3)/sum - dR(k,1)*dsum(2)/sum^2 - dR(k,2)*dsum(1)/sum^2 ...
             -R(k)*ddsum(3)/sum^2 + 2*R(k)*dsum(1)*dsum(2)/sum^3;
  
  dR(k,1) = dR(k,1)/sum - R(k)*dsum(1)/sum^2;
  dR(k,2) = dR(k,2)/sum - R(k)*dsum(2)/sum^2;
end