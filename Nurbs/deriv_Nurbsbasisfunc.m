function [R,dR] = deriv_Nurbsbasisfunc(p,i,u,U,q,j,v,V,CP)
 

if (i==0); i = findspan(u,U,length(CP(:,1,1)));  end
if (j==0); j = findspan(v,V,length(CP(1,:,1)));  end

ne = (p+1)*(q+1);       
N = deriv2(i,p,u,U);     
M = deriv2(j,q,v,V);     

 
k = 0;
sum = 0;
dsum = zeros(2,1);

for c = 0:q
  for b = 0:p
    k = k+1;
  
    R(k) = N(1,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    sum = sum + R(k);
            
 
    dR(k,1) = N(2,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    dsum(1)  = dsum(1) + dR(k,1);
    dR(k,2) = N(1,b+1)*M(2,c+1)*CP(i-p+b,j-q+c,4);
    dsum(2)  = dsum(2) + dR(k,2);
  end
end

 
for k = 1:ne
    dR(k,1) = dR(k,1)/sum - (R(k)*dsum(1))/(sum^2);
    dR(k,2) = dR(k,2)/sum - (R(k)*dsum(2))/(sum^2);
end