function R = refinement_vec_repeated(U,ref,p)

R=[];
k = 1;
if (p==2)
for i=1:length(U)-1
  if (U(i)~=U(i+1))
    for j = 1:ref-1;
      R(k) = j/ref*(U(i+1)-U(i))+U(i);
      R(k+1) = j/ref*(U(i+1)-U(i))+U(i);
      k =k+2;
    end 
  end
end
end

if (p==3)
for i=1:length(U)-1
  if (U(i)~=U(i+1))
    for j = 1:ref-1;
      R(k) = j/ref*(U(i+1)-U(i))+U(i);
      R(k+1) = j/ref*(U(i+1)-U(i))+U(i);
      R(k+2) = j/ref*(U(i+1)-U(i))+U(i);
      k =k+3;
    end 
  end
end
end 

if (p==4)
  for i=1:length(U)-1
  if (U(i)~=U(i+1))
    for j = 1:ref-1;
      R(k) = j/ref*(U(i+1)-U(i))+U(i);
      R(k+1) = j/ref*(U(i+1)-U(i))+U(i);
      R(k+2) = j/ref*(U(i+1)-U(i))+U(i);
      R(k+3) = j/ref*(U(i+1)-U(i))+U(i);
      k =k+4;
    end 
  end
  end
end