function [cp] = clusterprominence(Im)
%Im=rgb2gray(Im);

GLCM=graycomatrix(Im);
[x,y] = size(GLCM);
sum=0;
for i = 1:x 
    for j = 1:y
        sum=sum+GLCM(i,j);
    end
end
ux=0;
for i =1:x
    p(i)=0;
    for j =1:y
        p(i)=p(i)+(GLCM(i,j)/sum);
    end
    ux=ux+(i*p(i));
end
uy=0;
for j=1:y
    q(j)=0;
    for i = 1:x
        q(j) = q(j) + (GLCM(i,j)/sum);
    end
    uy=uy+(j*q(j));
end
cp=0;
for i=1:x
    for j=1:y
        cp=cp+(((i+j-ux-uy)^4)*(GLCM(i,j)/sum));
    end
end
end
