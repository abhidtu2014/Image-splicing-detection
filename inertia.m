function [in] = inertia(Im)
%Im=rgb2gray(Im);
GLCM=graycomatrix(Im);
[x,y] = size(GLCM);
sum=0;
in=0;
for i = 1:x 
    for j = 1:y
        sum=sum+GLCM(i,j);
    end
end
for i=1:x
    for j=1:y
        in = in + (((i-j)^2)*(GLCM(i,j)/sum));
    end
end
        