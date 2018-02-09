function [ idm ] = IDM_func( Im )
%Im=rgb2gray(Im);

% This function calculates the Inverse Difference Moment feature of a complete
% image matrix where Im is the input image matrix and m & n are the sizes of 
idm = 0;  % idm is the Inverse Difference Moment feature of the image
sum=0;

GLCM=graycomatrix(Im);
[x,y] = size(GLCM);
for i = 1:x 
    for j = 1:y
        sum=sum+GLCM(i,j);
    end
end

for i = 1:x
    for j = 1:y
            idm = idm +  (( GLCM(i,j) / ( 1 +((i-j)^2) ))/sum);
    end
end
end

