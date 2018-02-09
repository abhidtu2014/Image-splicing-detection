function [ con ] = Contrast_func( Im)
%Im=rgb2gray(Im);

% This function calculates the contrast feature of a complete image matrix 
% where Im is the input image matrix and m & n are the sizes of Im
con = 0;                        % con is the contrast feature of the image
GLCM=graycomatrix(Im);
[x,y] = size(GLCM);
sum=0;
for i = 1:x 
    for j = 1:y
        sum=sum+GLCM(i,j);
    end
end

% Calculation of contrast feaure for the image
for i = 1:x
    for j = 1:y
       con=con+((GLCM(i,j)/sum)*((abs(i-j))^2));
    end
end
end