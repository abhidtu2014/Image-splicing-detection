function [asm] = angularsecondmoment(Im)
%Im=rgb2gray(Im);

GLCM=graycomatrix(Im);
[x,y] = size(GLCM);
asm = 0;
sum=0;
for i = 1:x 
    for j = 1:y
        sum=sum+GLCM(i,j);
    end
end

for i = 1:x 
    for j = 1:y
        asm = asm+((GLCM(i,j)/sum)^2);
    end
end
end