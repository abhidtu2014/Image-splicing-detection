function [feature] = WLD(I);

[cA1,cH1,cV1,cD1] = dwt2(I,'bior3.7');

A1 = upcoef2('a',cA1,'bior3.7',1);
H1 = upcoef2('h',cH1,'bior3.7',1); 
V1 = upcoef2('v',cV1,'bior3.7',1);
D1 = upcoef2('d',cD1,'bior3.7',1);
 
% subplot(2,2,1); image(wcodemat(A1,192));
% title('Approximation A1')
% subplot(2,2,2); image(wcodemat(H1,192));
% title('Horizontal Detail H1')
% subplot(2,2,3); image(wcodemat(V1,192));
% title('Vertical Detail V1') 
% subplot(2,2,4); image(wcodemat(D1,192));
% title('Diagonal Detail D1')

feature = idwt2(cA1,cH1,cV1,cD1,'bior3.7');



