function [PS]=spatial_predict(image)

[rows,cols]=size(image);

PS=uint8(zeros(rows,cols)); % Predicted Image initialisation

for i=1: rows
    for j=1:cols
        if(j>1)&&(i>1)
            A=image(i-1,j);
            B=image(i,j-1);
            C=image(i-1,j-1);
            
            if C>=(max(A,B))
                PS(i,j)=uint8(min(A,B));
            elseif C<=min(A,B)
                PS(i,j)=uint8(max(A,B));
            else
                PS(i,j)=uint8(A+(B-C));
            end
            
        end
    end
end

for i=1: rows
    for j=1:cols
        e1(i,j)=image(i,j)-PS(i,j); %Calculate error
    end
end


end