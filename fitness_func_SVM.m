function [avg,class_result] = fitness_func_SVM(c)                  %total no. of images
%training_group = ones(1,19);
% Training feature matrix of class 1
for g=1:10
    str = 'AU_T_0';
     file=strcat(str,num2str(g),'.bmp');  
       y = imread(file);

        t=[];
             
%                 feature=WLD(y);
%                 [n,m]=size(feature);
%                  feature=reshape(feature,1,n*m);
%                  t= [t,feature];
                
%                  ed=Edge_detection(y);
%                  [n,m]=size(ed);
%                  ed=reshape(ed,1,n*m);
%                  t= [t,ed];
%                  
%                  spe=spatial_predict(y);
%                  [n,m]=size(spe);
%                  spe=reshape(spe,1,n*m);
%                  t= [t,spe];
                 
                ene = ILBP81ri(y);
                t= [ t, ene];

                ent = entropy(y);
                t= [ t, ent];

                con = Contrast_func(y);
                t= [ t, con];

               cs=hog_feature_vector(y);
                t= [ t, cs];

                idm = IDM_func(y);
                t= [ t, idm];

               in = inertia(y);
                t= [ t,in ];
 
                cp = clusterprominence(y);
                t= [ t, cp];

                asm = angularsecondmoment(y);
                t= [ t, asm];
        
        eval(['tr', num2str((g)),'= t;']);
end

for g=01:10
    str = 'SP_T_0';
     file=strcat(str,int2str(g),'.bmp');  
        y = imread(file);
        %imshow(a,map) % tmp is an indexed image file using 
        % color map "map"
        %y=double(rgb2gray(x));
        q=[];

%                 feature=WLD(y);
%                 [n,m]=size(feature);
%                  feature=reshape(feature,1,n*m);
%                  t= [t,feature];
                
%                  ed=Edge_detection(y);
%                  [n,m]=size(ed);
%                  ed=reshape(ed,1,n*m);
%                  t= [t,ed];
%                  
%                  spe=spatial_predict(y);
%                  [n,m]=size(spe);
%                  spe=reshape(spe,1,n*m);
%                  t= [t,spe];
                ene = ILBP81ri(y);
                q= [ q, ene];

                ent = entropy(y);
                q= [ q, ent];

                con = Contrast_func(y);
                q= [ q, con];
 
                cs=hog_feature_vector(y);
                q= [ q, cs];
  
                idm = IDM_func(y);
                q= [ q, idm];
  
                in = inertia(y);
                q= [ q,in ];

                cp = clusterprominence(y);
                q= [ q, cp];
  
                asm = angularsecondmoment(y);
                q= [ q, asm];
                
        eval(['tr', num2str((10+g)),'= q;']);
end
warning('off', 'Images:initSize:adjustingMag');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for g=01:10
    str = 'AU_TT_O_0';
     file=strcat(str,num2str(g),'.bmp');  
       y = imread(file);
       %imshow(x,map) 
       % tmp is an indexed image file using 
        % color map "map"
        %y=double(rgb2gray(x));
        t=[];
         
%                 feature=WLD(y);
%                 [n,m]=size(feature);
%                  feature=reshape(feature,1,n*m);
%                  t= [t,feature];
                
%                  ed=Edge_detection(y);
%                  [n,m]=size(ed);
%                  ed=reshape(ed,1,n*m);
%                  t= [t,ed];
%                  
%                  spe=spatial_predict(y);
%                  [n,m]=size(spe);
%                  spe=reshape(spe,1,n*m);
%                  t= [t,spe];
                ene = ILBP81ri(y);
                t= [ t, ene];

                ent = entropy(y);
                t= [ t, ent];

                con = Contrast_func(y);
                t= [ t, con];

                cs=hog_feature_vector(y);
                t= [ t, cs];

                idm = IDM_func(y);
                t= [ t, idm];

                in = inertia(y);
                t= [ t,in ];

                cp = clusterprominence(y);
                t= [ t, cp];

                asm = angularsecondmoment(y);
                t= [ t, asm];
       
        eval(['tr', num2str((20+g)),'= t;']);
end

for g=01:10
    str = 'SP_TT_O_0';
     file=strcat(str,int2str(g),'.bmp');  
        y = imread(file);
        %imshow(a,map) % tmp is an indexed image file using 
        % color map "map"
        %y=double(rgb2gray(x));
        q=[];
            
%                   feature=WLD(y);
%                 [n,m]=size(feature);
%                  feature=reshape(feature,1,n*m);
%                  t= [t,feature];
                
%                  ed=Edge_detection(y);
%                  [n,m]=size(ed);
%                  ed=reshape(ed,1,n*m);
%                  t= [t,ed];
%                  
%                  spe=spatial_predict(y);
%                  [n,m]=size(spe);
%                  spe=reshape(spe,1,n*m);
%                  t= [t,spe];
                 ene = ILBP81ri(y);
              
                q= [ q, ene];

                ent = entropy(y);
                q= [ q, ent];

                con = Contrast_func(y);
                q= [ q, con];

                cs=hog_feature_vector(y);
                q= [ q, cs];

                idm = IDM_func(y);
                q= [ q, idm];

                in = inertia(y);
                q= [ q,in ];

                cp = clusterprominence(y);
                q= [ q, cp];

                asm = angularsecondmoment(y);
                q= [ q, asm];
            
        
        eval(['tr', num2str((30+g)),'= q;']);
end
warning('off', 'Images:initSize:adjustingMag');

 for i=1:10
     training_group(1,i)=1;
 end
 for i=11:20
     training_group(1,i)=0;
 end
 
training=[];
for i =1:20
training = [training;eval(['tr',num2str(i),';'])];
end

testing = [];
for i =21:40
    testing = [testing;eval(['tr',num2str(i),';'])];
end

Train = training;
Test = testing;
group=training_group';
SVMStruct=svmtrain(Train,group);    %'method','svm'
class_result = svmclassify(SVMStruct, Test);

display(class_result);

%display(classperf(group));

%display(classperf(class_result)); 


TP=0;

for i=1:10
    if class_result(i)==1
        TP=TP+1;
    end
end

TN=0;
newcounter=10;
for i=1:10
    if class_result(i+newcounter)==0
        TN=TN+1;
    end
end

FP=10-TP;
FN=10-TN;
Sensitivity=(TP/(TP+FN))*100;
Specificity=(TN/(TN+FP))*100;
Accuracy=((TP+TN)/(TP+TN+FP+FN))*100;
display(Accuracy);
display(Specificity);
display(Sensitivity);
%ne=FP+FN;
avg = (Sensitivity + Specificity + Accuracy)/3;

class_result = class_result';
%display(class_result);

% ans =
end
