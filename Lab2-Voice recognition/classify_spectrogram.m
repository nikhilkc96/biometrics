%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Course: Biometrics A.A. 2020/2021
%SECOND LAB EXPERIENCE - VOICE RECOGNITION
%author: Simone Milani (simone.milani@dei.unipd.it)
%May 20, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

%labels
vet_voices={ 'bibi' 'madiba' 'jens' 'julia' 'margie'};

Nclasses=5;

%Trainig set
%find the minimum number of trainig samples per class (to balance)
ns=100000000;
for c=1:Nclasses
    eval(sprintf('load %s_train',vet_voices{c}));
    if (size(Xsound,3)<ns)
        ns=size(Xsound,3);
    end;
end;

%load the training samples
Xtrain=zeros(32,32,1,ns*5);
val_train=zeros(ns*Nclasses,1);
cnt=1;
for c=1:5
    eval(sprintf('load %s_train',vet_voices{c}));  
    x=reshape(Xsound,[32 32 1 size(Xsound,3)]);  %reshape the array of spectrograms
                                                 %into 32x32x1xns                 
    Xtrain(:,:,1,cnt:cnt+ns-1)=x(:,:,1,1:ns);
    cnt=cnt+ns;
    val_train(cnt:cnt+ns-1)=c;  %labels
end;
%remove last samples
val_train= val_train(1:cnt-1); 
Xtrain=Xtrain(:,:,1,1:cnt-1);

val_train=categorical(val_train);  %convert into categorical


%scramble data
ind=randperm(cnt-1);
Xtrain=Xtrain(:,:,1,ind);
val_train=val_train(ind);


%Validation set
%find the minimum number of trainig samples per class (to balance)
ns=100000000;
for c=1:Nclasses
    eval(sprintf('load %s_valid',vet_voices{c}));
    if (size(Xsound,3)<ns)
        ns=size(Xsound,3);
    end;
end;

%load the training samples
Xvalid=zeros(32,32,1,ns*5);
val_valid=zeros(ns*5,1);
cnt=1;
for c=1:Nclasses
    eval(sprintf('load %s_valid',vet_voices{c}));  
    x=reshape(Xsound,[32 32 1 size(Xsound,3)]);  %reshape the array of spectrograms
                                                 %into 32x32x1xns                 
    Xvalid(:,:,1,cnt:cnt+ns-1)=x(:,:,1,1:ns);
    cnt=cnt+ns;
    val_valid(cnt:cnt+ns-1)=c;  %labels
end;
%remove last samples
val_valid= val_valid(1:cnt-1); 
Xvalid=Xvalid(:,:,1,1:cnt-1);

val_valid=categorical(val_valid);  %convert into categorical
 
%Test set
%find the minimum number of test samples per class (to balance)
ns=100000000;
for c=1:Nclasses
    eval(sprintf('load %s_test',vet_voices{c}));
    if (size(Xsound,3)<ns)
        ns=size(Xsound,3);
    end;
end;

%load test
Xtest=zeros(32,32,1,ns*5);
val_test=zeros(ns*5,1);
cnt=1;
for c=1:Nclasses
    eval(sprintf('load %s_test',vet_voices{c}));
    x=reshape(Xsound,[32 32 1 size(Xsound,3)]); %reshape the array of spectrograms
                                                 %into 32x32x1xns
    Xtest(:,:,1,cnt:cnt+ns-1)=x(:,:,1,1:ns);
    cnt=cnt+ns;
    val_test(cnt:cnt+ns-1)=c;
end;
%remove last samples
val_test= val_test(1:cnt-1);
Xtest=Xtest(:,:,1,1:cnt-1);

val_test=categorical( val_test); %convert to categorical

%optimization parameters
miniBatchSize = 256;  %minibatch size
numValidationsPerEpoch = 2;  %number of validation per epoch
validationFrequency = floor(size(Xtrain,4)/miniBatchSize/numValidationsPerEpoch);

options = trainingOptions('sgdm',...
    'ExecutionEnvironment','cpu',...  %use since no GPU available
    'MaxEpochs',3, ...
    'ValidationData',{Xvalid,val_valid'},...  %(optional) validate the CNN during training
    'ValidationFrequency',30,...
    'VerboseFrequency',30,...
    'Verbose',true,...  %display messages
    'Plots','training-progress'); %display plots


%%%%%%%%%%%%%%TO DO 2 %%%%%%%%%%%%%%%%%
%modify the spectrum here


layers = [
    imageInputLayer([32 32 1])   %first input layer (image size 28 x 28)

    convolution2dLayer(3,8,'Padding',1)  %conv layer
    batchNormalizationLayer  %normalize
    reluLayer  %ReLU
    
    convolution2dLayer(3,16,'Padding',1)  %conv layer
    batchNormalizationLayer  %normalize
    reluLayer  %ReLU

    maxPooling2dLayer(2,'Stride',2)  %max pooling

    convolution2dLayer(3,64,'Padding',1)
    batchNormalizationLayer
    reluLayer

    maxPooling2dLayer(2,'Stride',2)

    convolution2dLayer(3,128,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
%     maxPooling2dLayer(2,'Stride',2)
% 
%     convolution2dLayer(3,128,'Padding',1)
%     batchNormalizationLayer
%     reluLayer

    
    fullyConnectedLayer(Nclasses)  %final layer
    softmaxLayer   %softmax
    classificationLayer];  %classify (take the output values with the highest prob.
 
%%%%%%%%%%%%%%%%%%%%%%%%END of modification%%%%%%%%%%%%%%%%%%%%%%

rng('default') % For reproducibility
net = trainNetwork(Xtrain,val_train,layers,options);  %train the network

%evaluate the accuracy
label_out = classify(net,Xtest,'ExecutionEnvironment','cpu','MiniBatchSize',miniBatchSize);

%average precision
accuracy = sum(label_out(:) == val_test(:))/numel(val_test)

%generate confusion matrix
conf_mat=zeros(Nclasses,Nclasses);
for c1=1:Nclasses
    ind1=find(double(val_test(:))==c1);  %select values that belong to c1 classes
    for c2=1:Nclasses
        ind2=find(double(label_out(:))==c2);  %elements classified as c2
        ind=intersect(ind1,ind2);  %c1 elements classified as c2
        conf_mat(c1,c2)=length(ind);  %save the number
    end;
end;
%normalize the numbers into percentages
conf_mat=conf_mat./(sum(conf_mat,2)*ones(1,Nclasses))*100;
disp(conf_mat);

plotconfusion(val_test, label_out)
