%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Course: Biometrics
%
%Title: First lab experience - Fingerprint matching
%
%Author: prof. Simone Milani (simone.milani@dei.unipd.it)
%
%Check the correspondence between couples of fingeprints
%Using FingerPrint Demo software by Florence Kussener, 2007
%https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/16728/versions/5/previews/FingerPrint/html/fingerprint.html#21
%Bug were fixed and software was modified in order to save FP minutiae in
%mat data

%clear workspaces
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINGERPRINT 1 (Query)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load database/A2.mat  %minutiae
I1=imread('database/A.bmp');  %Image
%load termination
min_term1=MinutiaFin;  
min_term1= [ min_term1 min_term1(:,3) min_term1(:,3) ];  %padding the final
                                                         %angle to make it
                                                         %compatible with
                                                         %bifurcation
%load bifurcation
min_bif1=MinutiaSep;
points1=[min_term1(:,1:5) ; min_bif1(:,1:5) ]; %add at the end
ind_type1=[zeros(size(min_term1,1),1) ; ones(size(min_bif1,1),1)  ];
ind=find(sum(isnan(points1),2)==0);  %select valid points  (not NaN)
points1=points1(ind,:);
N1=size(points1,1); %total minutiae number
ind_type1=ind_type1(ind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINGERPRINT 2 (Template)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load database/archive/00000.mat;   %minutiae
I2=imread('database/archive/00000.bmp');  %image
%load termination
min_term2=MinutiaFin;
min_term2= [ min_term2 min_term2(:,3) min_term2(:,3) ]; %paddiing (see above)
%load bifurcation
min_bif2=MinutiaSep;
points2=[min_term2(:,1:5) ;  min_bif2(:,1:5) ];  %concatenate
ind_type2=[zeros(size(min_term2,1),1) ; ones(size(min_bif2,1),1)  ];
ind=find(sum(isnan(points2),2)==0);   %find valid points
points2=points2(ind,:);
N2=size(points2,1);  %number of minutiae
ind_type2=ind_type2(ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Visualize FPs and Minutiae
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(2,2,1);  %first subplot
imshow(I1);
hold on;
plot(min_term1(:,1),min_term1(:,2),'bo',min_bif1(:,1),min_bif1(:,2),'gv');
hold off;
subplot(2,2,2);  %second subplot
imshow(I2);
hold on;
plot(min_term2(:,1),min_term2(:,2),'bo',min_bif2(:,1),min_bif2(:,2),'gv');
hold off;

%center the fingerprints (subtract the average minutiae coordinates)
points1(:,1:2)=points1(:,1:2)-ones(N1,1)*mean(points1(:,1:2));
points2(:,1:2)=points2(:,1:2)-ones(N2,1)*mean(points2(:,1:2));

%%%%%%%%%%%%%%%%%%%%%%
%HOUGH TRANSFORM
%%%%%%%%%%%%%%%%%%%%%%

%PARAMETER RANGE  (to be optimized)
Deltax=[-300:10:300];
Deltay=[-300:10:300];
DeltaT=[-5.8:0.05:5.8];

%Accumulator tensor
Amat=zeros(length(Deltax),length(Deltay),length(DeltaT));

for i=1:N1
    for j=1:N2
        %compute the distance between the angles
        d1=abs(points1(i,3:5) - points2(j,3:5));  %dist1
        d2=abs(points1(i,3:5) - (2*pi-points2(j,3:5))); %dist2 (360-angle)
        d12=min([d1; d2]);  %find minimum
        
        dT=mean(d12,2) ;  %average all over the angles
     
        %Delta x
        dx=points1(i,1)-points2(j,1)*cos(dT) + points2(j,2)*sin(dT);
        %Delta y
        dy=points1(i,2)-points2(j,1)*sin(dT) - points2(j,2)*cos(dT);

        %find the correponding coordinates into the accumulator
        [~,iT]=min(abs(dT*ones(1,length(DeltaT))-DeltaT));
        [~,ix]=min(abs(dx*ones(1,length(Deltax))-Deltax));
        [~,iy]=min(abs(dy*ones(1,length(Deltay))-Deltay));
        %update votes
        Amat(ix,iy,iT)=Amat(ix,iy,iT)+1;
    end;
end;

%Maximum frequency
[Mminutiae,idx]=max(Amat(:));
%Find index for the max value in the accumulator tensor
best_dx=mod((idx-1),length(Deltax))+1;
best_dy=mod(floor((idx-1)/length(Deltax)),length(Deltay))+1;
best_dT=mod(floor((idx-1)/(length(Deltax)*length(Deltay))),length(DeltaT))+1;

%Paramaters DeltaTheta, DeltaX, DeltaY for the most voted configuration
DT=DeltaT(best_dT);
DX=Deltax(best_dx);
DY=Deltay(best_dy);

%show the results
disp([DX DY DT Mminutiae]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POINTS MATCHING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rotation matrix
R=[cos(DT) -1*sin(DT) ;  sin(DT) cos(DT) ];

%match minutiae coordinates
points2(:,1:2)=points2(:,1:2)*R-ones(N2,1)*[DX DY];
points2(:,3:5)=points2(:,3:5)+DT;

%rotate and translate image
I2r=imrotate(I2,DT/pi*180);
I2m=imtranslate(I2r,[DX DT ]);

%Visualize 
figure(1);
subplot(2,2,3);  
imshow(I2m);

%Visualize 
figure(1);
subplot(2,2,4);  %plot rotated points
plot(points1(:,1),points1(:,2),'bo',points2(:,1),points2(:,2),'gv');


%counter (to check if pints were assigned or not)
fp1=zeros(N1,1);
fp2=zeros(N2,1);

%%%%%%%%%%%%THRESHOLDS (to be optimized) %%%%%%%%%
thD=10;      %th on distance
thA=0.05;     %th on angles

for i=1:N1   %loop over Query minutiae
    %find distance with respect to coords of template  minutiae 
    dplace=sqrt(mean((ones(N2,1)*points1(i,1:2)-points2(:,1:2)).^2,2));
    %find angle distance with respect to orientation of template  minutiae 
    dangle1=mean(abs(ones(N2,1)*points1(i,3:5)-points2(:,3:5)),2);
    dangle2=mean(abs(ones(N2,1)*(2*pi-points1(i,3:5))-points2(:,3:5)),2);
    dangle=min([dangle1' ; dangle2']);
    
    %check points that are lower than thresholds
    idplace=find(dplace<thD);
    idangle=find(dangle<thA);
    ival=intersect(idplace,idangle);  %index of valid points
    ival=intersect(ival,find(fp2==0));  %check if minutiae were not assigned
    ival=intersect(ival,find(ind_type2==ind_type1(i)));  %check if minutiae were not assigned
    
    
    if length(ival)>0
        [~,ir]=min(dplace(ival));  %if more than one candidate, choose the closest
        i1=ival(ir);
        fp1(i)=i1;  %set counters to 1
        fp2(i1)=1;
    end;
end;


%Display the score
score=sum(fp2)/((N1+N2)/2);
disp(score);
disp('--------')
