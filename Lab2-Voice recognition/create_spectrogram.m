%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Course: Biometrics A.A. 2020/2021
%SECOND LAB EXPERIENCE - VOICE RECOGNITION
%author: Simone Milani (simone.milani@dei.unipd.it)
%May 20, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

%shorter name to be used in the dataset
vet_voices={ 'bibi' 'madiba' 'jens' 'julia' 'margie'};

%name of the directories containing voice files
vet_voices_dir={ 'Benjamin_Netanyau' 'Nelson_Mandela' 'Jens_Stoltenberg' ...
   'Julia_Gillard' 'Margaret_Tatcher'  };

%flag to display graph
flag_display=0;

%parameter for FFT values
Nfft=128;  %length of FFT
woverlap=64;  %overlap between windoes (in samples)
nspect=32; %side of the spectrogram to be generated

minVal=-20;   %minimum value for FFT coeffs
maxVal=33;  %max values for 10*log10 FFT coeffs

for c=1:5  %loop on each voice/folder
    
    eval(sprintf('vet_files=dir(''16000_pcm_speeches/%s/*.wav'');',...
        vet_voices_dir{c}));
    
    Nsound0=length(vet_files); %#of wav files
    
    %Nsound0=min(Nsound0,500);   %enable to use less data!
    
    %divide into traiing/validation/testing files
    ind0=randperm(Nsound0);  
    ntrain=round(Nsound0*0.35);  %size of traiing set
    nvalid=round(Nsound0*0.3);   %size of validation
    ntest=Nsound0-ntrain-nvalid; %size of testing
    
    %%%%%%TO DO 3 %%%%
    %you can change the percentage of train/test/validation
    %%%%%%%%%%%%%%%%%%
    
    
    
    for run=0:2   %loop for training/validation/testing
        if run==0  %training
            ind=ind0(1:ntrain);
            Nsound=length(ind);
        elseif run==1    %validation
            ind=ind0(ntrain+1:ntrain+nvalid);
            Nsound=length(ind);
        else    %testing
            ind=ind0(ntrain+nvalid+1:end);
            Nsound=length(ind);
        end;
    
        %array where spectrograms are to be stored
        Xsound=zeros(32,32,Nsound*9);
    
        cnt_s=1;
        for f=1:Nsound    %loop over files
            %file reading
            str_sound=[ vet_files(ind(f)).folder '/' vet_files(ind(f)).name ];
            [x,Fs]=audioread(str_sound);
        
            %x: array
            %Fs: sampling frequency
            
            Nsamples=length(x);  %#samples
            
            wSpect=length(1:woverlap:Nsamples-Nfft); %array of initial steps
            
            s=zeros(Nfft,wSpect);   %stores the FFT
            
            cnt=1;
            for istep=1:woverlap:Nsamples-Nfft  %scan the file
                
                frame=x(istep:istep+Nfft-1);  %time frame
                
                F=20*log10(abs(fft(frame)));  %FFT
                
                F(find(F<minVal))=minVal;  %remove very negative values
                
                
                %%%%%%%%TO DO 1 %%%%%%%%%
                
                %modify the code here
                
                %%%%%%%%%%%%%%%%%%%%%%%%%
                
                %store FFT
                s(:,cnt)=F';
                cnt=cnt+1;
            end;
        
            if flag_display
            %plot the FFT coefficients
                figure(1);
                imagesc(s);
                title(sprintf('%05d',f));
                colorbar;
            end;
            
            %compute power in the window (select relevant parts of
            %spectrogram)
            for istep=1:nspect:wSpect-nspect
                avg_pow=mean(mean(s(1:nspect,istep:istep+nspect-1)));
                if (avg_pow>0.004)  %if power > 0.004 - store it
                    Xsound(:,:,cnt_s)=s(1:nspect,istep:istep+nspect-1);
                    cnt_s=cnt_s+1;
                end;
            end;
            
        end;
        
        %crop the set of valid spectrograms
        Xsound=Xsound(:,:,1:cnt_s-1);
        
        %store the different values
        if (run==0)
            eval(sprintf('save %s_train Xsound',vet_voices{c}));
        elseif (run==1)
            eval(sprintf('save %s_valid Xsound',vet_voices{c}));
        else
            eval(sprintf('save %s_test Xsound',vet_voices{c}));
        end;
    end;
end;