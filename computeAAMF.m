function []=computeAAMF(audio_data,dest_data,wav_list,feat_dest,FS,NC,B,ctx)
                           
%Milton Orlando sarria Paja
%INRS-EMT

%audio_data is the path where you can find the audio files 
%dest_data is the path to save the h5 files with the features
%wav_list is a list of audio files that are in  audio_data
%feat_dest is a list of names to save the h5 files
%file_pca is a file that contains a PCA projection already trained for MODS
%FS sampling rate, default 16e3
%NC = number of components default 40

if ~exist('B','var'),    B     = 80; end
if ~exist('ctx','var'),  ctx   = 0.2; end
if ~exist('FS','var'),   FS    = 16e3; end
if ~exist('NC','var'),   NC    = 40; end
%16 khz use 27 filters
%8 khz  use 20 filters
%4khz   use 17 filters

lowF=10; 
if FS==16e3
    file_pca = 'modelsPCA/PCA_16khz.h5';
    nFilters=27; 
    hiF=7000;
elseif FS==8e3
    file_pca = 'modelsPCA/PCA_8khz.h5';
    nFilters=20;   
    hiF=3528;
elseif FS==4e3
    file_pca = 'modelsPCA/PCA_4khz.h5';
    nFilters=17;
    hiF=2000;
end

fL = lowF/FS; 
fH = hiF/FS; 
nChan = nFilters; 

%16 khz we use 27 filters
%8 khz we use 20 filters
%4khz we use  17 filter
preemph = 0.97;
%max modulation domain 40 hz, with 200 ms context 
%% read list of files and start processing
if ~isempty(file_pca)
    U  = h5read(file_pca,'/U'); W=U(:,1:NC);
    MU = h5read(file_pca,'/Mean');
else
    feats(2)=0;
end
%% read list of files and start processing
for rr=1:length(wav_list)    
    
    name_wav    =[audio_data wav_list{rr}]; disp(name_wav)
    [s,fs]      =audioread(name_wav);         
    %change sampling rate of audio file
    if FS~=fs
        [b,a]   = cheby2(9,20,(FS)/(fs));
        s       =filter(b,a,s);
        s       =resample(s, FS, fs); %remuestrear              
    end    
   
   s=s/max(abs(s));
   s = filter([1 -preemph], 1, s);

%save features to a hdf5 file
   
   feat_file=[dest_data feat_dest{rr},'.h5']; disp(feat_file) 
   path_file = fileparts(feat_file);
   if ~(exist(path_file,'dir')==7), unix(['mkdir -p ' path_file]); end  
   %if any(VAD)
 
        msf=aamf(s,FS,B,ctx,fL,fH,nChan); 
        [nf,nc,nfr]=size(msf);
        if nfr~=0
            mx= reshape(msf,nf*nc,nfr); mx=mx';         
            mx=log10(abs(mx));
            mx=bsxfun(@minus, mx, MU);
            mx=mx*W;
            h5create(feat_file,'/aamf',size(mx)); 
            h5write(feat_file, '/aamf', mx);
            
        else    
            h5create(feat_file,'/aamf',size(0)); 
            h5write(feat_file, '/aamf', 0);
        end
  
   %end         
   
         
end
%%
function [Ei,fa,fm]=aamf(s,fs,B,M,fl,fh,nf)

%This function computes the modulation spectrum using the stft.
%Inputs:
%       s   input signal
%       fs  sampling rate, ( 16khz)
%       B   Hz, this is the bandwidth to determine what is the max
%       modulation frequency we want to detect (63 hz by default)
%       M   sec. this is the length of the windows in the freq domain (0.1
%       secs)
%       fl  lower frequency  to limit the range of freq, (0) 
%       fh  higher frequency (0.5 (fs/2) )
%Outputs:
%       Ei  detailed mod spectrum
%       fa acoustic frequency
%       fm modulation fequency
%B=63; M=0.2; fl=0; fh=0.5; nf=24; out=0;  s=filter([1 -0.97], 1, s);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('B','var'),   B   =63; end
if ~exist('M','var'),   M   =0.1; end
if ~exist('fl','var'),  fl  =0; end
if ~exist('fh','var'),  fh  =0.5; end
if ~exist('nf','var'),  nf  =27; end
%B=48
nf_m=8;                         %number of filters modulation domain
L=ceil(2*fs/B);                 %window size for the given B
nfft=1024;                      %number of fft bins (nfft/2)
fr=2*B;                         %frame rate, frames per second
inc1=ceil((fs-L)/(fr-1));       %calculate the increment 
sf=enframe(s,hamming(L),inc1);  %enframe the signal

%%%% delete silence frames, short utterances keep silence
%[~, ~, params]      = VQVAD; % Get default parameters to "params" struct
%params.frame_len    = L/fs;   % Frame duration, seconds
%params.frame_shift  = inc1/fs;   % Frame duration, seconds
%[speechInd]         = VQVAD(s, fs, params); VAD=speechInd==1;
%sf = sf(VAD,:);
%%%%

S=rfft(sf,nfft,2);              %real fft  
%process in the modulation domain
%compute window size and overlap

n=fix(M*fr); 
inc2=ceil(1/4*n);  
nfft2=256;                      %128 bins in the modulation domain
%%% group the bins in channels using filters
S=S';                                       %nfft/2+1 x number of frames, each column is a frame, each row a fft bin
[m,a,b,c]=melbankm(nf,nfft,fs,fl,fh,'u');   %filter bank
fa=a; a=b; b=c;                           
pw=S(a:b,:).*conj(S(a:b,:));                %power domain
pth=max(pw(:))*1E-20;  
mS=max(m*abs(S(a:b,:)),sqrt(pth))';  %
%ym=log(ym);
%process in the modulation domain
%to know the size of the output then the enframe is done for the first bin
sff=enframe(mS(:,1),hamming(n),inc2);           %inc1*inc2/fs
E2=zeros(nf, nfft2/2+1,size(sff,1));     
%calculate mod spectrum using nfft2 bins
for i=1:size(mS,2) %across columns - dft bins or filters
    sff=enframe(mS(:,i),hamming(n),inc2); %# frames x n points (M secs)
    Sf=rfft(sff,nfft2,2);                   %# frames (i-th bin) x nfft2/2 points (fft bins in mod frequency)   
    %Sf=abs(Sf);            %magnitude spectrum of i-th fft bin, |Sf|, each row a frame
     for j=1:size(Sf,1)
         E2(i,:,j)=Sf(j,:);
     end
end

fl=1/fr; fh=0.5;
[m,a,b,c]=melbankm(nf_m,nfft2,fr,fl,fh,'ul'); %filter bank 2 for mod spec
fm=a; a=b; b=c; %fc=find(10.^fc<20);
%fa=mel2frq(fa); %fai=find(fa>1e3);
Ei=zeros(nf, nf_m,size(sff,1));     
for i=1:size(E2,3)
    Sm=E2(:,:,i)';        %fft bins x number of filters in acoustic domain
    pm=Sm(a:b,:).*conj(Sm(a:b,:)); 
    pth=max(pm(:))*1E-20;   
    Ei(:,:,i)=max(m*pm,sqrt(pth))';    
end
%display it
fa=mel2frq(fa);
fm=10.^fm;
if nargout < 1
   surf(fm,fa,mean(Ei,3),'edgecolor','none'); axis tight; view(0,90);
end
