%Milton Orlando Sarria
%milton.sarria00@usc.edu.co
%papers:
%Milton Sarria-Paja, Tiago H. Falk, Fusion of auditory inspired amplitude modulation spectrum and cepstral features for whispered and normal speech speaker verification, Computer Speech & Language, Volume 45, 2017.
%Milton Sarria-Paja and Tiago H. Falk, Whispered Speech Detection in Noise Using Auditory-Inspired Modulation Spectrum Features. IEEE Signal Processing Letters, Vol. 20, No. 8, pp. 783-786, Aug. 2013. 
clear all
%AAMF: Auditory-inspired amplitude modulation features, this uses a
%frequency approach to compute the modulation components with the stft,
%represents the speech signal in 3-d array and reduces the dimensionality
%with pca.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%path where the algorithm searches for audio files and will save features
path_data='';

%EXAMPLE location: '/home/user/data/audio/' then path_data='/home/user/data/'
audio_data=[path_data 'audio/'];

%path where the algorithm saves the resulting features in hdf5 format, if
%the absolute path does not exist, it will be created
feats_data =[path_data 'features/'];
%% you can  have audio files and target feature files in a text file
% just put in one column the name of the audio file and in front of it, in a second column, the name
% of the target feature file separated with tab ('\t' character). 
%See the estructure of Audio2Feats.list file.
[wav_list,feat_list]=textread('Audio2Feats.list','%s\t%s');



%add toolboxes
addpath('voicebox');

%first it is necessary to compute the raw modulation spectrum
%representation. This is a high dimensional array (acoustic filters x
%modulation filters)
%(audio_data,dest_data,wav_list,feat_dest,FS)
B   = 80; %hz this is the maximum modulation domain frequency we want to see
ctx = 0.2;%context in seconds. After having the stft, we need to use time contexts to compute the mod spectrum
            %ctx is the length of the time context in seconds
%the window lentght and the overlap of the stft are selected on the basis of these two
%parameters (B and ctx)
%Here you should change the sampling rate to be used internally, based on
%this, the PCA file and number of filters are selected
FS=16e3;

computeFeatsRaw(audio_data,feats_data,wav_list,feat_list,FS,B,ctx)

%next, you can add any pos-processing, for instance I implemented a feature
%selection algorithm to discard channels. In this example I go straight  
%to compute the pca model, if you are going to use a big database, then 
%you should use a different algorithm to compute the pca model. we used an 
%iterative alorithm to compute the mean and the covariance matrix before 
%applying singular value decomposition SVD (I am attaching the algorithm in
%a separate script, computePCAmodel_itera)

if FS==16e3
    file_pca = 'PCA_16khz.h5';
elseif FS==8e3
    file_pca = 'PCA_8khz.h5';
elseif FS==4e3
    file_pca = 'PCA_4khz.h5';
end
%compute the PCA model
model_destination=[path_data 'modelsPCA/' file_pca];
computePCAmodel(feats_data,feat_list,model_destination)

%number of principal components to be selected
NC=40;
%% start the process to compute the aamf features using the PCA model
computeAAMF(audio_data,feats_data,wav_list,feat_list,FS,NC,B,ctx)




