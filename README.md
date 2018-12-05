# AmplitudeModulationFeatures
AAMF: Auditory-inspired amplitude modulation features

Milton Orlando Sarria

Universidad Santiago de Cali


## papers:
Milton Sarria-Paja, Tiago H. Falk, Fusion of auditory inspired amplitude modulation spectrum and cepstral features for whispered and normal speech speaker verification, Computer Speech & Language, Volume 45, 2017.

Milton Sarria-Paja and Tiago H. Falk, Whispered Speech Detection in Noise Using Auditory-Inspired Modulation Spectrum Features. IEEE Signal Processing Letters, Vol. 20, No. 8, pp. 783-786, Aug. 2013. 


this code uses a frequency approach to compute the modulation components with the stft,
represents the speech signal in 3-d array and reduces the dimensionality with pca.

you can  have audio files and target feature files in a text file
just put in one column the name of the audio file and in front of it, in a second column, the name
of the target feature file separated with tab ('\t' character). 
See the estructure of Audio2Feats.list file.

first it is necessary to compute the raw modulation spectrum representation. This is a high dimensional array (acoustic filters x modulation filters)

B   = 80; hz this is the maximum modulation frequency we want to see
ctx = 0.2;context in seconds. After having the stft, we need to use time contexts to compute the mod spectrum. ctx is the length of the time context in seconds

the window lentght and the overlap of the stft are selected on the basis of these two
parameters (B and ctx)

next, you can add any pos-processing, for instance we implemented a feature
selection algorithm to discard channels. In this example I go straight  
to compute the pca model, if you are going to use a big database, then 
you should use a different algorithm to compute the pca model. we used an 
iterative alorithm to compute the mean and the covariance matrix before 
applying singular value decomposition SVD (computePCAmodel_itera.m)


run example_MAIN.m 


