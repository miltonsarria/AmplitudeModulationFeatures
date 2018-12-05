function [] = computePCAmodel(feat_data ,feat_list,model_destination)

field='/mods';
alog=1;

DATA=cell(length(feat_list),1);

for ii=1:length(feat_list) 
        
    x=h5read([feat_data feat_list{ii} '.h5'],field);
    [nf,nc,nfr]=size(x);
    x=reshape(x,nf*nc,nfr);      

    if alog
        x=log10(abs(x'));
    end
    
    DATA{ii}=x;
    clear x
end
%compute the mean and covariance matrix 
DATA=cell2mat(DATA);
MU=mean(DATA);
%center data
DATA=bsxfun(@minus, DATA, MU);
Sigma = cov(DATA);

%perform SVD
[U,S,V] = svd(Sigma);

%varac=cumsum(diag(S)/sum(diag(S)));
%varac=~(varac>0.999);
path_file = fileparts(model_destination);
if ~(exist(path_file,'dir')==7), unix(['mkdir -p ' path_file]); end  

h5create(model_destination,'/U',size(U)); 
h5write(model_destination, '/U', U);

h5create(model_destination,'/Mean',size(MU)); 
h5write(model_destination, '/Mean', MU);

h5create(model_destination,'/Sigma',size(Sigma)); 
h5write(model_destination, '/Sigma', Sigma);

h5create(model_destination,'/S',size(S)); 
h5write(model_destination, '/S', S);


