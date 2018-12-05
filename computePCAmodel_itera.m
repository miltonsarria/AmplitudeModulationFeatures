function [] = computePCAmodel_itera(feat_data ,feat_list,model_destination)

n=0;
field='/mods';
alog=1;
%%%%%%%%%%%%%
%we need to estimate the mean an covariance matrix
for ii=1:length(feat_list) 
    
        x=h5read([feat_data feat_list{ii} '.h5'],field);
        [nf,nc,nfr]=size(x);
        x=reshape(x,nf*nc,nfr);      
        if alog
            x=log10(abs(x'));
        end
              
        n=n+nfr;
        
        if ii==1
            sigma_hat_n  = cov(x);
            mean_hat_n   = mean(x); mean_hat_n=mean_hat_n';
            N=nfr;
            
        else
           for jj=1:size(x,1)
                y=(x(jj,:))';               
                mean_hat_1 =mean_hat_n + 1/(N+1)*(y-mean_hat_n);
                sigma_hat_1=(N-1)/N*sigma_hat_n + 1/(N+1)*(y-mean_hat_n)*(y-mean_hat_n)';
                %update the old parameters
                N=N+1;
                mean_hat_n=mean_hat_1;
                sigma_hat_n=sigma_hat_1;
           end
        end
    clear x
end
%Having the covariance matrix and the mean, we can perform svd and save
%results
MU=mean_hat_n';
Sigma=sigma_hat_n;
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



