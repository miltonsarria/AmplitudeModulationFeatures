function [frq] = wss2frq(wss)
%FRQ2BARK  Convert wss  frequency scale to Hertz
%       frq = wss2frq(wss) 
%       to the corresponding values on the wss scale.
% Inputs: ww  vector of frequencies in wss
%
% Outputs: frq  hertz values
%REFERENCES
%201O international Conference on Information, Networking and Automation (ICINA)
%Noise Reduction in Whisper Speech Based On the Auditory Masking Model
%Zhi Tao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Milton Orlando Sarria
%INRS-EMT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lim=frq2wss(2000);
a=2478; b=1220^4;

frq=zeros(size(wss));
wl=wss(wss<lim); 
frq(wss<lim)=((-wl*b)./(wl-a)).^(1/4);

wh=wss(wss>=lim); 

a=4100; b=2000; c=3000; d=310;
frq(wss>=lim)=log((wh+b-a)./(a-wh))*d+c;


if ~nargout
    plot(wss,frq,'-x');
    xlabel(['Frequency (' xticksi 'wss)']);
    ylabel(['Frequency (' yticksi 'Hz)']);
end
