function [wss] = frq2wss(frq)
%FRQ2BARK  Convert Hertz to wss frequency scale 
%       wss = frq2wss(frq) converts a vector of frequencies (in Hz)
%       to the corresponding values on the wss scale.
% Inputs: f  vector of frequencies in Hz
%
% Outputs: wss  wss values
%REFERENCES
%201O international Conference on Information, Networking and Automation (ICINA)
%Noise Reduction in Whisper Speech Based On the Auditory Masking Model
%Zhi Tao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Milton Orlando Sarria
%INRS-EMT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fl=frq(frq<2e3);
fh=frq(frq>=2e3);
wss=zeros(size(frq));
wfl=(2478*fl.^4)./(1220^4 + fl.^4);
wfh=4100-2000./(1+exp((fh-3000)/310));
wss(frq<2e3)=wfl;
wss(frq>=2e3)=wfh;
if ~nargout
    plot(frq,wss,'-x');
    xlabel(['Frequency (' xticksi 'Hz)']);
    ylabel(['Frequency (' yticksi 'wss)']);
end
