function [fexp] = frq2exp(frq)
%FRQ2exp  Convert Hertz to exponencial scale exp=(FRQ)
%	The relationship between exp and frq is given by:
%
%	fexp = 6375*(10^(f/50000)-1) 
%
%  	This means that m(1000) = 1000
%
%	References:
%
%[6] X. Fan and J. Hansen, "Speaker identi?cation within whispered speech audio streams", Audio,
%Speech, and Language Processing, IEEE Transactions on, vol. 19, pp. 1408 �1421, july 2011.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
af=abs(frq);
%fs2=4e3;

c=10610;
k=50e3;

%c=(1127*log(1+fs2/700))/(10^(fs2/k)-1);
%k=fs2/log10((2595*log10(1+fs2/700)+c)/c);



fexp = sign(frq).*(c*(10.^(af/k)-1)); 
if ~nargout
    plot(frq,fexp,'-x');
    xlabel(['Frequency (' xticksi 'Hz)']);
    ylabel(['Frequency (' yticksi 'Mel)']);
end
