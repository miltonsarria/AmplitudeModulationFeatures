function [frq] = exp2frq(mel)
%exp2FRQ  Convert exp frequency scale to Hertz FRQ=(MEL)
%	The relationship between exp and frq is given by:
%
%	fexp = 6375*(10^(f/50000)-1)
%
%[6] X. Fan and J. Hansen, "Speaker identi?cation within whispered speech audio streams,� Audio,
%Speech, and Language Processing, IEEE Transactions on, vol. 19, pp. 1408 �1421, july 2011.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fs2=4e3;


%c=(1127*log(1+fs2/700))/(10^(fs2/k)-1)
c=10610;
k=50e3;

frq=k*(sign(mel).*log10(abs(mel)/c+1));


if ~nargout
    plot(mel,frq,'-x');
    xlabel(['Frequency (' xticksi 'Mel)']);
    ylabel(['Frequency (' yticksi 'Hz)']);
end
