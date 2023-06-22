function [P1,f] = fncFFT(X,Fs)
%fourier
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%if mode ==1
    f = Fs*(0:(L/2))/L;
%     plot(f,P1/max(P1));
    
    %xlim([0,75])
%end
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
end