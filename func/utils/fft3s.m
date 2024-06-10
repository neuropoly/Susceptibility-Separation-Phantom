function out = fft3s(data)
% function out = fft3s(data)
% 
% Do 3D FFT on each slice of data and fftshift

out = fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(data,1),2),3),[],1),[],2),[],3),1),2),3); 
