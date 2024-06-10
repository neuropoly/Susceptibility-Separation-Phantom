function out = ifft3s(data)
% function out = ifft3s(data)
% 
% Do 3D FFT on each slice of data and fftshift

out = ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(fftshift(fftshift(fftshift(data,1),2),3),[],1),[],2),[],3),1),2),3); 
