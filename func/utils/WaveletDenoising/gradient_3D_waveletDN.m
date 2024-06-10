function [dx, dy, dz] = gradient_3D_waveletDN(I,varargin)
%[dx, dy, dz] = gradient_3D_waveletDN(I,varargin)
%[dx, dy, dz] = gradient_3D_waveletDN(I,res,method,noiselevel,L,wname,)
% varargin{1} is the resolution 
% varargin{2} corresponds to the gradientmethod, the default is 0 which corresponds
% to the mid gradient definition
% varargin{3} is the noise level 
% varargin{4} is number of decomposition levels 
% varargin{5} is the wavelet name  
% varargin{6} is the mask  


if nargin==1
    res=[1,1,1];
else
    if isempty(varargin{1})
    res=[1,1,1];
else
    res=varargin{1};
    end;
end

if nargin>=3
    gradientmethod = varargin{2};
else
    gradientmethod = 0;
end;

if nargin>=4
    threshold = varargin{3};
else
    threshold = 1;
end;

if nargin>=5
    L= varargin{4};
else
    L = 4;
end;
if nargin>=6
    h= varargin{5};
else
    h = 'sym4';
    h = 'bior6.8';
end;
if nargin>=7
    mask= varargin{6};
else
    mask=[];
end;

outmat=Wavedec3Denoising(I,threshold,L,h,mask);
[dx, dy, dz] = gradient_3D(outmat,res,method);

dthreshold=threshold*(1./res)*sqrt(2);
dx=Wavedec3Denoising(dx,dthreshold,L,h,mask);
dy=Wavedec3Denoising(dy,dthreshold,L,h,mask);
dz=Wavedec3Denoising(dz,dthreshold,L,h,mask);
 UNFINISHED UNTESTED CODE