function S=smooth3D(X,FWHM,vox)

% will apply a 3D gaussian smoothing filter to a 3D image, call with
% S=smooth(X,FWFHM,voxdims)
% FWHM is in mm

%zero fill to required size
dimsorig=size(X);
if length(dimsorig)>=3
    vol=1
else
    vol=0;
    if dimsorig(1)==1
        X=repmat(X,[8 1 8]);
    end;
    if dimsorig(2)==1
        X=repmat(X,[1 8 8]);
    end;    
end;


dim=size(X);
%create gaussian smoothing matrix

x1=((-dim(1)/2:dim(1)/2-1))*vox(1);
x2=((-dim(2)/2:dim(2)/2-1))*vox(2);
x3=((-dim(3)/2:dim(3)/2-1))*vox(3);


[X1,X2,X3]=ndgrid(x1,x2,x3);

dev=(FWHM/2.35);

GAUSS=exp(-(X1.^2+X2.^2+X3.^2)/(2*dev^2));
  GAUSS=GAUSS/(sum(abs(GAUSS(:))));


S=zeros(dim);

S=(ifft3s(bsxfun(@times,fft3s(X),fft3s(GAUSS))));


if length(dimsorig)<3
   if dimsorig(1)==1
        S=S(end/2,:,end/2);
    end;
    if dimsorig(2)==1
        S=S(:,end/2,end/2);
    end;    

end





