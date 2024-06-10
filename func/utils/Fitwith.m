function [FIT3D, b,Residuals]=Fitwith(Volume,Mask,varargin)
% function [FIT3D, b,Residuals]=Fitwith(Volumenii,Mask,varargin)

dim=size(Volume);
size(varargin);
Indices=find(Mask(:,:,:));
[x1,y1,z1]=ind2sub(size(Volume(:,:,:)),Indices);
R=Volume(Indices);


model(1:length(Indices),1)=1;
for k=1:length(varargin)
  model(1:length(Indices),k+1)=varargin{k}(Indices);  
end;

% b=model\R;
 b=pinv(model)*R;

% Fit=model*b;
% clear model
% FIT3D=Volumenii;
% FIT3D.img=zeros(size(Volumenii.img(:,:,:)));
% for pos=1:length(x1)
%    FIT3D.img(x1(pos),y1(pos),z1(pos))=Fit(pos);
% end;
% clear Fit
% FIT3D.img(isnan(FIT3D.img))=0;
% Residuals=FIT3D;
% Residuals.img=(Volumenii.img-FIT3D.img).*Mask.img;
% 
% 






model(1:(prod(size(Mask))),1)=1;
for k=1:length(varargin)
  model(1:(prod(size(Mask))),k+1)=varargin{k}(:);  
end;

Fit=model*b;
clear model
FIT3D=Volume;
% FIT3D.img=zeros(size(Volumenii.img(:,:,:)));
   FIT3D=reshape(Fit,size(Mask));
clear Fit
FIT3D(isnan(FIT3D))=0;
Residuals=FIT3D;
Residuals=(Volume-FIT3D).*Mask;




Residuals=(Volume-FIT3D);

