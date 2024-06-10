function outmat=WaveletDenoising3D(inmat,threshold,L,h,method,varargin);

 [y_r,L3D] = mdwt_complex_3D(inmat,h,L);
%  [y_r] = wavedec3(inmat,L,h);


if nargin>5
    if ~isempty(varargin{1})
        mask1=varargin{1};
        [y_mask,L2D] = mdwt_complex_3D(mask1,h,L);
        [y_all,L2D] = mdwt_complex_3D(ones(size(mask1)),h,L);
%         [y_mask,L3D] = wavedec3(mask1,L,h);
%         [y_all,L3D] = wavedec3(ones(size(mask1)),L,h);
%         clear y_r3
%         subplot(121)
%         hold off
%         plot(y_r)
%         hold on
%         plot(y_r2)
%         plot(y_r2,'r')
        w_mask(abs(y_all)>1.1*abs(y_mask))=0;
        w_mask(abs(y_all)<=1.1*abs(y_mask))=1;
        w_mask(abs(y_mask)<0.01)=1;
        w_mask=reshape(w_mask,size(y_r));
%         plot(w_mask,'g')
%         [outmat] = waverec2(w_mask,L2D,h);
%         subplot(122)
%         imagesc(abs(outmat))
        
%         indexes_mask=find(abs(y_r)<threshold)
    else
     w_mask=ones(size(y_r));
       
    end;
else
    w_mask=ones(size(y_r));
%     indexes_mask=1:prod(size(inmat));

end
% keyboard
 indexes_dn=find(and(abs(y_r)<threshold,w_mask));
indexes=find(~and(abs(y_r)<threshold,w_mask));


% indexes_dn=find((abs(y_r)<threshold));
% indexes=find((abs(y_r)>=threshold));

if strcmp(method,'hard')
    
    y_r(indexes_dn)=0;
    
    
elseif  strcmp(method,'soft')
    
    
    y_r(indexes_dn)=0;
    y_r(indexes)=y_r(indexes)-exp(1i*angle(y_r(indexes)))*threshold;
    
    
    
    
elseif  strcmp(method,'verysoft')
    
    y_r(indexes_dn)=0;
    y_r(indexes)=max(abs(y_r(:)))/(max(abs(y_r(:)))-threshold).*(y_r(indexes)-exp(1i*angle(y_r(indexes)))*threshold);
    
    
    
end;

%   [outmat] = conj(waverec3(y_r,L3D,h));
 [outmat] = conj(midwt_complex_3D(y_r,h,L3D));
