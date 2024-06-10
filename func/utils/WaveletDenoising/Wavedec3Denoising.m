function outmat=Wavedec3Denoising(inmat,threshold,L,h,method,varargin);
% function outmat=Wavedec3Denoising(inmat,threshold,L,h,method,varargin);
% varargin{1} is a mask
% varargin{2} is noise weighting term


if nargin>6
    if ~isempty(varargin{2})
        [y_r] = wavedec3(inmat.*varargin{2},L,h);
    else
        [y_r] = wavedec3(inmat,L,h);
        
    end
else
    [y_r] = wavedec3(inmat,L,h);
    
end


if nargin>5
    if ~isempty(varargin{1})
        mask1=varargin{1};
        [y_mask] = wavedec3(mask1,L,h);
        [y_all] = wavedec3(ones(size(mask1)),L,h);
        
        for k=1:length(y_mask.dec)
            w_mask.dec{k}=0*y_all.dec{k};
            w_mask.dec{k}(abs(y_all.dec{k})>1.1*abs(y_mask.dec{k}))=0;
            w_mask.dec{k}(abs(y_all.dec{k})<=1.1*abs(y_mask.dec{k}))=1;
            w_mask.dec{k}(abs(y_mask.dec{k})<0.01)=1;
            w_mask.dec{k}=reshape(w_mask.dec{k},size(y_mask.dec{k}));
        end;
    else
        for k=1:length(y_r.dec)
            w_mask.dec{k}=ones(size(y_r.dec{k}));
        end;
    end;
else
    for k=1:length(y_r.dec)
        w_mask.dec{k}=ones(size(y_r.dec{k}));
    end;
    %     indexes_mask=1:prod(size(inmat));
    
end
%  keyboard
% size(w_mask.dec{1});
for k=1:length(y_r.dec)
% size(w_mask.dec{k});
    indexes_dn2.dec{k}=find(and(abs(y_r.dec{k})<3*threshold,w_mask.dec{k}));
    indexes_dn.dec{k} =find(and(abs(y_r.dec{k})<  threshold,w_mask.dec{k}));
    indexes.dec{k}=find(~and(abs(y_r.dec{k})<threshold,w_mask.dec{k}));
end;


% indexes_dn=find((abs(y_r)<threshold));
% indexes=find((abs(y_r)>=threshold));

if strcmp(method,'hard')
    
    for k=1:length(y_r.dec)
        y_r.dec{k}(indexes_dn.dec{k})=0;
    end;
    
    
elseif  strcmp(method,'soft')
    
    
    for k=1:length(y_r.dec)
        y_r.dec{k}(indexes_dn.dec{k})=0;
        % regions that were not truncatedget the phase they had
        % before but a constant amplitude (threshold) subtracted from it
        y_r.dec{k}(indexes.dec{k})=y_r.dec{k}(indexes.dec{k})-exp(1i*angle(y_r.dec{k}(indexes.dec{k})))*threshold;
    end;
    
    
    
    
elseif  strcmp(method,'verysoft')
    
    for k=1:length(y_r.dec)
        y_r.dec{k}(indexes_dn.dec{k})=0;
        y_r.dec{k}(indexes.dec{k})=max(abs(y_r.dec{k}(:)))/(max(abs(y_r.dec{k}(:)))-threshold).*(y_r.dec{k}(indexes.dec{k})-exp(1i*angle(y_r.dec{k}(indexes.dec{k})))*threshold);
    end;
    
    
elseif  strcmp(method,'verysoft2')
    % it only applies the very soft threshold from threshold to 3*threshold
    % and the rest of the values remain unchanged
    for k=1:length(y_r.dec)
                maximum=3*threshold;
        y_r.dec{k}(indexes_dn2.dec{k})=maximum/(maximum-threshold).*(y_r.dec{k}(indexes_dn2.dec{k})-exp(1i*angle(y_r.dec{k}(indexes_dn2.dec{k})))*threshold);
        y_r.dec{k}(indexes_dn.dec{k})=0;
    end;
    
    
    
end;

[outmat] = (waverec3(y_r));
if nargin>6
    if ~isempty(varargin{2})
        outmat= outmat./varargin{2};
    end    
end

