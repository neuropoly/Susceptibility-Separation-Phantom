function [dx, dy, dz] = gradient_3D(I,varargin)
% varargin{2} corresponds to the gradientmethod, the default is 0 which corresponds
% to the mid gradient definition
%defined based on div_op by Gilles Puy
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


gradientmethod;

switch gradientmethod
    case 1
%defined based on gradient_op by Gilles Puy - it is a forward differential
dx = cat(1, I(2:end,:,:)-I(1:end-1,:,:) , zeros(1, size(I, 2), size(I, 3))) ;
dy = cat(2, I(:,2:end,:)-I(:,1:end-1,:) , zeros(size(I, 1), 1, size(I, 3))) ;
dz = cat(3, I(:,:,2:end)-I(:,:,1:end-1) , zeros(size(I, 1), size(I, 2), 1));
    case 0
%defined based on my favorite way - midpoint differenttial
dx = cat(1, I(2,:,:)-I(1,:,:) , 0.5*(I(3:end,:,:)-I(1:end-2,:,:)) ,I(end,:,:)-I(end-1,:,:)) ;
dy = cat(2, I(:,2,:)-I(:,1,:) , 0.5*(I(:,3:end,:)-I(:,1:end-2,:)) ,I(:,end,:)-I(:,end-1,:)) ;
dz = cat(3, I(:,:,2)-I(:,:,1) , 0.5*(I(:,:,3:end)-I(:,:,1:end-2)) ,I(:,:,end)-I(:,:,end-1));
    case -1
%defined based on gradient_op by Gilles Puy - it is a backwards differential
dx = cat(1, zeros(1, size(I, 2), size(I, 3)), I(2:end,:,:)-I(1:end-1,:,:) ) ;
dy = cat(2, zeros(size(I, 1), 1, size(I, 3)), I(:,2:end,:)-I(:,1:end-1,:) ) ;
dz = cat(3, zeros(size(I, 1), size(I, 2), 1), I(:,:,2:end)-I(:,:,1:end-1) );
    case 2
%          keyboard
%defined based on my favorite way - but very large midpoint differenttial
dx = cat(1, I(2:3,:,:)-I(1:2,:,:) , 0.25*(I(5:end,:,:)-I(1:end-4,:,:)) ,I((end-1):(end),:,:)-I((end-2):(end-1),:,:)) ;
dy = cat(2, I(:,2:3,:)-I(:,1:2,:) , 0.25*(I(:,5:end,:)-I(:,1:end-4,:)) ,I(:,(end-1):(end),:)-I(:,(end-2):(end-1),:)) ;
dz = cat(3, I(:,:,2:3)-I(:,:,1:2) , 0.25*(I(:,:,5:end)-I(:,:,1:end-4)) ,I(:,:,(end-1):(end))-I(:,:,(end-2):(end-1)));

    otherwise
        disp('gradient is not performed');
end
    dx = dx /res(1);
    dy = dy /res(2);
    dz = dz /res(3);

end
% %testing code
% [X,Y,Z]=meshgrid([-5:1:5],[-5:1:5],[-5:1:5]);
% [gy, gx, gz]=gradient(X,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(X,res(1),res(2)*10,res(3));
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))
% [gy, gx, gz]=gradient(Y,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(Y,res(1),res(2)*10,res(3));
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))
% [gy, gx, gz]=gradient(Z,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(Z,res(1),res(2)*10,res(3));
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))
% 
% [X,Y,Z]=meshgrid([-5:1:5],[-5:1:5],[-5:1:5]);
% [gy, gx, gz]=gradient(X,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(X,res(1),res(2)*10,res(3),-1);
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))
% [gy, gx, gz]=gradient(Y,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(Y,res(1),res(2)*10,res(3),-1);
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))
% [gy, gx, gz]=gradient(Z,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(Z,res(1),res(2)*10,res(3),-1);
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))
% 
% [X,Y,Z]=meshgrid([-5:1:5],[-5:1:5],[-5:1:5]);
% [gy, gx, gz]=gradient(X,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(X,res(1),res(2)*10,res(3),1);
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))
% [gy, gx, gz]=gradient(Y,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(Y,res(1),res(2)*10,res(3),1);
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))
% [gy, gx, gz]=gradient(Z,res(2)*10,res(1),res(3));
% [dx, dy, dz]=gradient_3D(Z,res(1),res(2)*10,res(3),1);
% image_view3(cat(4,gx,dx,gy,dy,gz,dz))

