function D = create_dipole_kernel(B0_dir, voxel_size, N, createInKspace)

if createInKspace  
    disp('Dipole kernel created directly in k-space - it does support rotations')
    [ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);

    kx = (kx / max(abs(kx(:)))) / (2*voxel_size(1));
    ky = (ky / max(abs(ky(:)))) / (2*voxel_size(2));
    kz = (kz / max(abs(kz(:)))) / (2*voxel_size(3));

    k2 = kx.^2 + ky.^2 + kz.^2;

%     R_tot = eye(3);
%     R_tot (3,1:3)=B0_dir;
    D = fftshift( 1/3 - (kx * B0_dir(1) + ky * B0_dir(2) + kz * B0_dir(3)).^2 ./ (k2 + eps) );       
else
    disp('Dipole kernel created in the object domain')
    % dipole kernel: (3*cos(theta)^2 - 1) / (4*pi*r^3)

    [Y,X,Z] = meshgrid(-N(2)/2:(N(2)/2-1),...
        -N(1)/2:(N(1)/2-1),...
        -N(3)/2:(N(3)/2-1));

    X = X * voxel_size(1);
    Y = Y * voxel_size(2);
    Z = Z * voxel_size(3);

    d = prod(voxel_size) * ( 3 * (X*B0_dir(1) + Y*B0_dir(2) + Z*B0_dir(3)).^2 - X.^2 - Y.^2 - Z.^2 ) ./ (4 * pi * (X.^2 + Y.^2 + Z.^2).^2.5);

    d(isnan(d)) = 0;
    D = fftn(fftshift(d));
end