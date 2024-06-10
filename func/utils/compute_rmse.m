function [ rmse ] = compute_rmse( chi_recon, chi_true )


rmse = 100 * norm( chi_recon(:) - chi_true(:) ) / norm(chi_true(:));


end

