
% Calculate correlation coefficient matrix for matrix X
function corrcoef_PC = corrcoef_mat(X)
    [coef, score, latent] = princomp(X);
    corrcoef_PC = (coef)*diag(latent)*inv(coef);
end