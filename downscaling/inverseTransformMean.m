function EOds_hat = inverseTransformMean(ESgOhat, mu_gO, O_t, sigma2, name_var)
%INVERSETRANSFORMMEAN Inverse transform of expected downscaled values
%
% This function reverses the preprocessing applied during downscaling,
% but specifically for the **expected values** (not stochastic realizations).
%
% Inputs:
%   - ESgOhat : Expected downscaled values in the transformed space 
%               (after PCR regression, anomalies)
%   - mu_gO   : Temporal mean of the transformed original observations
%   - O_t     : Translation constant applied to avoid log(0) in precipitation
%   - sigma2  : Variance of the Gaussian residuals (only used for 'pr')
%   - name_var: Variable name
%                 'tas' → temperature
%                 'pr'  → precipitation
%
% Outputs:
%   - EOds_hat : Expected downscaled values transformed back 
%                to the original scale
%
% Notes:
%     For 'tas' (temperature), the transform was centering, so the inverse
%     simply adds back the mean and subtracts O_t (which is defined 0 anyway).
%     For 'pr' (precipitation), the transform was a log-transform with
%     translation. Because the residuals are Gaussian, the expected value
%     requires adding exp(sigma^2/2).

switch name_var
    case 'tas'
        % Temperature: inverse of centering
        EOds_hat = ESgOhat + mu_gO - O_t;

    case 'pr'
        % Precipitation: inverse of log-transform + lognormal correction
        EOds_hat = exp(ESgOhat + mu_gO + sigma2/2) - O_t;
end

end
