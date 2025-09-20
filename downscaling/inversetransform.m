function Ods_hat = inversetransform(SgOhat, mu_gO, O_t, name_var)
%INVERSETRANSFORM Inverse transformation of downscaled data
%
% This function reverses the preprocessing transformations applied during
% downscaling for temperature ('tas') or precipitation ('pr').
%
% Inputs:
%   - SgOhat  : Downscaled anomaly values (after PCR, on transformed data)
%   - mu_gO   : Temporal mean of the transformed original observations used for centering
%   - O_t     : Translation term applied to avoid log(0) for precipitation
%   - name_var: Variable name ('tas' for temperature, 'pr' for precipitation)
%   
% Outputs:
%   - Ods_hat : Downscaled data transformed back to the original scale

switch name_var
    case 'tas'
        % Temperature ('tas') was centered by subtracting mu_gO.
        % Inverse transform simply adds back the mean and subtracts O_t (usually 0).
        Ods_hat = SgOhat + mu_gO - O_t;
        
    case 'pr'
        % Precipitation ('pr') was log-transformed with a translation O_t to avoid zeros.
        % Inverse transform exponentiates the anomaly, adds the mean, and subtracts the translation.
        Ods_hat = exp(SgOhat + mu_gO) - O_t;
end

end
