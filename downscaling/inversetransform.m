function Y_hat = inversetransform(SY, mu_Y, Y_t, name_var)
% This function performs the inverse transformation of the downscaled data.
% Depending on the variable name (specified by 'name_var'), it reverses the
% scaling and transformation applied during the preprocessing step.
%
% Inputs:
%   - SY: Downscaled value(s) that need to be transformed back to the original scale
%   - mu_Y: Mean of the original data used for centering
%   - Y_t: Scaling factor used for the original data
%   - name_var: The name of the variable being processed ('tas' or 'pr'), 
%     which determines the specific inverse transformation to apply
%
% Outputs:
%   - Y_hat: The transformed value(s) back to the original scale based on 'name_var'

switch name_var
    case 'tas'
        % If the variable is 'tas' (temperature), apply the inverse transformation
        % Temperature data is centered using mu_Y and scaled by Y_t (which would be 0 anyway), so:
        % Y_hat = Downscaled value (SY) + Mean value (mu_Y) - Scaling factor (Y_t)
        Y_hat = SY + mu_Y - Y_t;
        
    case 'pr'
        % If the variable is 'pr' (precipitation), the data was log-transformed 
        % before scaling. To reverse this, we first take the exponential of SY
        % and then apply the scaling and shifting:
        % Y_hat = exp(Downscaled value (SY) + Mean value (mu_Y)) - Scaling factor (Y_t)
        Y_hat = exp(SY + mu_Y) - Y_t;
        
end

end
