function [alpha_surface,alpha_i_local,alpha_o_local,R] = surface_albedo_fractional(R,par)
% Grid-cell surface albedo for fractional sea-ice cover.

R = validate_fractional_ice_cover(R,par,'surface_albedo_fractional');
alpha_i_local = expand_albedo(par.alpha_i,par,'alpha_i');
alpha_o_local = expand_albedo(par.alpha_o,par,'alpha_o');

% Area-weight ice and open-ocean albedos over the fractional grid cell,
% following the Pollard and Kasting (2005) fractional-cell treatment.
alpha_surface = R .* alpha_i_local + (1 - R) .* alpha_o_local;
end

function alpha = expand_albedo(alphaIn,par,name)
if ~isnumeric(alphaIn) || isempty(alphaIn) || any(~isfinite(alphaIn(:)))
    error('surface_albedo_fractional:invalidAlbedo', ...
        'par.%s must be finite numeric albedo values.',name);
end

if isscalar(alphaIn)
    alpha = alphaIn .* ones(par.nj,1);
else
    alpha = alphaIn(:);
    if numel(alpha) ~= par.nj
        error('surface_albedo_fractional:albedoDimensionMismatch', ...
            'par.%s must be scalar or contain par.nj=%d elements, got %d.', ...
            name,par.nj,numel(alpha));
    end
end
end
