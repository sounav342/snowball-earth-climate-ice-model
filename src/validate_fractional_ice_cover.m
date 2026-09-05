function R = validate_fractional_ice_cover(R,par,callerName)
% Validate and orient the existing fractional ice-cover field.

if nargin < 3
    callerName = 'validate_fractional_ice_cover';
end

if ~isnumeric(R) || isempty(R)
    error('validate_fractional_ice_cover:invalidType', ...
        '%s expected numeric, non-empty fractional ice cover R.',callerName);
end

R = R(:);
if ~isfield(par,'nj') || numel(R) ~= par.nj
    error('validate_fractional_ice_cover:dimensionMismatch', ...
        '%s expected R to contain par.nj=%d elements, got %d.', ...
        callerName,par.nj,numel(R));
end

if any(~isfinite(R))
    error('validate_fractional_ice_cover:nonFinite', ...
        '%s received non-finite values in R.',callerName);
end

roundoffTol = 1.0e-12;
minR = min(R);
maxR = max(R);
if minR < -roundoffTol || maxR > 1 + roundoffTol
    error('validate_fractional_ice_cover:outOfRange', ...
        '%s received R outside [0,1]: min(R)=%g, max(R)=%g.', ...
        callerName,minR,maxR);
end

if minR < 0 || maxR > 1
    R = min(1,max(0,R));
end
end
