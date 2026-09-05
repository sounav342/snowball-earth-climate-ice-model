function Ts = fill_ice_surface_temperature(Ts,R,par)
% Fill missing ice-surface temperatures over ice-covered cells.
% Newly glaciated cells can inherit NaN Ts from previously open ocean cells.
% Use the nearest finite ice-surface temperature, preferring the poleward
% direction in each hemisphere, so full-snowball handoffs remain finite.

Ts = Ts(:);
R = R(:);

missing_ice = find(R>0 & isnan(Ts));
if isempty(missing_ice)
    return;
end

valid_ice = find(R>0 & ~isnan(Ts));
if isempty(valid_ice)
    Ts(missing_ice) = par.T_f;
    return;
end

for k = 1:numel(missing_ice)
    j = missing_ice(k);

    if j<par.EQ
        candidates = valid_ice(valid_ice<j);
    elseif j>par.EQ
        candidates = valid_ice(valid_ice>j);
    else
        candidates = valid_ice;
    end

    if isempty(candidates)
        candidates = valid_ice;
    end

    [~,idx] = min(abs(candidates-j));
    Ts(j) = Ts(candidates(idx));
    valid_ice = sort([valid_ice; j]);
end

end
