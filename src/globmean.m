% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Take the global area-weighted mean
%     % Ensure f is a numerical array
%     % Ensure par.c is defined and is a numerical array

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Take the global area-weighted mean

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Take the global area-weighted mean

function gm = globmean(f, par)
    f = f(:);
    if isfield(par,'dc')
        dc = par.dc(:);
    else
        dc = diff(par.c(:));
    end
    gm = sum((f(2:end) + f(1:end-1)) .* dc / 2) / sum(dc);
end

