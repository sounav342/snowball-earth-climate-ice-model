function R = find_R(h,par)
%% finds the partial ice coverage field R of a height array, h or delta h
R = NaN(length(h),1);
idx_full = h >= par.Hcr;
idx_partial = (h >= 0) & (h < par.Hcr);
R(idx_full) = 1;
R(idx_partial) = h(idx_partial) / par.Hcr;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% finds the partial ice coverage field R of a height array, h or delta h
% % Find indices where h is greater than or equal to Hcr
% % Find indices where h is between 0 and Hcr
% % Set R values for total ice cover
% % Set R values for partial or zero ice cover
