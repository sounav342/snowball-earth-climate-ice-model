function [h,R,par] = sgrid(sgi,h,R,v_np1,par,S,n)

%% Subgrid Parametrization: Advance and Retreat %%
% Note: at any given time step, any additional "melt" beyond the subgrid cell will be
% categorized as heat input into the ocean

for i = 1:length(sgi)
    j = sgi(i);
    if j < par.EQ
        delta_hN = par.dt*h(par.Nedge,2)*v_np1(par.Nedge)/(par.R*par.dtheta_rad);
        if delta_hN < 0
            delta_hN = 0;
        end
        delta_h = par.dt*S(j) + delta_hN;
        % Since 0 < R(j,2) < 1, the volume of subgrid j is
        h_sg = R(j,2)*par.Hcr;
        if delta_h <= 0
            % Retreat of Ice Shelf
            if abs(delta_h) > h_sg
                % subgrid cell j melts entirely
                % additional melt categorized as heat input into the ocean
                h(j,3) = 0; % subgrid melts
                R(j,3) = 0; % subgrid melts
                % apply remainder melt to ocean:

            elseif abs(delta_h) < h_sg
                % subgrid cell j only partially melts
                h(j,3) = 0;
                h_sg = h_sg - abs(delta_h);  % find new subgrid height
                R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
            end
        elseif delta_h > 0
            % Advance of Ice Shelf
            h_sg_2R1 = par.Hcr - h_sg;

            % If advance into next cell, identify the cell type of j+1
            if R(j+1,2)==0
                % ocean cell
                h_adv_2fill_jp1 = par.Hcr;
            elseif R(j+1,2)==1
                % grid cell (isolated shelves)
                h_adv_2fill_jp1 = 0;

                % do we build up ice in j?
            elseif (R(j+1,2)>0) & (R(j+1,2)<1)
                % subgrid cell
                h_adv_2fill_jp1 = par.Hcr - R(j+1,2)*par.Hcr;
                % additional ice mass goes where?
            end

            if delta_h > (h_sg_2R1 + h_adv_2fill_jp1)
                % subgrid cells j & j+1 fill to par.Hcr (limits on advance)
                R(j,3) = 1; % subgrid fills/converts to grid
                h(j,3) = par.Hcr; % subgrid fills/converts to grid
                R(j+1,3) = 1; % next ocean subgrid fills/converts to grid
                h(j+1,3) = par.Hcr; % next ocean subgrid fills/converts to grid
            elseif delta_h > h_sg_2R1
                % subgrid cell j fills entirely and partly fills next subgrid cell j+1
                R(j,3) = 1; % subgrid fills/converts to grid
                h(j,3) = par.Hcr; % subgrid fills/converts to grid
                % find adv height of ice to subgrid cell j+1
                delta_h = delta_h - h_sg_2R1;
                % update R value in subgrid cell j+1`
                R(j+1,3) = R(j+1,2)+find_R(delta_h,par);
            elseif delta_h < h_sg_2R1
                % subgrid cell j partially fills
                h(j,3) = 0;
                h_sg = h_sg + delta_h; % partially fill subgrid
                R(j,3) = find_R(h_sg,par); % update subgrid field
            end
        end
    elseif j > par.EQ
        delta_hS = -par.dt*h(par.Sedge,2)*v_np1(par.Sedge)/(par.R*par.dtheta_rad);
        if delta_hS < 0
            delta_hS = 0;
        end
        delta_h = par.dt*S(j) + delta_hS;
        % Since 0 < R(j,2) < 1, the volume of subgrid j is
        h_sg = R(j,2)*par.Hcr;
        if delta_h <= 0
            % Retreat of Ice Shelf
            if abs(delta_h) > h_sg
                % subgrid cell j melts entirely
                % additional melt categorized as heat input into the ocean
                h(j,3) = 0; % subgrid melts
                R(j,3) = 0; % subgrid melts
                % apply remainder melt to ocean:

                % par.Sedge remains unchanged
            elseif abs(delta_h) < h_sg
                % subgrid cell j only partially melts
                h(j,3) = 0;
                h_sg = h_sg - abs(delta_h);  % find new subgrid height
                R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
            end
        elseif delta_h > 0
            % Advance of Ice Shelf
            h_sg_2R1 = par.Hcr - h_sg;

            % If advance into next cell, identify the cell type of j-1
            if R(j-1,2)==0
                % ocean cell
                h_adv_2fill_jp1 = par.Hcr;
            elseif R(j-1,2)==1
                % grid cell (isolated shelves)
                h_adv_2fill_jp1 = 0;

                % do we build up ice in j?
            elseif (R(j-1,2)>0) & (R(j-1,2)<1)
                % subgrid cell
                h_adv_2fill_jp1 = par.Hcr - R(j-1,2)*par.Hcr;
                % additional ice mass goes where?
            end

            if delta_h > (h_sg_2R1 + h_adv_2fill_jp1)
                % subgrid cells j & j+1 fill to par.Hcr (limits on advance)
                R(j,3) = 1; % subgrid fills/converts to grid
                h(j,3) = par.Hcr; % subgrid fills/converts to grid
                R(j-1,3) = 1; % next ocean subgrid fills/converts to grid
                h(j-1,3) = par.Hcr; % next ocean subgrid fills/converts to grid
            elseif delta_h > h_sg_2R1
                % subgrid cell j fills entirely and partly fills next subgrid cell j-1
                R(j,3) = 1; % subgrid fills/converts to grid
                h(j,3) = par.Hcr; % subgrid fills/converts to grid
                % find adv height of ice in ocean subgrid cell j-1
                delta_h = delta_h - h_sg_2R1;
                % update value in ocean subgrid cell j+1
                R(j-1,3) = R(j-1,2)+find_R(delta_h,par);
            elseif delta_h < h_sg_2R1
                % subgrid cell j partially fills
                h(j,3) = 0;
                h_sg = h_sg + delta_h; % partially fill subgrid
                R(j,3) = find_R(h_sg,par); % update subgrid field
            end
        end
    elseif j==par.EQ
        %% Note: symmetric melt/advance might be problematic
        if ((par.Nedge==(par.EQ-1)) && (par.Sedge==(par.EQ+1)))
            %% equatorial mass balance involves both hemispheres
            delta_hN = par.dt*h(par.Nedge,2)*v_np1(par.Nedge)/(par.R*par.dtheta_rad);
            delta_hS = -par.dt*h(par.Sedge,2)*v_np1(par.Sedge)/(par.R*par.dtheta_rad);
            if delta_hN < 0
                delta_hN = 0;
            end
            if delta_hS < 0
                delta_hS = 0;
            end
            delta_h = par.dt*S(j) + delta_hN + delta_hS;
            % Since 0 < R(j,2) < 1, the volume of subgrid j is
            h_sg = R(j,2)*par.Hcr;

            if delta_h <= 0
                %% Retreat of Ice Shelf
                if abs(delta_h) > h_sg
                    % subgrid cell j melts entirely and
                    % additional melt categorized as heat input into the ocean
                    h(j,3) = 0; % subgrid melts
                    R(j,3) = 0; % subgrid melts
                    % apply remainder melt to ocean:

                    % unchanged terminuses
                elseif abs(delta_h) < h_sg
                    % subgrid cell j only partially melts
                    h(j,3) = 0;
                    h_sg = h_sg - abs(delta_h); % find new subgrid height
                    R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
                end
            elseif delta_h > 0
                % Advance of Ice Shelf in regime of subgrid @ equator
                % symmetric terminuses @ EQ-1 on either side
                h_sg_2R1 = par.Hcr - h_sg;
                if delta_h > h_sg_2R1
                    % subgrid fills entirely, but growth limited to Hcr
                    h(j,3) = par.Hcr; % subgrid fills/converts to grid
                    R(j,3) = 1; % subgrid fills/converts to grid
                    % transition from partial to global
                elseif delta_h < h_sg_2R1
                    % subgrid cell j partially fills
                    h(j,3) = 0;
                    h_sg = h_sg + delta_h; % partially fill subgrid
                    R(j,3) = find_R(h_sg,par); % update subgrid field
                end
            end

        elseif ((par.Nedge==(par.EQ-1)) && (par.Sedge~=(par.EQ+1)))
            %% equatorial mass balance involves NH
            if par.Sedge==0
                fprintf(2,'***error: par.Sedge cannot equal zero\n');
                error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
            end
            delta_hN = par.dt*h(par.Nedge,2)*v_np1(par.Nedge)/(par.R*par.dtheta_rad);
            if delta_hN < 0
                delta_hN = 0;
            end
            delta_h = par.dt*S(j) + delta_hN;
            % Since 0 < R(j,2) < 1, the volume of subgrid j is
            h_sg = R(j,2)*par.Hcr;

            if delta_h <= 0
                % Retreat of Ice Shelf
                if abs(delta_h) > h_sg
                    % subgrid cell j melts entirely and
                    % additional melt cateogrized as heat input into the ocean
                    h(j,3) = 0; % subgrid melts
                    R(j,3) = 0; % subgrid melts
                    % apply remainder melt to ocean:

                    % par.Nedge remains unchanged
                elseif abs(delta_h) < h_sg
                    % subgrid cell j only partially melts
                    h(j,3) = 0;
                    h_sg = h_sg - abs(delta_h);  % find new subgrid height
                    R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
                end
            elseif delta_h > 0
                % Equatorial Advance of Ice Shelf from NH in subgrid scale
                h_sg_2R1 = par.Hcr - h_sg;
                if delta_h > h_sg_2R1
                    % subgrid fills entirely, but growth limited to Hcr
                    h(j,3) = par.Hcr; % subgrid fills/converts to grid
                    R(j,3) = 1; % subgrid fills/converts to grid
                    % transition from partial to global hemisphere
                elseif delta_h < h_sg_2R1
                    % subgrid cell j partially fills
                    h(j,3) = 0;
                    h_sg = h_sg + delta_h; % partially fill subgrid
                    R(j,3) = find_R(h_sg,par); % update subgrid field
                end
            end
        elseif((par.Nedge~=(par.EQ-1)) && (par.Sedge==(par.EQ+1)))
            %% equatorial mass balance involves SH
            if par.Nedge==0
                fprintf(2,'***error:par.Nedge cannot equal zero');
                error('Snowball:InvariantViolation','Model invariant failed; inspect the preceding diagnostic.');
            end
            delta_hS = -par.dt*h(par.Sedge,2)*v_np1(par.Sedge)/(par.R*par.dtheta_rad);
            if delta_hS < 0
                delta_hS = 0;
            end
            delta_h = par.dt*S(j) + delta_hS;
            h_sg = R(j,2)*par.Hcr;

            if delta_h <= 0
                % Retreat of Ice Shelf
                if abs(delta_h) > h_sg
                    % subgrid cell j melts entirely and
                    % additional melt categorized as heat input into the ocean
                    h(j,3) = 0; % subgrid melts
                    R(j,3) = 0; % subgrid melts
                    % apply remainder melt to the ocean:

                    % par.Sedge remains unchanged
                elseif abs(delta_h) < h_sg
                    % subgrid cell j only partially melts
                    h(j,3) = 0;
                    h_sg = h_sg - abs(delta_h);  % find new subgrid height
                    R(j,3) = find_R(h_sg,par); % update R in subgrid cell j
                end
            elseif delta_h > 0
                % Equatorial Advance of Ice Shelf from SH in subgrid scale
                h_sg_2R1 = par.Hcr - h_sg;
                if delta_h > h_sg_2R1
                    % subgrid fills entirely, but growth limited to Hcr
                    h(j,3) = par.Hcr; % subgrid fills/converts to grid
                    R(j,3) = 1; % subgrid fills/converts to grid
                    % transition from partial to global hemisphere
                elseif delta_h < h_sg_2R1
                    % subgrid cell j partially fills
                    h(j,3) = 0;
                    h_sg = h_sg + delta_h; % partially fill subgrid
                    R(j,3) = find_R(h_sg,par); % update subgrid field
                end
            end
        end
    end
end
