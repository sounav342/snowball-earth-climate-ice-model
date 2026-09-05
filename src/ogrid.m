function [h,R,par] = ogrid(ogi,h,R,v_np1,par,S)

for i = 1:length(ogi)
    j = ogi(i);

    % ogi ~ ocean grid index: for terminuses with no adjacent subgrid, ogi gives the
    % index of the adjacent ocean grid cell
    % Note that 'negative' advection from terminus is negated in delta_h
    if j < par.EQ
        delta_hN = par.dt*h(par.Nedge,2)*v_np1(par.Nedge)/(par.R*par.dtheta_rad);
        if delta_hN < 0
            delta_hN=0;
        end
        delta_h = par.dt*S(j) + delta_hN;

        if delta_h > 0
            % Advance of Ice Shelf into Next Ocean Grid Cell
            if delta_h >= par.Hcr
                % grid scale advection
                h(j,3) = delta_h;
                R(j,3) = 1;
            elseif delta_h < par.Hcr
                % subgrid scale advection
                h(j,3) = 0;
                R(j,3) = find_R(delta_h,par);
            end
        elseif delta_h <= 0
            % No change to ocean grid cell ~ no ice advected from terminus
            h(j,3) = 0;
            R(j,3) = 0;
        end

    elseif j > par.EQ
        delta_hS=-par.dt*h(par.Sedge,2)*v_np1(par.Sedge)/(par.R*par.dtheta_rad);
        if delta_hS < 0
            delta_hS = 0;
        end
        delta_h = par.dt*S(j) + delta_hS;

        if delta_h > 0
            % Advance of Ice Shelf into Next Ocean Grid Cell
            if delta_h >= par.Hcr
                % grid scale advection
                h(j,3) = delta_h;
                R(j,3) = 1;
            elseif delta_h < par.Hcr
                % subgrid scale advection
                h(j,3) = 0;
                R(j,3) = find_R(delta_h,par);
            end
        elseif delta_h <= 0
            % No change to ocean grid cell ~ no ice advected from terminus
            h(j,3) = 0;
            R(j,3) = 0;
        end

    elseif j == par.EQ
        if ((par.Nedge==(par.EQ-1)) && (par.Sedge==(par.EQ+1)))
            delta_hN = par.dt*h(par.Nedge,2)*v_np1(par.Nedge)/(par.R*par.dtheta_rad);
            delta_hS = -par.dt*h(par.Sedge,2)*v_np1(par.Sedge)/(par.R*par.dtheta_rad);
            if delta_hN < 0
                delta_hN = 0;
            end
            if delta_hS < 0
                delta_hS = 0;
            end
            delta_h = par.dt*S(j) + delta_hN + delta_hS;

            if delta_h > 0
                if delta_h >= par.Hcr
                    % grid scale advection
                    h(j,3) = delta_h;
                    R(j,3) = 1;
                    % equatorial gap closed
                elseif delta_h < par.Hcr
                    % subgrid scale advection
                    h(j,3) = 0;
                    R(j,3) = find_R(delta_h,par);
                end
            elseif delta_h <= 0
                % No change to ocean grid cell
                h(j,3) = 0;
                R(j,3) = 0;
            end
        end
    end
end
