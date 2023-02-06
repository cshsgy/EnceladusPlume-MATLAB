Tb = 273.15; % Bottom temperature in K
depth = 1500; % in meters
width = 0.075; % in meters

% For stability here we use dichotomy
r_l = 1e-5;
r_r = 1-1e-5;
iter = 0;
while(abs(r_l-r_r)>1e-4)
    iter = iter+1;
    display(iter);
    r_m = (r_l+r_r)/2;
    [phi_to_zero,rho_to_zero,mach_top] = solve_function(Tb,depth,width,r_m);
    if(mach_top==0) % has not reached top
        if(rho_to_zero)<(phi_to_zero)
            % rho becomes zero first.
            % rho_0 increases with r, but phi_0 decreases with r
            % So need to increase r
            r_l = r_m;
        else
            r_r = r_m;
        end
    else
        if(mach_top>1)
            % increase r, increase rho, decrease v
            r_l = r_m;
        else
            % decrease r, decrease rho, increase v
            r_r = r_m;
        end
    end
end

r = (r_l+r_r)/2;
display(r);