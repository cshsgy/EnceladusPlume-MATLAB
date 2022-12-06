function main_func(folder_name,file_name)
% /home/shc/Enceladus/2021_25Jul_ParallelCrackDyn/2021_Dec14_SmallWidthRatio/
%% (Hyper-)Parameters
% Naming convention: the first alphabet for an input parameter is
% capicalized...
% Crack dynamics input file (consider running it too?)
Inp_file = strcat(folder_name,file_name);
% Heat conduction discretization
Dx = 0.05;
% Heat conduction solve up to. Note: characteristic length~0.34m
Max_x = 0.4;
% Number of Vertical Levels 
Nz = 50;
% Effective surface temperature
Te = 68;
% Periods to run. sqrt(Kappa*T)~3m here, good enough to reach eq
Periods_Run = 300;

% physical parameters
Kappa = 1e-6; % Diffusivity
Period = 118800;
Lv = 2.84E6; % Sublimation latent heat
K = 3; % Conductivity
G = 0.113; % Gravity
Cd = 0.002; % Drag Coefficient

% Derived Parameters
load(Inp_file,'t_rec','h_rec','w_rec','L');
hmax = max(h_rec);
[t_rec,h_rec,w_rec] = SimplifyCrackDynamicsData(t_rec,Period,h_rec,w_rec);
D = L/10; % Crack depth taken as 10% of the full length
if max(h_rec)>D
	disp(strcat(file_name,' encounters water to the surface, PROGRAM TERMINATED WITHOUT SAVING.'));
	return;
end

Dt = 0.01/Kappa*Dx^2;

%% Variables initialization
% Solve up to 1 grid from the surface to protect the heat conduction val
z_wet = linspace(min(h_rec),D,Nz+1);
z_wet = z_wet(1:end-1);
dz = (max(z_wet)-min(z_wet))/(Nz-1);

% Solve at most 1.5m into the wall
xs = 0.0:Dx:Max_x;
nx = length(xs);
% Initialize to close to water temp
T = ones(Nz,nx)*273.15;
% Not wetted walls initialized to be colder
T(z_wet>max(h_rec),:) = T(z_wet>max(h_rec),:)-5;
for i=1:Nz
    Ts = SurfaceTemperature(Te,T(i,1),D-z_wet(i)+max(xs),K);
    if Ts<0
        error('NUMERICAL ERROR: No Solution to the surface temperature.');
    end
    T(i,:) = T(i,:)-xs*(T(i,1)-Ts)*2/pi/(D-z_wet(i)+max(xs));
end

% Recording variables
phi_rec = [];
r_rec = [0.6]; % just some random number put up here
time_rec = [];

Ev_rec = {};
Tw_rec = {};

% For debugging purposes
T_rec = {};
width_rec = [];
zs_rec = {};

%% Run Main Solver 
t = 0; iter = 1; tic;
while(t<Periods_Run*Period)
    time_rec(end+1) = t;
    if mod(iter,50)==0
        % Printing diagnostics every several steps
        disp('Number of iteration:');
        disp(iter);
        tmp = toc;
		fprintf('Estimated time left (h): %.2f\n',t/(Periods_Run*Period));
    end
    % Interpolate to the above-liquid grid
    height = interp1(t_rec,h_rec,mod(t,Period));
    depth = D-height;
    width = interp1(t_rec,w_rec,mod(t,Period));
    zs = linspace(0,depth,Nz);
    T2_o = T(:,2);
    T2 = interp1(z_wet,T2_o,zs+height,'linear','extrap');
    
    % Calculate the gas
    [phi_rec(end+1),Tw,Ev,r_rec(end+1)] =
    GasDynamicsMarchInTimeConstrained(width,zs,T2,K,Lv,G,Cd,Dx,r_rec(end),hmax);
    T_rec{end+1} = T;
    width_rec(end+1) = width;
    zs_rec{end+1} = zs;
    
    % Interpolate back
    Tw = interp1(zs+height,Tw,z_wet,'linear','extrap');
    Ev = interp1(zs+height,Ev,z_wet,'linear','extrap');
    Tw(z_wet<height) = 273.15;
    Ev(z_wet<height) = 0;
    
    % Take record
    Ev_rec{end+1} = Ev;
    Tw_rec{end+1} = Tw;
    
    % Calculate heat diffusion, parallelized
    for ii=1:Nz
        T(ii,:) = HeatDiffusionMarchInTime(z_wet(ii),xs,height,depth,...
            T(ii,:),Tw(ii),Kappa,K,Te,Dx,Dt);
    end
    
    t = t+Dt; iter = iter+1;
end
save(strcat('./run_results/',file_name));
end
