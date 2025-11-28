% CP2_RF_1005_2025
% 1-D FDTD of a three-region slab with three different biological tissues
% Region 2 (tissue) is between interfaces at Hy = 500.5 and 600.5
% => Ex indices 501:600, Hy indices 500:600

clear; clc; close all;

%% Physical constants
c0   = 299792458;
mu0  = 4*pi*1e-7;              % H/m
eps0 = 1/(mu0*c0^2);           % F/m

%% Grid parameters
Nz  = 1000;                    % number of Ex nodes (1..1000)
dz  = 1.5e-2;                  % spatial step [m]
dt  = 2.5e-11;                 % time step [s]
Nt  = 3000;                    % total time steps (~75 ns, enough for 22.5–62.5 ns window)

% Time arrays
tE = ((0:Nt-1)+0.5)*dt;        % Ex at half time steps
tH = (0:Nt-1)*dt;              % Hy at integer steps (not strictly needed)

%% Gaussian source at Ex(k=1)
Tg = 1e-9;                     % pulse width parameter T = 1 ns
gaussian = @(tt) (tt >= 0 & tt <= 6*Tg) .* exp(-((tt - 3*Tg).^2)/(Tg^2)); %create anonymous gaussian fcn with input tt

%% Tissue properties (Region 2)
% 1) Bone marrow @ 64MHz (lossless for assignment)
tissue(1).name  = 'Bone marrow';
tissue(1).eps_r = 7.2103;
tissue(1).sigma = 0;

% 2) Fat @ 128MHz
tissue(2).name  = 'Fat';
tissue(2).eps_r = 5.9215;
tissue(2).sigma = 0.03687;     % S/m

% 3) Bladder @ 300MHz
tissue(3).name  = 'Bladder';
tissue(3).eps_r = 20.093;
tissue(3).sigma = 0.31684;     % S/m

%% Probe locations (Ex indices)
kProbe = [480, 501, 560, 620];      % as specified

% Time window for plotting: 22.5 ns to 62.5 ns
tStart = 22.5e-9;
tEnd   = 62.5e-9;
nStart = round(tStart / dt);        % should be ~900
nEnd   = round(tEnd   / dt);        % should be ~2500
nRange = nStart:nEnd;
nKeep  = numel(nRange);

% Storage: Ex_hist(case, probeIndex, timeIndex) for plotting only
Ex_hist = zeros(3, numel(kProbe), nKeep);

% For Part 3 (bone marrow only) we store the FULL time history at k=560,620
Ex560_full = zeros(1, Nt);          % Ex at k=560, case 1
Ex620_full = zeros(1, Nt);          % Ex at k=620, case 1

%% ---------------- MAIN LOOP OVER TISSUE CASES ----------------
for icase = 1:3
    
    %---------------- Material profiles at Ex nodes ----------------
    eps_z   = eps0 * ones(1, Nz);       % default: free space
    sigma_z = zeros(1, Nz);             % default: lossless
    
    % Region 2 (tissue): Ex indices 501:600 (between Hy 500.5 and 600.5)
    eps_z(501:600)   = tissue(icase).eps_r * eps0;
    sigma_z(501:600) = tissue(icase).sigma;
    
    % Precompute Ex constant coefficients (from derived equation 4)
    % E^{n+0.5}_k = Ca(k) * E^{n-0.5}_k - Cb(k) * (Hy_k^n - Hy_{k-1}^n)
    Ca = zeros(1, Nz);      % (1- ...) term
    Cb = zeros(1, Nz);      %
    for k = 1:Nz
        alpha = dt * sigma_z(k) / (2 * eps_z(k));  % α_k = dt σ / (2 ε)
        Ca(k) = (1 - alpha) / (1 + alpha);
        Cb(k) = (dt / (eps_z(k) * dz)) / (1 + alpha);
    end
    
    %---------------- Initialize fields ----------------
    Ex_old = zeros(1, Nz);      % E^{(n-0.5)}
    Ex_new = zeros(1, Nz);      % E^{(n+0.5)}
    Hy     = zeros(1, Nz-1);    % H_y(i) at (i+0.5)Δz, i=1..Nz-1
    
    keepIdx = 1;                % index for Ex_hist storage
    
    %---------------- Time stepping ----------------
    for n = 0:Nt-1
        
        t_half = (n + 0.5) * dt;    % time for E^{(n+0.5)}
        
        % ---- 1) Update Ex everywhere except boundaries ----
        for k = 2:Nz-1
            Ex_new(k) = Ca(k) * Ex_old(k) ...
                      - Cb(k) * (Hy(k) - Hy(k-1));
            % sparate the derived equation intwo two terms with a common denominator
            % allows to rearrange such that:
            % Ca = (1-...) / (1+...)            (precomputed)
            % Cb = -dt/(eps * dz) / (1+...)     (precomputed)
        end
        
        % Left boundary: poll Gaussian source at k=1 @ t_half
        Ex_new(1) = gaussian(t_half);
        
        % Right boundary: PEC at k=Nz
        Ex_new(Nz) = 0.0;
        
        % ---- 2) Update Hy from n to n+1 (derived equation 3) ----
        % Hy_i^{n+1} = Hy_i^n - dt/(mu0*dz) * (Ex_{i+1}^{n+0.5} - Ex_i^{n+0.5})
        for i = 1:Nz-1
            Hy(i) = Hy(i) ...
                  - (dt / (mu0 * dz)) * (Ex_new(i+1) - Ex_new(i));
        end
        
        % ---- 3) Store Ex at probes for plotting window ----
        if n >= nStart && n <= nEnd
            for p = 1:numel(kProbe)
                Ex_hist(icase, p, keepIdx) = Ex_new(kProbe(p));
            end
            keepIdx = keepIdx + 1;
        end
        
        % ---- 4) For bone marrow case, store full history at 560 & 620 ----
        if icase == 1
            Ex560_full(n+1) = Ex_new(560);
            Ex620_full(n+1) = Ex_new(620);
        end
        
        % Advance fields in time
        Ex_old = Ex_new;
    end
    
    fprintf('Finished tissue case %d: %s\n', icase, tissue(icase).name);
end

%% ---------------- PLOTTING: 12 CURVES ----------------
t_plot    = tE(nRange);         % time vector for stored window (s)
t_plot_ns = t_plot * 1e9;       % in ns

styles = {'k-','r--','b:'};     % bone, fat, bladder

for p = 1:numel(kProbe)
    figure;
    hold on;
    for icase = 1:3
        plot(t_plot_ns, squeeze(Ex_hist(icase,p,:)), styles{icase}, 'LineWidth', 1.3);
    end
    hold off;
    grid on;
    xlabel('Time [ns]');
    ylabel(sprintf('E_x at k = %d', kProbe(p)));
    legend({tissue(1).name, tissue(2).name, tissue(3).name}, 'Location', 'Best');
    title(sprintf('Electric field at index k = %d', kProbe(p)));
end

%% ---------------- PART 3: Bone marrow reflection/transmission ----------------
% Use FULL time history at k=560 and k=620 for bone marrow (case 1)

bm = 1;  % index of bone marrow case in Ex_hist

% Windowed time vector (same one used for plotting)
t_win = t_plot;          % seconds, length = nKeep

% Ex at k = 560 and 620 for bone marrow over 22.5–62.5 ns window
Ex560_win = squeeze(Ex_hist(bm, 3, :));   % k = 560
Ex620_win = squeeze(Ex_hist(bm, 4, :));   % k = 620

% Choose a split time between the first and second pulses.
% Midpoint of the window is simple and works with your plots:
tSplit = 0.5*(22.5e-9 + 62.5e-9);   % 42.5 ns

% Early and late parts of the signal at k = 560
early  = t_win <= tSplit;
late   = t_win >  tSplit;

A_560_1 = max(Ex560_win(early));    % first pulse (inside slab)
A_560_2 = max(Ex560_win(late));     % second pulse (reflected inside slab)
A_620    = max(Ex620_win);          % transmitted pulse after slab

fprintf('\nBone marrow peaks from simulation (window 22.5–62.5 ns):\n');
fprintf('  First peak at k=560  ~ %.4f\n', A_560_1);
fprintf('  Second peak at k=560 ~ %.4f\n', A_560_2);
fprintf('  Peak at k=620        ~ %.4f\n', A_620);

% Theory for comparison (single-interface coefficients)
% for losless media: eta = sqrt(mu/eps)
eps_r = 7.2103;
eta0  = sqrt(mu0/eps0);               
eta2  = sqrt(mu0/(eps_r*eps0));

% at interface of media: gamma = (eta2-eta1)/(eta2+eta1)
% at boundary: tau = 2eta2(eta2+eta1)
Gamma12 = (eta2 - eta0)/(eta2 + eta0);      % free->bone
Tau12   = 2*eta2/(eta2 + eta0);
Gamma23 = (eta0 - eta2)/(eta0 + eta2);      % bone->free
Tau23   = 2*eta0/(eta0 + eta2);

fprintf('\nTheory (for comparison):\n');
fprintf('  tau12              = %.4f\n', Tau12);
fprintf('  |Gamma23*tau12|    = %.4f\n', abs(Gamma23*Tau12));
fprintf('  |Tau23*tau12|      = %.4f\n', abs(Tau23*Tau12));