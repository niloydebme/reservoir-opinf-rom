%==========================================================================
%  case2_finemesh.m  —  MRST 2025b
%
%  Highly heterogeneous reservoir — FINE MESH / HIGH RESOLUTION version
%  of Case 2 (case2_heterogeneous.m), created for computational complexity
%  benchmarking in Part 2 of the MSc thesis.
%
%  PURPOSE
%  -------
%  Demonstrate the computational cost of high-resolution simulation versus
%  the coarse reference (Case 2), and show that POD achieves similar
%  reconstruction accuracy on BOTH grids, validating the method's
%  mesh-independence and its value as a compression / ROM tool.
%
%  RESOLUTION COMPARISON
%  ----------------------
%  |  Parameter        |  Coarse (Case 2)     |  Fine (this script)  |
%  |-------------------|----------------------|----------------------|
%  |  Grid             |  50 x 50 cells       |  200 x 200 cells     |
%  |  Cell size dx=dy  |  30 m                |  7.5 m               |
%  |  N_cells          |  2,500               |  40,000  (x16)       |
%  |  Timesteps        |  350 (log-spaced)    |  1,500 (log-spaced)  |
%  |  Snapshot matrix  |  2,500 x 350         |  40,000 x 1,500      |
%  |  Domain size      |  1500 m x 1500 m     |  1500 m x 1500 m     |
%  |  Time range       |  0.001 – 5000 hr     |  0.001 – 5000 hr     |
%
%  PERMEABILITY FIELD
%  ------------------
%  The coarse 50x50 log-normal field (rng=42, sigma_lnk=2.0, k_geo=20 mD)
%  is generated first, then bilinearly interpolated onto the 200x200 fine
%  grid. This preserves the spatial heterogeneity structure so that both
%  grids sample the same underlying random field.
%
%  All fluid/well/rock parameters are IDENTICAL to Case 2:
%    phi = 0.18     ct = 1.2e-5 psi-1   mu = 1.2 cp    Bo = 1.15
%    h   = 150 ft   rw = 0.35 ft        S  = 2          C  = 0.003 bbl/psi
%    q   = 800 STB/day   pi = 4500 psia
%    sigma_lnk = 2.0  (V_DP = 0.865)  k_geo = 20 mD
%
%  Output: data/01_synthetic/case2_finemesh/
%    snapshots_psia.mat, well_data.mat, rock_field.mat,
%    grid_field.mat, params.json, timing.json
%==========================================================================

%% ── 0. Initialise MRST ──────────────────────────────────────────────────
mrst_root = 'E:\Datasets Module\MRST2025';
if isempty(which('mrstPath'))
    addpath(mrst_root);
    run(fullfile(mrst_root, 'startup.m'));
end
mrstModule add ad-core ad-blackoil ad-props
clear; clc;

fprintf('\n');
fprintf('==========================================================\n');
fprintf(' Case 2 FINE MESH — Highly Heterogeneous Reservoir\n');
fprintf(' Grid: 200x200 | 1500 timesteps | For Complexity Benchmark\n');
fprintf(' sigma_lnk=2.0  V_DP=0.865  k_geo=20 mD\n');
fprintf('==========================================================\n\n');

SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR, 'utils'));
fu = si_to_field();

OUT = fullfile(SCRIPT_DIR, '..', 'data', '01_synthetic', 'case2_finemesh');
if ~exist(OUT, 'dir'), mkdir(OUT); end
fprintf('Output  : %s\n\n', OUT);

% Start wall-clock timer for total simulation cost reporting
t_total_wall = tic;

%% ── 1. Fluid and rock parameters (identical to Case 2 coarse) ───────────
phi       = 0.18;
ct_SI     = 1.740e-9;              % Pa-1  (= 1.2e-5 psi-1)
mu_SI     = 1.2e-3;                % Pa.s  (1.2 cp)
Bo        = 1.15;
rho_SI    = 780;                   % kg/m3
pi_SI     = 310.3 * barsa;         % ~4500 psia
h_m       = 45.72;                 % m  (= 150 ft)
rw_m      = 0.1067;                % m  (= 0.35 ft)
S_skin    = 2.0;
C_bbl_psi = 0.003;
C_SI      = C_bbl_psi * 0.158987 / 6894.76;  % m3/Pa
q_STBd    = 800;
q_SI      = -(q_STBd * Bo * 0.158987) / 86400;

%% ── 2. FINE grid : 200x200, dx=dy=7.5 m (same 1500 m x 1500 m domain) ──
%
%  Coarse: 50 cells * 30 m/cell = 1500 m
%  Fine  : 200 cells *  7.5 m/cell = 1500 m  [unchanged domain]
%  Refinement factor: 4x in each direction → 16x total cells
%
nx_fine = 200;  ny_fine = 200;
dx_fine = 7.5;  dy_fine = 7.5;    % m (= 1500 m / 200)

G = cartGrid([nx_fine, ny_fine, 1], [nx_fine*dx_fine, ny_fine*dy_fine, h_m]);
G = computeGeometry(G);

re_m  = (nx_fine * dx_fine) / 2;  % drainage radius (unchanged: 750 m)
re_ft = re_m * fu.m_to_ft;

fprintf('Fine grid       : %dx%d = %d cells\n', nx_fine, ny_fine, G.cells.num);
fprintf('Cell size       : %.2f m x %.2f m\n',   dx_fine, dy_fine);
fprintf('Domain          : %.0f m x %.0f m (unchanged)\n', ...
        nx_fine*dx_fine, ny_fine*dy_fine);

%% ── 3. Permeability field: interpolate coarse 50x50 → fine 200x200 ───────
%
%  Strategy: generate the identical coarse 50x50 log-normal field (same
%  rng seed = 42, same sigma_lnk = 2.0, same k_geo = 20 mD) and then
%  bilinearly interpolate onto the fine Cartesian grid. This ensures both
%  grids sample the SAME underlying heterogeneous structure, making the
%  computational comparison physically meaningful.
%
nx_c = 50;  ny_c = 50;   % coarse dimensions (Case 2)
dx_c = 30;  dy_c = 30;   % coarse cell size [m]

k_geo_mD  = 20.0;
k_geo_SI  = k_geo_mD * milli * darcy;
sigma_lnk = 2.0;
mu_lnk    = log(k_geo_SI);

% ── Generate SAME coarse log-normal field (rng=42, identical to Case 2) ──
rng(42);
ln_k_coarse = mu_lnk + sigma_lnk * randn(nx_c * ny_c, 1);

% Reshape to 2-D matrix (MATLAB column-major: rows=y-cells, cols=x-cells)
ln_k_coarse_2D = reshape(ln_k_coarse, nx_c, ny_c);   % [nx_c x ny_c]

% Cell-centre coordinates on the coarse grid [m]
xc_c = (0.5:nx_c) * dx_c;   % 15, 45, 75 … 1485  (50 values)
yc_c = (0.5:ny_c) * dy_c;   % 15, 45, 75 … 1485  (50 values)

% Cell-centre coordinates on the fine grid [m]
xc_f = (0.5:nx_fine) * dx_fine;   % 3.75, 11.25 … 1496.25  (200 values)
yc_f = (0.5:ny_fine) * dy_fine;

% Bilinear interpolation: coarse field → fine grid
% interp2 operates on [Y,X] convention when input is a grid
[Xc_grid, Yc_grid] = meshgrid(xc_c, yc_c);   % [50 x 50] each
[Xf_grid, Yf_grid] = meshgrid(xc_f, yc_f);   % [200 x 200] each

% ln_k_coarse_2D in MATLAB layout: (i = x-index, j = y-index)
% interp2 expects (X, Y, V, Xq, Yq) where V rows correspond to Y
ln_k_fine_2D = interp2(xc_c, yc_c, ln_k_coarse_2D', Xf_grid, Yf_grid, 'linear');
% Extrapolate boundary with nearest-neighbour to avoid NaN at edges
ln_k_fine_2D = fillmissing(ln_k_fine_2D, 'nearest');

% Reshape back to column vector in MATLAB cell ordering (x varies fastest)
ln_k_fine = reshape(ln_k_fine_2D', nx_fine * ny_fine, 1);
k_field_SI = exp(ln_k_fine);

% ── Clamp to physical bounds (same as Case 2) ────────────────────────────
k_min_SI   = 1    * milli * darcy;
k_max_SI   = 2000 * milli * darcy;
k_field_SI = max(k_min_SI, min(k_max_SI, k_field_SI));

% ── Statistics ───────────────────────────────────────────────────────────
k_field_mD   = k_field_SI * fu.m2_to_mD;
k_geo_actual = exp(mean(log(k_field_mD)));
sigma_actual = std(log(k_field_mD));
V_DP         = 1 - exp(-sigma_actual);
k_ratio      = max(k_field_mD) / min(k_field_mD);

fprintf('\nPermeability field (interpolated from coarse rng=42):\n');
fprintf('  sigma_lnk (target)  = %.2f\n',    sigma_lnk);
fprintf('  sigma_lnk (actual)  = %.4f\n',    sigma_actual);
fprintf('  k_geo     (target)  = %.1f mD\n', k_geo_mD);
fprintf('  k_geo     (actual)  = %.2f mD\n', k_geo_actual);
fprintf('  V_DP                = %.3f\n',     V_DP);
fprintf('  k range             = %.2f – %.1f mD\n', ...
        min(k_field_mD), max(k_field_mD));
fprintf('  Contrast ratio      = %.0f\n',     k_ratio);

%% ── 4. Rock, fluid, model ────────────────────────────────────────────────
rock   = makeRock(G, k_field_SI, phi);
fluid  = initSimpleADIFluid('phases','W','mu',mu_SI,'rho',rho_SI, ...
                            'c',ct_SI,'pRef',pi_SI);
model  = GenericBlackOilModel(G, rock, fluid, ...
         'water',true,'oil',false,'gas',false);
state0 = initResSol(G, pi_SI, 1.0);
fprintf('\nModel   : %s\n', class(model));

%% ── 5. Well (same location, rate, skin, WBS as Case 2) ──────────────────
%  Well placed at domain centre; map coarse cell index (1338, 0-based)
%  to fine grid. Coarse centre = (25,25) → fine centre = (100,100).
wc = (round(ny_fine/2) - 1) * nx_fine + round(nx_fine/2);

W = addWell([], G, rock, wc,       ...
    'Type',   'rate',               ...
    'Val',    q_SI,                 ...
    'Radius', rw_m,                 ...
    'Skin',   S_skin,               ...
    'Cs',     C_SI,                 ...
    'Name',   'PROD1',              ...
    'Comp_i', 1);

fprintf('Well    : cell %d | q = %d STB/day | S = %.0f | C = %.3f bbl/psi\n', ...
        wc, q_STBd, S_skin, C_bbl_psi);

%% ── 6. Schedule : 1500 log-spaced steps, t = 0.001 to 5000 hr ───────────
%  Same time range as Case 2; 4.3x more timesteps for higher temporal
%  resolution and increased computational cost.
%
t_start_hr = 0.001;
t_end_hr   = 5000;
n_steps    = 1500;
t_nodes_s  = logspace(log10(t_start_hr * 3600), ...
                       log10(t_end_hr   * 3600), n_steps)';
dt_vec_s   = diff([0; t_nodes_s]);
schedule   = simpleSchedule(dt_vec_s, 'W', W);

fprintf('Schedule: %d steps | t = %.4f to %.0f hr\n', ...
        n_steps, t_start_hr, t_end_hr);

%% ── Analytical regime boundaries (based on k_geo, identical to Case 2) ──
kh_mDft   = k_geo_mD * (h_m * fu.m_to_ft);
mu_cp     = mu_SI  * fu.Pas_to_cp;
ct_psi1   = ct_SI  * fu.Pam1_to_psi1;
pi_psia   = pi_SI  * fu.Pa_to_psia;

t_wbs_hr  = 170000 * C_bbl_psi * mu_cp * exp(0.14*S_skin) / kh_mDft;
t_pss_hr  = (phi * mu_cp * ct_psi1 * re_ft^2) / ...
            (0.000264 * k_geo_mD) * 0.1;
dp_prime  = 70.6  * q_STBd * mu_cp * Bo / kh_mDft;
m_slope   = 162.6 * q_STBd * mu_cp * Bo / kh_mDft;

fprintf('\nRegime boundaries (k_geo = %.1f mD):\n', k_geo_mD);
fprintf('  WBS ends    : t ~ %.2f hr\n',  t_wbs_hr);
fprintf('  PSS starts  : t ~ %.1f hr\n',   t_pss_hr);
fprintf('  Bourdet flat: dp'' = %.4f psi\n', dp_prime);
fprintf('  MTR slope   : m   = %.4f psi/cycle\n', m_slope);

%% ── Memory estimate ──────────────────────────────────────────────────────
mem_GB = (G.cells.num * n_steps * 8) / 1e9;
fprintf('\nSnapshot matrix  : %d x %d  (%.2f GB at float64)\n', ...
        G.cells.num, n_steps, mem_GB);

%% ── 7. Simulate ──────────────────────────────────────────────────────────
fprintf('\nRunning FINE MESH simulation (this may take several minutes)...\n');
t_sim_wall = tic;
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
sim_wall_s     = toc(t_sim_wall);
n_t            = numel(states);
fprintf('  Completed %d / %d timesteps in %.1f s  (%.1f min)\n', ...
        n_t, n_steps, sim_wall_s, sim_wall_s/60);

%% ── 8. Extract well solution ─────────────────────────────────────────────
time_hr   = cumsum(dt_vec_s(1:n_t)) * fu.s_to_hr;
bhp_psia  = zeros(1, n_t);
rate_STBd = zeros(1, n_t);
for k_idx = 1:n_t
    ws = states{k_idx}.wellSol(1);
    bhp_psia(k_idx)  = ws.bhp     * fu.Pa_to_psia;
    rate_STBd(k_idx) = abs(ws.qWs) * fu.m3s_to_STBd;
end
delta_p_psi = pi_psia - bhp_psia;

fprintf('\nBHP     : %.1f – %.1f psia\n',  min(bhp_psia),    max(bhp_psia));
fprintf('Delta-p : %.4f – %.2f psi\n',   min(delta_p_psi), max(delta_p_psi));

n_wbs = sum(time_hr < t_wbs_hr);
n_mtr = sum(time_hr >= t_wbs_hr & time_hr < t_pss_hr);
n_bdf = sum(time_hr >= t_pss_hr);
fprintf('WBS (%d)  MTR (%d)  BDF (%d)\n', n_wbs, n_mtr, n_bdf);

%% ── 9. Export ────────────────────────────────────────────────────────────
fprintf('\nExporting outputs...\n');
t_export_wall = tic;

% Snapshot matrix (large – requires HDF5/v7.3)
X_psia = zeros(G.cells.num, n_t);
for k_idx = 1:n_t
    X_psia(:,k_idx) = states{k_idx}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUT,'snapshots_psia.mat'),'X_psia','-v7.3');
fprintf('  snapshots_psia.mat  [%d x %d]  %.2f GB\n', ...
        G.cells.num, n_t, (G.cells.num*n_t*8)/1e9);

save(fullfile(OUT,'well_data.mat'), ...
     'bhp_psia','delta_p_psi','rate_STBd','time_hr','-v7.3');
fprintf('  well_data.mat       [%d timesteps]\n', n_t);

perm_mD = rock.perm * fu.m2_to_mD;
poro    = rock.poro;
save(fullfile(OUT,'rock_field.mat'),'perm_mD','poro','-v7.3');

centroids_ft = G.cells.centroids * fu.m_to_ft;
save(fullfile(OUT,'grid_field.mat'),'centroids_ft','-v7.3');

export_wall_s = toc(t_export_wall);
total_wall_s  = toc(t_total_wall);

%% ── 10. params.json ──────────────────────────────────────────────────────
p.description        = 'Case 2 FINE MESH: 200x200, 1500 steps, same domain/physics as coarse';
p.p_init_psia        = pi_psia;
p.k_mD               = k_geo_mD;
p.k_geo_mD           = k_geo_actual;
p.k_min_mD           = min(k_field_mD);
p.k_max_mD           = max(k_field_mD);
p.sigma_lnk          = sigma_actual;
p.k_well_cell_mD     = perm_mD(wc);
p.V_DP               = V_DP;
p.phi                = phi;
p.ct_psi1            = ct_psi1;
p.mu_cp              = mu_cp;
p.Bo                 = Bo;
p.h_ft               = h_m  * fu.m_to_ft;
p.rw_ft              = rw_m * fu.m_to_ft;
p.S_skin             = S_skin;
p.C_bbl_psi          = C_bbl_psi;
p.q_STBd             = q_STBd;
p.kh_mDft            = kh_mDft;
p.re_ft              = re_ft;
p.nx                 = nx_fine;
p.ny                 = ny_fine;
p.dx_ft              = dx_fine * fu.m_to_ft;
p.n_timesteps        = n_t;
p.t_start_hr         = t_nodes_s(1)   / 3600;
p.t_end_hr           = t_nodes_s(end) / 3600;
p.t_wbs_end_hr       = t_wbs_hr;
p.t_pss_start_hr     = t_pss_hr;
p.dp_prime_psi       = dp_prime;
p.m_slope_psi_cycle  = m_slope;
p.well_cell_matlab   = wc;
% Fine-mesh specific entries for comparison script
p.mesh_type          = 'fine';
p.coarse_nx          = nx_c;
p.coarse_ny          = ny_c;
p.refinement_factor  = nx_fine / nx_c;   % = 4
p.n_cells            = G.cells.num;
p.n_timesteps_target = n_steps;
p.perm_interp_method = 'bilinear_from_coarse_rng42';

fid = fopen(fullfile(OUT,'params.json'),'w');
fn  = fieldnames(p);
fprintf(fid,'{\n');
for i = 1:numel(fn)
    v = p.(fn{i});
    if ischar(v), fprintf(fid,'  "%s": "%s"', fn{i}, v);
    else,         fprintf(fid,'  "%s": %.10g', fn{i}, v);
    end
    if i < numel(fn), fprintf(fid,','); end
    fprintf(fid,'\n');
end
fprintf(fid,'}\n');
fclose(fid);
fprintf('  params.json\n');

%% ── 11. timing.json ──────────────────────────────────────────────────────
%  Separate timing record for the complexity comparison script.
timing.mesh_type           = 'fine';
timing.nx                  = nx_fine;
timing.ny                  = ny_fine;
timing.n_cells             = G.cells.num;
timing.n_timesteps         = n_t;
timing.sim_wall_s          = sim_wall_s;
timing.export_wall_s       = export_wall_s;
timing.total_wall_s        = total_wall_s;
timing.sim_per_step_s      = sim_wall_s / n_t;
timing.snapshot_matrix_GB  = (G.cells.num * n_t * 8) / 1e9;

fid2 = fopen(fullfile(OUT,'timing.json'),'w');
fn2  = fieldnames(timing);
fprintf(fid2,'{\n');
for i = 1:numel(fn2)
    v = timing.(fn2{i});
    if ischar(v), fprintf(fid2,'  "%s": "%s"', fn2{i}, v);
    else,         fprintf(fid2,'  "%s": %.6g',  fn2{i}, v);
    end
    if i < numel(fn2), fprintf(fid2,','); end
    fprintf(fid2,'\n');
end
fprintf(fid2,'}\n');
fclose(fid2);
fprintf('  timing.json\n');

%% ── Final summary ────────────────────────────────────────────────────────
fprintf('\n==========================================================\n');
fprintf(' Case 2 FINE MESH — Complete\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Grid         : %dx%d = %d cells (coarse: 50x50 = 2500)\n', ...
        nx_fine, ny_fine, G.cells.num);
fprintf(' Refinement   : 4x spatial  |  %.1fx temporal\n', n_steps/350);
fprintf(' Timesteps    : %d (coarse: 350)\n', n_t);
fprintf(' Sim. time    : %.1f s  (%.1f min)\n', sim_wall_s, sim_wall_s/60);
fprintf(' Total time   : %.1f s  (%.1f min)\n', total_wall_s, total_wall_s/60);
fprintf(' Output size  : %.2f GB\n', (G.cells.num*n_t*8)/1e9);
fprintf('----------------------------------------------------------\n');
fprintf(' k_geo=%.1f mD  sigma_lnk=%.2f  V_DP=%.3f\n', ...
        k_geo_mD, sigma_actual, V_DP);
fprintf(' WBS  : 0.001 – %.2f hr  (%d steps)\n', t_wbs_hr, n_wbs);
fprintf(' MTR  : %.2f – %.1f hr  (%d steps)\n',  t_wbs_hr, t_pss_hr, n_mtr);
fprintf(' BDF  : %.1f – %.0f hr  (%d steps)\n',  t_pss_hr, t_end_hr, n_bdf);
fprintf(' Output: %s\n', OUT);
fprintf('==========================================================\n');