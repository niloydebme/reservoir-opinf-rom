%==========================================================================
%  case1_homogeneous.m  —  MRST 2025b
%
%  Homogeneous radial reservoir, single-phase, constant-rate drawdown.
%  Three flow regimes are present in the dataset:
%
%    1. Wellbore Storage (WBS)   : t = 0.001 – ~9.3 hr
%       Wellbore unloading dominates. BHP barely changes; Bourdet
%       derivative rises on a unit slope on the log-log plot.
%
%    2. Middle Time Region (MTR) : t = ~9.3 – ~1189 hr  (2.1 decades)
%       Infinite-acting radial flow. Semi-log plot shows a straight line
%       with slope m = 59.84 psi/cycle. Bourdet derivative is flat at
%       dp' = 25.98 psi. This is the primary PTA diagnostic window.
%
%    3. Boundary-Dominated Flow  : t = ~1189 – 5000 hr
%       Closed boundary reached. BHP falls steeply below the MTR line.
%       Bourdet derivative rises above the flat level (slope -> 1).
%
%  Reservoir (all quantities exported in oilfield field units):
%    k   = 20 mD      phi = 0.18     ct = 1.2e-5 psi-1   mu = 1.2 cp
%    h   = 150 ft     rw  = 0.35 ft  Bo = 1.15 RB/STB
%    q   = 800 STB/d  S   = 2        C  = 0.003 bbl/psi
%    pi  = 4500 psia
%    Grid: 50x50, dx=dy=30 m  ->  1500 m x 1500 m
%    Schedule: 350 log-spaced steps, t = 0.001 to 5000 hr
%
%  Output files (DATA directory):
%    snapshots_psia.mat   — X_psia  [N_cells x N_t]  psia
%    well_data.mat        — bhp_psia, delta_p_psi, rate_STBd, time_hr
%    rock_field.mat       — perm_mD, poro
%    grid_field.mat       — centroids_ft
%    params.json          — all scalars in field units
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
fprintf(' Case 1 — Homogeneous Radial Reservoir\n');
fprintf(' WBS  |  Radial Flow  |  Boundary-Dominated Flow\n');
fprintf('==========================================================\n\n');

SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR, 'utils'));
fu = si_to_field();

OUT = fullfile(SCRIPT_DIR, '..', 'data', '01_synthetic', 'case1_homogeneous');
if ~exist(OUT, 'dir'), mkdir(OUT); end
fprintf('Output  : %s\n\n', OUT);

%% ── 1. Reservoir parameters (all defined in SI, converted at export) ────

%  Permeability and rock
k_SI    = 20    * milli * darcy;    %  20 mD
phi     = 0.18;                     %  18 % porosity
ct_SI   = 1.740e-9;                 %  Pa-1  (= 1.2e-5 psi-1)

%  Fluid
mu_SI   = 1.2e-3;                   %  1.2 cp  (moderate viscosity oil)
Bo      = 1.15;                     %  1.15 RB/STB
rho_SI  = 780;                      %  kg/m3

%  Geometry
h_m     = 45.72;                    %  m  (= 150 ft)
rw_m    = 0.1067;                   %  m  (= 0.35 ft)
pi_SI   = 310.3 * barsa;            %  Pa  (= 4500 psia)

%  Well control
q_STBd  = 800;                      %  STB/day surface rate
S_skin  = 2.0;                      %  skin factor

%  Wellbore storage  C = 0.003 bbl/psi  ->  SI  [m3/Pa]
C_bbl_psi = 0.003;
C_SI      = C_bbl_psi * 0.158987 / 6894.76;   %  = 6.918e-8  m3/Pa

fprintf('Reservoir parameters:\n');
fprintf('  k   = %.0f mD      phi  = %.2f       ct   = %.2e psi-1\n', ...
        k_SI*fu.m2_to_mD, phi, ct_SI*fu.Pam1_to_psi1);
fprintf('  mu  = %.1f cp      Bo   = %.2f       h    = %.0f ft\n', ...
        mu_SI*fu.Pas_to_cp, Bo, h_m*fu.m_to_ft);
fprintf('  rw  = %.3f ft     S    = %.0f         C    = %.3f bbl/psi\n', ...
        rw_m*fu.m_to_ft, S_skin, C_bbl_psi);
fprintf('  q   = %d STB/day  pi   = %.0f psia\n', q_STBd, pi_SI*fu.Pa_to_psia);

%% ── 2. Grid : 50x50 Cartesian, dx=dy=30 m  ->  1500 m x 1500 m ─────────
nx = 50;  ny = 50;
dx = 30;  dy = 30;   % metres per cell
G  = cartGrid([nx, ny, 1], [nx*dx, ny*dy, h_m]);
G  = computeGeometry(G);

re_m  = (nx * dx) / 2;              % effective drainage radius
re_ft = re_m * fu.m_to_ft;

fprintf('\nGrid    : %dx%d = %d cells | %.0fx%.0f m | re = %.0f ft\n', ...
        nx, ny, G.cells.num, nx*dx, ny*dy, re_ft);

%% ── 3. Rock and fluid ───────────────────────────────────────────────────
rock  = makeRock(G, k_SI, phi);

fluid = initSimpleADIFluid( ...
    'phases', 'W',     ...
    'mu',     mu_SI,   ...
    'rho',    rho_SI,  ...
    'c',      ct_SI,   ...
    'pRef',   pi_SI);

%% ── 4. Model and initial conditions ─────────────────────────────────────
model  = GenericBlackOilModel(G, rock, fluid, ...
         'water', true, 'oil', false, 'gas', false);
state0 = initResSol(G, pi_SI, 1.0);
fprintf('Model   : %s\n', class(model));

%% ── 5. Well : centre cell, constant rate, skin + wellbore storage ────────
%
%  Surface rate in SI  (negative sign = production)
%    q [m3/s] = q [STB/day] * Bo [RB/STB] * 0.158987 [m3/RB] / 86400 [s/day]
q_SI = -(q_STBd * Bo * 0.158987) / 86400;

%  Central grid cell (MATLAB 1-indexed)
wc = (round(ny/2) - 1)*nx + round(nx/2);

W = addWell([], G, rock, wc,       ...
    'Type',   'rate',               ...
    'Val',    q_SI,                 ...
    'Radius', rw_m,                 ...
    'Skin',   S_skin,               ...  skin factor S = 2
    'Cs',     C_SI,                 ...  wellbore storage [m3/Pa]
    'Name',   'PROD1',              ...
    'Comp_i', 1);

fprintf('\nWell    : cell %d | q = %d STB/day | S = %.0f | C = %.3f bbl/psi\n', ...
        wc, q_STBd, S_skin, C_bbl_psi);

%% ── 6. Schedule : 350 log-spaced elapsed times, 0.001 to 5000 hr ────────
%
%  KEY DESIGN DECISION:
%    t_nodes  = the ELAPSED simulation time at each output step [seconds]
%    dt_vec   = diff([0; t_nodes]) = the actual timestep SIZE for each step
%
%  This is NOT the same as logspace(dt_min, dt_max, n).
%  Using diff(elapsed) guarantees:
%    (a) the first timestep covers exactly 0.001 hr from t=0
%    (b) the final cumulative time equals t_end_hr exactly
%    (c) timestep sizes grow smoothly (no sudden jumps)
%
%  The cumsum(logspace(...)) approach that appears in some MRST tutorials
%  produces a different end-time from what was intended and should be avoided.

t_start_hr  = 0.001;   % hours  — deep inside WBS regime
t_end_hr    = 5000;    % hours  — well into BDF
n_steps     = 350;

t_nodes_s   = logspace(log10(t_start_hr*3600), ...
                        log10(t_end_hr  *3600), n_steps)';   % [s]
dt_vec_s    = diff([0; t_nodes_s]);                          % [s]

schedule    = simpleSchedule(dt_vec_s, 'W', W);

fprintf('Schedule: %d steps | t = %.4f to %.0f hr\n', ...
        n_steps, t_start_hr, t_end_hr);
fprintf('         dt_min = %.5f hr | dt_max = %.1f hr\n', ...
        dt_vec_s(1)/3600, dt_vec_s(end)/3600);

%% ── Analytical regime boundaries (pre-simulation verification) ───────────
kh_mDft     = (k_SI * fu.m2_to_mD) * (h_m * fu.m_to_ft);
mu_cp       = mu_SI  * fu.Pas_to_cp;
ct_psi1     = ct_SI  * fu.Pam1_to_psi1;
pi_psia     = pi_SI  * fu.Pa_to_psia;

t_wbs_hr    = 170000 * C_bbl_psi * mu_cp * exp(0.14*S_skin) / kh_mDft;
t_pss_hr    = (phi * mu_cp * ct_psi1 * re_ft^2) / ...
              (0.000264 * (k_SI * fu.m2_to_mD)) * 0.1;
dp_prime    = 70.6  * q_STBd * mu_cp * Bo / kh_mDft;
m_slope     = 162.6 * q_STBd * mu_cp * Bo / kh_mDft;
mtr_decades = log10(t_pss_hr / t_wbs_hr);

fprintf('\n--- Analytical flow regime boundaries ---\n');
fprintf('  kh              = %.0f mD-ft\n', kh_mDft);
fprintf('  re              = %.0f ft\n',    re_ft);
fprintf('\n  WBS ends at     t ~ %.2f hr\n',  t_wbs_hr);
fprintf('  MTR window      t = %.2f – %.1f hr  (%.2f decades)\n', ...
        t_wbs_hr, t_pss_hr, mtr_decades);
fprintf('  PSS starts at   t ~ %.1f hr  (%.1f days)\n', ...
        t_pss_hr, t_pss_hr/24);
fprintf('\n  Bourdet flat    dp'' = %.4f psi\n', dp_prime);
fprintf('  MTR slope       m   = %.4f psi/cycle\n', m_slope);
fprintf('  kh from slope   kh_est = 162.6*q*mu*Bo/m = %.0f mD-ft\n', kh_mDft);
fprintf('-----------------------------------------\n\n');

%% ── 7. Run simulation ────────────────────────────────────────────────────
fprintf('Running simulation...\n');
t_wall = tic;
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
t_elapsed = toc(t_wall);
n_t = numel(states);
fprintf('  Completed %d / %d timesteps in %.1f s\n', n_t, n_steps, t_elapsed);

%% ── 8. Extract well solution ─────────────────────────────────────────────
time_hr   = cumsum(dt_vec_s(1:n_t)) * fu.s_to_hr;
bhp_psia  = zeros(1, n_t);
rate_STBd = zeros(1, n_t);

for k_idx = 1:n_t
    ws = states{k_idx}.wellSol(1);
    bhp_psia(k_idx)  =  ws.bhp     * fu.Pa_to_psia;
    rate_STBd(k_idx) =  abs(ws.qWs) * fu.m3s_to_STBd;
end

delta_p_psi = pi_psia - bhp_psia;

fprintf('\nWell solution summary:\n');
fprintf('  BHP       : %.1f – %.1f psia\n',  min(bhp_psia),   max(bhp_psia));
fprintf('  Delta-p   : %.4f – %.2f psi\n',   min(delta_p_psi),max(delta_p_psi));
fprintf('  Rate (avg): %.1f STB/day\n',       mean(rate_STBd));

% Sanity check: count steps in each regime
n_wbs = sum(time_hr < t_wbs_hr);
n_mtr = sum(time_hr >= t_wbs_hr & time_hr < t_pss_hr);
n_bdf = sum(time_hr >= t_pss_hr);
fprintf('\n  Timesteps per regime:\n');
fprintf('  WBS (%d)   : t < %.2f hr\n',    n_wbs, t_wbs_hr);
fprintf('  MTR (%d)  : t = %.2f – %.1f hr\n', n_mtr, t_wbs_hr, t_pss_hr);
fprintf('  BDF (%d)   : t > %.1f hr\n',   n_bdf, t_pss_hr);

%% ── 9. Export ────────────────────────────────────────────────────────────
fprintf('\nExporting data...\n');

% (A) Snapshot matrix  [N_cells x N_t]  in psia
X_psia = zeros(G.cells.num, n_t);
for k_idx = 1:n_t
    X_psia(:, k_idx) = states{k_idx}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUT, 'snapshots_psia.mat'), 'X_psia', '-v7.3');
fprintf('  snapshots_psia.mat  [%d x %d] psia\n', G.cells.num, n_t);

% (B) Well data
save(fullfile(OUT, 'well_data.mat'), ...
     'bhp_psia', 'delta_p_psi', 'rate_STBd', 'time_hr', '-v7.3');
fprintf('  well_data.mat       [%d timesteps]\n', n_t);

% (C) Rock
perm_mD = rock.perm * fu.m2_to_mD;
poro    = rock.poro;
save(fullfile(OUT, 'rock_field.mat'), 'perm_mD', 'poro', '-v7.3');

% (D) Grid centroids in ft
centroids_ft = G.cells.centroids * fu.m_to_ft;
save(fullfile(OUT, 'grid_field.mat'), 'centroids_ft', '-v7.3');

% (E) params.json  — every scalar in oilfield field units
p.description        = 'Case 1: homogeneous radial, WBS+MTR+BDF, field units';
p.p_init_psia        = pi_psia;
p.k_mD               = k_SI    * fu.m2_to_mD;
p.phi                = phi;
p.ct_psi1            = ct_psi1;
p.mu_cp              = mu_cp;
p.Bo                 = Bo;
p.h_ft               = h_m     * fu.m_to_ft;
p.rw_ft              = rw_m    * fu.m_to_ft;
p.S_skin             = S_skin;
p.C_bbl_psi          = C_bbl_psi;
p.q_STBd             = q_STBd;
p.kh_mDft            = kh_mDft;
p.re_ft              = re_ft;
p.nx                 = nx;
p.ny                 = ny;
p.dx_ft              = dx       * fu.m_to_ft;
p.n_timesteps        = n_t;
p.t_start_hr         = t_nodes_s(1)   / 3600;
p.t_end_hr           = t_nodes_s(end) / 3600;
p.t_wbs_end_hr       = t_wbs_hr;
p.t_pss_start_hr     = t_pss_hr;
p.dp_prime_psi       = dp_prime;
p.m_slope_psi_cycle  = m_slope;
p.well_cell_matlab   = wc;

fid = fopen(fullfile(OUT, 'params.json'), 'w');
fn  = fieldnames(p);
fprintf(fid, '{\n');
for i = 1:numel(fn)
    v = p.(fn{i});
    if ischar(v)
        fprintf(fid, '  "%s": "%s"', fn{i}, v);
    else
        fprintf(fid, '  "%s": %.10g', fn{i}, v);
    end
    if i < numel(fn), fprintf(fid, ','); end
    fprintf(fid, '\n');
end
fprintf(fid, '}\n');
fclose(fid);
fprintf('  params.json         [all scalars in field units]\n');

%% ── 10. Summary ──────────────────────────────────────────────────────────
fprintf('\n==========================================================\n');
fprintf(' Simulation complete\n');
fprintf('==========================================================\n');
fprintf(' WBS    : t = 0.001  – %.2f hr  (%d steps)\n', t_wbs_hr, n_wbs);
fprintf(' MTR    : t = %.2f  – %.1f hr  (%d steps)\n', t_wbs_hr, t_pss_hr, n_mtr);
fprintf(' BDF    : t = %.1f  – %.0f hr  (%d steps)\n', t_pss_hr, t_end_hr, n_bdf);
fprintf(' kh     : %.0f mD-ft  |  dp'' = %.4f psi  |  m = %.4f psi/cycle\n', ...
        kh_mDft, dp_prime, m_slope);
fprintf(' BHP at t_end : %.1f psia  (pi = %.0f, Dp = %.1f psi)\n', ...
        min(bhp_psia), pi_psia, max(delta_p_psi));
fprintf(' Output : %s\n', OUT);
fprintf('==========================================================\n');