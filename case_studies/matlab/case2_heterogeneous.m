%==========================================================================
%  case2_heterogeneous.m  —  MRST 2025b
%
%  Highly heterogeneous radial reservoir, single-phase, constant-rate
%  drawdown. Identical fluid, well and grid parameters to Case 1.
%  The permeability field is strongly heterogeneous: log-normal with
%  sigma_lnk = 2.0  (V_DP = 0.865 — strongly heterogeneous).
%
%  Permeability field:
%    Distribution : ln k ~ N(mu_lnk, sigma_lnk^2)
%    k_geo        = 20 mD  (same as Case 1 homogeneous value)
%    sigma_lnk    = 2.0
%    V_DP         = 1 - exp(-2.0) = 0.865
%    k range      ~ 0.4 – 1100 mD  (contrast ~ 3000)
%
%  All other parameters identical to Case 1:
%    phi = 0.18     ct = 1.2e-5 psi-1   mu = 1.2 cp    Bo = 1.15
%    h   = 150 ft   rw = 0.35 ft        S  = 2          C  = 0.003 bbl/psi
%    q   = 800 STB/day   pi = 4500 psia
%    Grid: 50x50, dx=dy=30 m  ->  1500 m x 1500 m
%    Schedule: 350 log-spaced steps, t = 0.001 to 5000 hr
%
%  Flow regime boundaries (based on k_geo):
%    WBS    :  t = 0.001 – ~
%    MTR    :  t = ~  – ~297 hr   (effective radial, 1.5 decades)
%    BDF    :  t = ~297  – 5000 hr
%
%  Output: data/01_synthetic/case2_heterogeneous/
%    snapshots_psia.mat, well_data.mat, rock_field.mat,
%    grid_field.mat, params.json
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
fprintf(' Case 2 — Highly Heterogeneous Reservoir\n');
fprintf(' sigma_lnk=2.0  V_DP=0.865  k_geo=20 mD\n');
fprintf('==========================================================\n\n');

SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR, 'utils'));
fu = si_to_field();

OUT = fullfile(SCRIPT_DIR, '..', 'data', '01_synthetic', 'case2_heterogeneous');
if ~exist(OUT, 'dir'), mkdir(OUT); end
fprintf('Output  : %s\n\n', OUT);

%% ── 1. Fluid and rock parameters (identical to Case 1) ──────────────────
phi       = 0.18;
ct_SI     = 1.740e-9;              % Pa-1  (= 1.2e-5 psi-1)
mu_SI     = 1.2e-3;                % 1.2 cp
Bo        = 1.15;
rho_SI    = 780;                   % kg/m3
pi_SI     = 310.3 * barsa;         % ~4500 psia
h_m       = 45.72;                 % m  (= 150 ft)
rw_m      = 0.1067;                % m  (= 0.35 ft)
S_skin    = 2.0;
C_bbl_psi = 0.003;
C_SI      = C_bbl_psi * 0.158987 / 6894.76;   % 6.918e-8 m3/Pa
q_STBd    = 800;
q_SI      = -(q_STBd * Bo * 0.158987) / 86400;

%% ── 2. Grid : 50x50, dx=dy=30 m (identical to Case 1) ───────────────────
nx = 50;  ny = 50;
dx = 30;  dy = 30;
G  = cartGrid([nx, ny, 1], [nx*dx, ny*dy, h_m]);
G  = computeGeometry(G);
re_m  = (nx*dx)/2;
re_ft = re_m * fu.m_to_ft;

%% ── 3. Highly heterogeneous permeability field ───────────────────────────
%
%  Log-normal distribution:
%    k_geo = 20 mD  (geometric mean = Case 1 value, for fair comparison)
%    sigma_lnk = 2.0  =>  V_DP = 0.865  (strongly heterogeneous)
%    k range ~ 0.4 to 1100 mD  (contrast ~ 3000)
%
%  Seed = 42 for reproducibility.
%
k_geo_mD  = 20.0;
k_geo_SI  = k_geo_mD * milli * darcy;
sigma_lnk = 2.0;
mu_lnk    = log(k_geo_SI);

rng(42);
ln_k_field = mu_lnk + sigma_lnk * randn(G.cells.num, 1);
k_field_SI = exp(ln_k_field);

%  Clamp permeability to physical bounds [1 mD, 2000 mD]
k_min_SI   = 1    * milli * darcy;
k_max_SI   = 2000 * milli * darcy;
k_field_SI = max(k_min_SI, min(k_max_SI, k_field_SI));

%  Statistics
k_field_mD  = k_field_SI * fu.m2_to_mD;
k_geo_actual = exp(mean(log(k_field_mD)));
sigma_actual = std(log(k_field_mD));
V_DP         = 1 - exp(-sigma_actual);
k_ratio      = max(k_field_mD) / min(k_field_mD);

fprintf('Permeability field (highly heterogeneous):\n');
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

%% ── 5. Well (identical to Case 1: rate, skin, WBS) ──────────────────────
wc = (round(ny/2) - 1)*nx + round(nx/2);

W = addWell([], G, rock, wc,      ...
    'Type',   'rate',              ...
    'Val',    q_SI,                ...
    'Radius', rw_m,                ...
    'Skin',   S_skin,              ...
    'Cs',     C_SI,                ...
    'Name',   'PROD1',             ...
    'Comp_i', 1);

fprintf('Well    : cell %d | q = %d STB/day | S = %.0f | C = %.3f bbl/psi\n', ...
        wc, q_STBd, S_skin, C_bbl_psi);

%% ── 6. Schedule : 350 log-spaced elapsed times, 0.001 to 5000 hr ─────────
t_start_hr = 0.001;
t_end_hr   = 5000;
n_steps    = 350;
t_nodes_s  = logspace(log10(t_start_hr*3600), ...
                       log10(t_end_hr  *3600), n_steps)';
dt_vec_s   = diff([0; t_nodes_s]);
schedule   = simpleSchedule(dt_vec_s, 'W', W);

fprintf('Schedule: %d steps | t = %.4f to %.0f hr\n', ...
        n_steps, t_start_hr, t_end_hr);

%% ── Analytical regime boundaries (based on k_geo) ────────────────────────
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

%% ── 7. Simulate ──────────────────────────────────────────────────────────
fprintf('\nRunning simulation...\n');
t_wall = tic;
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
n_t = numel(states);
fprintf('  Completed %d / %d timesteps in %.1f s\n', ...
        n_t, n_steps, toc(t_wall));

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
fprintf('\nExporting...\n');

X_psia = zeros(G.cells.num, n_t);
for k_idx = 1:n_t
    X_psia(:,k_idx) = states{k_idx}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUT,'snapshots_psia.mat'),'X_psia','-v7.3');
fprintf('  snapshots_psia.mat  [%d x %d] psia\n', G.cells.num, n_t);

save(fullfile(OUT,'well_data.mat'), ...
     'bhp_psia','delta_p_psi','rate_STBd','time_hr','-v7.3');
fprintf('  well_data.mat       [%d timesteps]\n', n_t);

perm_mD = rock.perm * fu.m2_to_mD;
poro    = rock.poro;
save(fullfile(OUT,'rock_field.mat'),'perm_mD','poro','-v7.3');

centroids_ft = G.cells.centroids * fu.m_to_ft;
save(fullfile(OUT,'grid_field.mat'),'centroids_ft','-v7.3');

%  params.json — same keys as Case 1 + heterogeneity additions
p.description        = 'Case 2: highly heterogeneous, sigma_lnk=2.0, V_DP=0.865';
p.p_init_psia        = pi_psia;
p.k_mD               = k_geo_mD;        % geometric mean (Case 1 equivalent)
p.k_geo_mD           = k_geo_actual;
p.k_min_mD           = min(k_field_mD);
p.k_max_mD           = max(k_field_mD);
p.sigma_lnk          = sigma_actual;
p.k_well_cell_mD = perm_mD(wc);   % actual permeability at well cell [mD]
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
p.nx                 = nx;
p.ny                 = ny;
p.dx_ft              = dx * fu.m_to_ft;
p.n_timesteps        = n_t;
p.t_start_hr         = t_nodes_s(1)   / 3600;
p.t_end_hr           = t_nodes_s(end) / 3600;
p.t_wbs_end_hr       = t_wbs_hr;
p.t_pss_start_hr     = t_pss_hr;
p.dp_prime_psi       = dp_prime;
p.m_slope_psi_cycle  = m_slope;
p.well_cell_matlab   = wc;

fid = fopen(fullfile(OUT,'params.json'),'w');
fn  = fieldnames(p);
fprintf(fid,'{\n');
for i = 1:numel(fn)
    v = p.(fn{i});
    if ischar(v), fprintf(fid,'  "%s": "%s"',fn{i},v);
    else,         fprintf(fid,'  "%s": %.10g',fn{i},v); end
    if i < numel(fn), fprintf(fid,','); end
    fprintf(fid,'\n');
end
fprintf(fid,'}\n');
fclose(fid);
fprintf('  params.json\n');

fprintf('\n==========================================================\n');
fprintf(' Case 2 complete\n');
fprintf(' k_geo=%.1f mD  sigma_lnk=%.2f  V_DP=%.3f\n', ...
        k_geo_mD, sigma_actual, V_DP);
fprintf(' WBS  : 0.001 – %.2f hr  (%d steps)\n', t_wbs_hr, n_wbs);
fprintf(' MTR  : %.2f – %.1f hr  (%d steps)\n',  t_wbs_hr, t_pss_hr, n_mtr);
fprintf(' BDF  : %.1f – %.0f hr  (%d steps)\n',  t_pss_hr, t_end_hr, n_bdf);
fprintf(' Output: %s\n', OUT);
fprintf('==========================================================\n');