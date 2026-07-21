%==========================================================================
%  case4_variable_rate.m  —  MRST 2025b
%
%  Multi-rate variable production, homogeneous reservoir.
%  Same reservoir, fluid and well as Case 1 (homogeneous k=20 mD) for
%  direct comparison.  The rate schedule introduces three rate changes
%  and a final shut-in (pressure buildup), enabling Rate-Transient
%  Analysis (RTA).
%
%  WHY variable rate (departure from Cases 1–3):
%    Cases 1–3 used constant rate to isolate individual flow regimes.
%    Case 4 tests the ROM under realistic, non-constant excitation.
%
%  Rate schedule (4 periods, 350 total log-spaced steps):
%    Period 1: q =  800 STB/day,  0.001–1000 hr  (WBS + RF + BDF onset)
%    Period 2: q =  400 STB/day,  1000–2500 hr   (rate drop, BDF)
%    Period 3: q = 1200 STB/day,  2500–4000 hr   (rate surge, BDF)
%    Period 4: q =    0 STB/day,  4000–5000 hr   (shut-in buildup)
%
%  Reservoir and fluid parameters identical to Case 1:
%    k=20 mD  phi=0.18  ct=1.2e-5 psi-1  mu=1.2 cp  Bo=1.15
%    h=150 ft  rw=0.35 ft  S=2  C=0.003 bbl/psi
%    Grid: 50x50, dx=dy=30 m  ->  1500x1500 m
%    pi=4500 psia
%
%  Regime boundaries (from Case 1 analysis, k=20 mD):
%    t_WBS ~  9.28 hr   t_PSS ~ 297 hr
%
%  Output: data/01_synthetic/case4_variable_rate/
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
fprintf(' Case 4 — Multi-Rate Variable Production\n');
fprintf(' Homogeneous k=20 mD | 4 periods | RTA\n');
fprintf('==========================================================\n\n');

SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR, 'utils'));
fu = si_to_field();

OUT = fullfile(SCRIPT_DIR,'..','data','01_synthetic','case4_variable_rate');
if ~exist(OUT,'dir'), mkdir(OUT); end
fprintf('Output  : %s\n\n', OUT);

%% ── 1. Parameters (identical to Case 1) ─────────────────────────────────
phi       = 0.18;
ct_SI     = 1.740e-9;
mu_SI     = 1.2e-3;
Bo        = 1.15;
rho_SI    = 780;
pi_SI     = 310.3 * barsa;
h_m       = 45.72;
rw_m      = 0.1067;
S_skin    = 2.0;
C_bbl_psi = 0.003;
C_SI      = C_bbl_psi * 0.158987 / 6894.76;

%  Rate values for each period [m3/s, negative = production]
q1_STBd =  800;   q1_SI = -(q1_STBd * Bo * 0.158987) / 86400;
q2_STBd =  400;   q2_SI = -(q2_STBd * Bo * 0.158987) / 86400;
q3_STBd = 1200;   q3_SI = -(q3_STBd * Bo * 0.158987) / 86400;
q4_SI   = -1e-10;   % or -1e-12;     % shut-in

%% ── 2. Grid : 50x50, dx=dy=30 m (identical to Cases 1 & 2) ─────────────
nx = 50;  ny = 50;
dx = 30;  dy = 30;
G  = cartGrid([nx, ny, 1], [nx*dx, ny*dy, h_m]);
G  = computeGeometry(G);
re_m  = (nx*dx) / 2;
re_ft = re_m * fu.m_to_ft;

%% ── 3. Homogeneous rock (identical to Case 1) ────────────────────────────
k_mD      = 20.0;
k_SI      = k_mD * milli * darcy;
rock      = makeRock(G, k_SI, phi);

fprintf('Rock    : k=%.0f mD  phi=%.2f (homogeneous, identical to Case 1)\n', ...
        k_mD, phi);

%% ── 4. Fluid, model, initial state ──────────────────────────────────────
fluid  = initSimpleADIFluid('phases','W','mu',mu_SI,'rho',rho_SI, ...
                            'c',ct_SI,'pRef',pi_SI);
model  = GenericBlackOilModel(G, rock, fluid, ...
         'water',true,'oil',false,'gas',false);
state0 = initResSol(G, pi_SI, 1.0);
fprintf('Model   : %s\n', class(model));

%% ── 5. Well (identical to Case 1) ───────────────────────────────────────
wc = (round(ny/2)-1)*nx + round(nx/2);

%  Four well controls — one per period
W1 = addWell([],G,rock,wc,'Type','rate','Val',q1_SI,'Radius',rw_m, ...
             'Skin',S_skin,'Cs',C_SI,'Name','PROD','Comp_i',1);
W2 = addWell([],G,rock,wc,'Type','rate','Val',q2_SI,'Radius',rw_m, ...
             'Skin',S_skin,'Cs',C_SI,'Name','PROD','Comp_i',1);
W3 = addWell([],G,rock,wc,'Type','rate','Val',q3_SI,'Radius',rw_m, ...
             'Skin',S_skin,'Cs',C_SI,'Name','PROD','Comp_i',1);
W4 = addWell([],G,rock,wc,'Type','rate','Val',q4_SI,'Radius',rw_m, ...
             'Skin',S_skin,'Cs',C_SI,'Name','PROD','Comp_i',1);
% Period 4 uses BHP control at pi to simulate shut-in (q->0)

fprintf('Well    : cell %d | S=%.0f | C=%.3f bbl/psi\n', wc, S_skin, C_bbl_psi);

%% ── 6. Variable-rate schedule ────────────────────────────────────────────
%
%  Each period has its own log-spaced timestep vector.
%  Period end-time of previous period becomes start of next period.
%
%  Period 1: 0.001  – 1000 hr  (150 steps)  q = 800 STB/day
%  Period 2: 1000   – 2500 hr  ( 75 steps)  q = 400 STB/day
%  Period 3: 2500   – 4000 hr  ( 75 steps)  q = 1200 STB/day
%  Period 4: 4000   – 5000 hr  ( 50 steps)  q = 0 (shut-in)
%
periods_hr = [0.001,  1000, 150;    % [t_start, t_end, n_steps]
              1000,   2500,  75;
              2500,   4000,  75;
              4000,   5000,  50];

n_periods   = size(periods_hr, 1);
ctrl_idx    = zeros(1, 0);   % control index for each timestep
dt_all      = zeros(1, 0);
elapsed_all = zeros(1, 0);

q_vals = [q1_STBd, q2_STBd, q3_STBd, 0];

for pp = 1:n_periods
    t_s   = periods_hr(pp, 1);
    t_e   = periods_hr(pp, 2);
    n_st  = periods_hr(pp, 3);

    t_nodes = logspace(log10(t_s*3600), log10(t_e*3600), n_st+1)';
    dt_p    = diff(t_nodes);
    dt_all  = [dt_all; dt_p];
    ctrl_idx = [ctrl_idx, pp*ones(1, n_st)];
    fprintf('Period %d : %6g – %5g hr  (%d steps)  q=%4d STB/day\n', ...
            pp, t_s, t_e, n_st, ...
            q_vals(pp));
end

%  Build MRST schedule manually
W_ctrl = {W1, W2, W3, W4};
step_ctrl = num2cell(ctrl_idx');   % each step points to a control
schedule.control = struct('W', W_ctrl);
schedule.step.val     = dt_all;
schedule.step.control = ctrl_idx';

fprintf('\nTotal   : %d timesteps  |  t_end = %.0f hr\n', ...
        numel(dt_all), sum(dt_all)/3600);

%% ── 7. Analytical regime references (k=20 mD, Case 1 basis) ─────────────
mu_cp    = mu_SI  * fu.Pas_to_cp;
ct_psi1  = ct_SI  * fu.Pam1_to_psi1;
pi_psia  = pi_SI  * fu.Pa_to_psia;
h_ft     = h_m    * fu.m_to_ft;
kh_mDft  = k_mD   * h_ft;

t_wbs_hr = 170000 * C_bbl_psi * mu_cp * exp(0.14*S_skin) / kh_mDft;
t_pss_hr = (phi * mu_cp * ct_psi1 * re_ft^2) / (0.000264*k_mD) * 0.1;
dp_flat  = 70.6  * q1_STBd * mu_cp * Bo / kh_mDft;
m_slope  = 162.6 * q1_STBd * mu_cp * Bo / kh_mDft;

fprintf('\nRegime boundaries (k=%.0f mD, Case 1 basis):\n', k_mD);
fprintf('  t_WBS = %.2f hr   t_PSS = %.1f hr\n', t_wbs_hr, t_pss_hr);
fprintf('  Bourdet flat (q1=%d): dp'' = %.4f psi\n', q1_STBd, dp_flat);
fprintf('  MTR slope    (q1=%d): m   = %.4f psi/cycle\n', q1_STBd, m_slope);

% Rate-normalised flat levels for each period
dp_flat_q2 = 70.6 * q2_STBd * mu_cp * Bo / kh_mDft;
dp_flat_q3 = 70.6 * q3_STBd * mu_cp * Bo / kh_mDft;
fprintf('  Bourdet flat (q2=%d): %.4f psi\n', q2_STBd, dp_flat_q2);
fprintf('  Bourdet flat (q3=%d): %.4f psi\n', q3_STBd, dp_flat_q3);

%% ── 8. Simulate ──────────────────────────────────────────────────────────
fprintf('\nRunning simulation...\n');
t_wall = tic;
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
n_t = numel(states);
fprintf('  Completed %d / %d timesteps in %.1f s\n', ...
        n_t, numel(dt_all), toc(t_wall));

%% ── 9. Extract well data ─────────────────────────────────────────────────
time_hr   = cumsum(dt_all(1:n_t)) * fu.s_to_hr;
bhp_psia  = zeros(1, n_t);
rate_STBd = zeros(1, n_t);
for k_idx = 1:n_t
    ws = states{k_idx}.wellSol(1);
    bhp_psia(k_idx)  = ws.bhp     * fu.Pa_to_psia;
    rate_STBd(k_idx) = abs(ws.qWs) * fu.m3s_to_STBd;
end

delta_p_psi   = pi_psia - bhp_psia;
cum_prod_STB  = cumtrapz(time_hr, rate_STBd);    % Np [STB]

% Superposition time (multi-rate): t_sup for rate-normalised PTA
% Using Agarwal equivalent time for last period before buildup
rate_schedule = [q1_STBd, q2_STBd, q3_STBd, 0];
period_start_hr = [periods_hr(:,1)]';

fprintf('\nBHP     : %.1f – %.1f psia\n', min(bhp_psia), max(bhp_psia));
fprintf('Cum prod: %.1f STB\n', cum_prod_STB(end));
fprintf('Rate STBd: %.1f – %.1f\n', min(rate_STBd), max(rate_STBd));

%% ── 10. Export ───────────────────────────────────────────────────────────
fprintf('\nExporting...\n');

X_psia = zeros(G.cells.num, n_t);
for k_idx = 1:n_t
    X_psia(:, k_idx) = states{k_idx}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUT,'snapshots_psia.mat'),'X_psia','-v7.3');
fprintf('  snapshots_psia.mat  [%d x %d]\n', G.cells.num, n_t);

% Well data — includes rate array (key difference from Cases 1–3)
save(fullfile(OUT,'well_data.mat'), ...
     'bhp_psia','delta_p_psi','rate_STBd','time_hr','cum_prod_STB', ...
     'rate_schedule','period_start_hr','-v7.3');
fprintf('  well_data.mat       [%d timesteps, rate history included]\n', n_t);

perm_mD = rock.perm * fu.m2_to_mD;
poro    = rock.poro;
save(fullfile(OUT,'rock_field.mat'),'perm_mD','poro','-v7.3');

centroids_ft = G.cells.centroids * fu.m_to_ft;
save(fullfile(OUT,'grid_field.mat'),'centroids_ft','-v7.3');

% Peaceman correction (identical to Case 1)
dx_ft = dx * fu.m_to_ft;
r_eq_ft_p = 0.2 * dx_ft;

p.description         = 'Case 4: multi-rate variable production, k=20 mD homogeneous';
p.p_init_psia         = pi_psia;
p.k_mD                = k_mD;
p.k_well_cell_mD      = k_mD;       % homogeneous: same as k
p.phi                 = phi;
p.ct_psi1             = ct_psi1;
p.mu_cp               = mu_cp;
p.Bo                  = Bo;
p.h_ft                = h_ft;
p.rw_ft               = rw_m * fu.m_to_ft;
p.S_skin              = S_skin;
p.C_bbl_psi           = C_bbl_psi;
p.q1_STBd             = q1_STBd;    % reference rate (Period 1)
p.q2_STBd             = q2_STBd;
p.q3_STBd             = q3_STBd;
p.kh_mDft             = kh_mDft;
p.re_ft               = re_ft;
p.nx                  = nx;
p.ny                  = ny;
p.dx_ft               = dx_ft;
p.n_timesteps         = n_t;
p.t_wbs_end_hr        = t_wbs_hr;
p.t_pss_start_hr      = t_pss_hr;
p.dp_prime_psi        = dp_flat;    % Bourdet flat for q1
p.dp_prime_q2_psi     = dp_flat_q2;
p.dp_prime_q3_psi     = dp_flat_q3;
p.m_slope_psi_cycle   = m_slope;
p.t_period1_end_hr    = periods_hr(1,2);
p.t_period2_end_hr    = periods_hr(2,2);
p.t_period3_end_hr    = periods_hr(3,2);
p.t_period4_end_hr    = periods_hr(4,2);
p.t_shutin_hr         = periods_hr(4,1);
p.well_cell_matlab    = wc;

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
fprintf(' Case 4 complete\n');
fprintf(' k=%.0f mD  phi=%.2f  S=%.0f  C=%.3f bbl/psi\n', k_mD,phi,S_skin,C_bbl_psi);
fprintf(' Period 1: q=%d STB/d  0.001–1000 hr\n', q1_STBd);
fprintf(' Period 2: q=%d STB/d  1000–2500 hr\n', q2_STBd);
fprintf(' Period 3: q=%d STB/d  2500–4000 hr\n', q3_STBd);
fprintf(' Period 4: shut-in      4000–5000 hr\n');
fprintf(' Output: %s\n', OUT);
fprintf('==========================================================\n');
