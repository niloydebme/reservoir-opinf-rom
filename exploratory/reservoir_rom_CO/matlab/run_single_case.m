function run_single_case(k_mD, S_skin, schedule_id, out_root, mrst_root)
%RUN_SINGLE_CASE  Parameterised MRST simulation for reservoir ROM training.
%
%  USAGE
%    run_single_case(k_mD, S_skin, schedule_id)
%    run_single_case(k_mD, S_skin, schedule_id, out_root)
%    run_single_case(k_mD, S_skin, schedule_id, out_root, mrst_root)
%
%  VARIABLE INPUTS
%    k_mD        : permeability [mD]          e.g. 5, 10, 20, 50, 100, 200, 500
%    S_skin      : skin factor  [-]           e.g. 0, 2, 5, 8, 10, 12, 15
%    schedule_id : rate schedule ID  (1–7)    see utils/make_rate_schedule.m
%
%  FIXED PARAMETERS (identical across all training runs):
%    phi = 0.18        ct  = 1.2e-5 psi-1    mu  = 1.2 cp
%    Bo  = 1.15        h   = 150 ft           rw  = 0.35 ft
%    C   = 0.003 bbl/psi                      pi  = 4500 psia
%    Grid: 50×50 cells, dx = dy = 30 m  →  1500 m × 1500 m
%    350 log-spaced output steps, t = 0.001 to 5000 hr
%
%  OUTPUT folder:  <out_root>/k<kk>_s<ss>_sched<id>/
%    snapshots_psia.mat  —  X_psia [2500 × N_t]  (psia)
%    well_data.mat       —  bhp_psia, delta_p_psi, rate_STBd, time_hr
%    params.json         —  all scalars and arrays in field units

%% ── 0. Defaults ───────────────────────────────────────────────────────
if nargin < 4 || isempty(out_root)
    out_root = fullfile(fileparts(mfilename('fullpath')), ...
                        '..', 'data', '01_training_CO');
end
if nargin < 5 || isempty(mrst_root)
    mrst_root = 'E:\Datasets Module\MRST2025';
end

%% ── 0b. Initialise MRST (safe to call even if already loaded) ─────────
if isempty(which('mrstPath'))
    addpath(mrst_root);
    run(fullfile(mrst_root, 'startup.m'));
end
mrstModule add ad-core ad-blackoil ad-props

SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR, 'utils'));
fu = si_to_field();

%% ── 0c. Output directory ──────────────────────────────────────────────
%  Naming convention:  k<3-digit mD>_s<2-digit skin>_sched<id>
%  Examples:  k005_s00_sched1   k100_s08_sched4   k500_s15_sched7
case_name = sprintf('k%03d_s%02d_sched%d', round(k_mD), round(S_skin), schedule_id);
OUT       = fullfile(out_root, case_name);
if ~exist(OUT, 'dir'), mkdir(OUT); end

fprintf('\n');
fprintf('==========================================================\n');
fprintf(' run_single_case : %s\n', case_name);
fprintf('  k = %.1f mD  |  S = %.0f  |  Schedule ID = %d\n', ...
        k_mD, S_skin, schedule_id);
fprintf('==========================================================\n\n');

%% ── 1. Fixed reservoir parameters ────────────────────────────────────
phi       = 0.18;
ct_SI     = 1.740e-9;              % Pa-1   (= 1.2e-5 psi-1)
mu_SI     = 1.2e-3;               % Pa·s   (= 1.2 cp)
Bo        = 1.15;                  % RB/STB
rho_SI    = 780;                   % kg/m3
h_m       = 45.72;                 % m      (= 150 ft)
rw_m      = 0.1067;               % m      (= 0.35 ft)
pi_SI     = 310.3 * barsa;        % Pa     (= 4500 psia)

C_bbl_psi = 0.003;
C_SI      = C_bbl_psi * 0.158987 / 6894.76;   % m3/Pa  = 6.918e-8

%% ── 2. Variable parameters ────────────────────────────────────────────
k_SI       = k_mD   * milli * darcy;
S_skin_val = S_skin;

%% ── 3. Grid (identical for all cases) ────────────────────────────────
nx = 50;  ny = 50;
dx = 30;  dy = 30;   % metres per cell

G  = cartGrid([nx, ny, 1], [nx*dx, ny*dy, h_m]);
G  = computeGeometry(G);

re_m  = (nx * dx) / 2;
re_ft = re_m * fu.m_to_ft;

fprintf('Grid    : %d×%d = %d cells | re = %.0f ft\n', ...
        nx, ny, G.cells.num, re_ft);

%% ── 4. Rock, fluid, model ────────────────────────────────────────────
%  rock.perm controls both transmissibility AND well index (via addWell),
%  so a single makeRock call propagates k through the entire simulation.
rock  = makeRock(G, k_SI, phi);

fluid = initSimpleADIFluid( ...
    'phases', 'W',   ...
    'mu',     mu_SI, ...
    'rho',    rho_SI,...
    'c',      ct_SI, ...
    'pRef',   pi_SI);

model  = GenericBlackOilModel(G, rock, fluid, ...
          'water', true, 'oil', false, 'gas', false);
state0 = initResSol(G, pi_SI, 1.0);

%% ── 5. Rate schedule ──────────────────────────────────────────────────
sched_info = make_rate_schedule(schedule_id);   % uses defaults: 5000 hr, 350 steps
n_steps    = sched_info.n_steps;
dt_vec_s   = sched_info.dt_vec_s;
ctrl_vec   = sched_info.ctrl_vec;

fprintf('Schedule: [%d] %s\n', schedule_id, sched_info.label);
fprintf('  %d steps | t = %.4f to %.0f hr\n', ...
        n_steps, sched_info.t_start_hr, sched_info.t_end_hr);
if ~isempty(sched_info.t_change_hr)
    fprintf('  Rate changes at: ');
    fprintf('%.0f hr  ', sched_info.t_change_hr);
    fprintf('\n');
end
fprintf('  Rates: ');
fprintf('%d STB/d  ', sched_info.rates_STBd);
fprintf('\n');

%% ── 6. Build one MRST well per control period ────────────────────────
%  addWell computes well index automatically from rock.perm,
%  so changing k automatically gives the correct WI without any manual formula.
wc = (round(ny/2) - 1)*nx + round(nx/2);   % central cell (MATLAB 1-indexed)

W_controls = cell(sched_info.n_controls, 1);
for ic = 1 : sched_info.n_controls
    q_rate = sched_info.rates_STBd(ic);
    % Convert surface rate to reservoir m3/s (negative = production)
    q_SI_ctrl = -(q_rate * Bo * 0.158987) / 86400;
    W_controls{ic} = addWell([], G, rock, wc, ...
    'Type',   'rate',        ...
    'Val',    q_SI_ctrl,     ...
    'Radius', rw_m,          ...
    'Skin',   S_skin_val,    ...
    'Name',   'PROD1',       ...
    'Comp_i', 1);

    W_controls{ic}.Cs = 0 ;   % wellbore storage is neglected in this simulation
    fprintf('  Control %d: q = %d STB/day  (WI = %.4e m3/(Pa.s))\n', ...
            ic, q_rate, W_controls{ic}.WI);
end

%% ── 7. Assemble MRST schedule struct ─────────────────────────────────
%  schedule.step.val(j)     = timestep size for step j       [s]
%  schedule.step.control(j) = which W to use at step j       [1-indexed]
%  schedule.control(k).W    = well configuration for control k
schedule.step.val     = dt_vec_s;
schedule.step.control = ctrl_vec;
for ic = 1 : sched_info.n_controls
    schedule.control(ic).W = W_controls{ic};
end

%% ── 8. Analytical diagnostics (for verification only) ────────────────
kh_mDft = k_mD * (h_m * fu.m_to_ft);
mu_cp   = mu_SI  * fu.Pas_to_cp;
ct_psi1 = ct_SI  * fu.Pam1_to_psi1;
pi_psia = pi_SI  * fu.Pa_to_psia;

q_ref     = sched_info.rates_STBd(1);   % use first period for diagnostics
t_wbs_hr  = 170000 * C_bbl_psi * exp(2*S_skin_val) / kh_mDft;
t_pss_hr  = (phi * mu_cp * ct_psi1 * re_ft^2) / (0.000264 * k_mD) * 0.1;
dp_prime  = 70.6  * q_ref * mu_cp * Bo / kh_mDft;
m_slope   = 162.6 * q_ref * mu_cp * Bo / kh_mDft;

fprintf('\n--- Analytical flow regime diagnostics ---\n');
fprintf('  kh = %.0f mD-ft\n', kh_mDft);
fprintf('  WBS ends  ~ %.4f hr\n',  t_wbs_hr);
fprintf('  MTR ends  ~ %.1f hr\n',  t_pss_hr);
fprintf('  dp prime  = %.4f psi  (Bourdet flat level)\n', dp_prime);
fprintf('  MTR slope = %.4f psi/cycle\n', m_slope);
if t_wbs_hr > sched_info.t_end_hr
    fprintf('  WARNING: WBS end time (%.1f hr) exceeds simulation end.\n', t_wbs_hr);
    fprintf('           High skin dominates — MTR may not be visible.\n');
end
fprintf('------------------------------------------\n\n');

%% ── 9. Run simulation ────────────────────────────────────────────────
fprintf('Running MRST simulation...\n');
t_wall = tic;
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
t_elapsed = toc(t_wall);
n_t = numel(states);
fprintf('  Completed %d / %d steps in %.1f s\n', n_t, n_steps, t_elapsed);

if n_t < n_steps
    warning('run_single_case: only %d of %d steps completed for %s', ...
            n_t, n_steps, case_name);
end

%% ── 10. Extract well solution ────────────────────────────────────────
time_hr   = cumsum(dt_vec_s(1:n_t)) / 3600;
bhp_psia  = zeros(1, n_t);
rate_STBd = zeros(1, n_t);

for k_idx = 1 : n_t
    ws = states{k_idx}.wellSol(1);
    bhp_psia(k_idx)  =  ws.bhp     * fu.Pa_to_psia;
    rate_STBd(k_idx) =  abs(ws.qWs) * fu.m3s_to_STBd;
end
delta_p_psi = pi_psia - bhp_psia;

% Sanity checks
if any(bhp_psia > pi_psia + 1)
    warning('run_single_case: %d steps have BHP > pi. Check inputs.', ...
            sum(bhp_psia > pi_psia));
end
if any(bhp_psia < 0)
    warning('run_single_case: %d steps have BHP < 0 psia. Simulation may be unphysical.', ...
            sum(bhp_psia < 0));
end

fprintf('\nWell solution summary:\n');
fprintf('  BHP range  : %.1f – %.1f psia\n', min(bhp_psia),   max(bhp_psia));
fprintf('  Delta-p    : %.4f – %.2f psi\n',  min(delta_p_psi),max(delta_p_psi));
fprintf('  Rate range : %.1f – %.1f STB/day\n', min(rate_STBd), max(rate_STBd));

% Regime step counts (using first-period diagnostics)
n_wbs = sum(time_hr < t_wbs_hr);
n_mtr = sum(time_hr >= t_wbs_hr & time_hr < t_pss_hr);
n_bdf = sum(time_hr >= t_pss_hr);
fprintf('  Regime steps:  WBS=%d  MTR=%d  BDF=%d\n', n_wbs, n_mtr, n_bdf);

%% ── 11. Export ───────────────────────────────────────────────────────
fprintf('\nExporting data...\n');

% (A) Snapshot matrix  [N_cells × N_t]  in psia
X_psia = zeros(G.cells.num, n_t);
for k_idx = 1 : n_t
    X_psia(:, k_idx) = states{k_idx}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUT, 'snapshots_psia.mat'), 'X_psia', '-v7.3');
fprintf('  snapshots_psia.mat  [%d × %d]\n', G.cells.num, n_t);

% (B) Well data
save(fullfile(OUT, 'well_data.mat'), ...
     'bhp_psia', 'delta_p_psi', 'rate_STBd', 'time_hr', '-v7.3');
fprintf('  well_data.mat       [%d timesteps]\n', n_t);

% (C) params.json  —  everything needed by Python notebooks
p = struct();
p.case_name          = case_name;
p.k_mD               = k_mD;
p.S_skin             = S_skin_val;
p.schedule_id        = schedule_id;
p.schedule_label     = sched_info.label;
p.rates_STBd         = sched_info.rates_STBd(:)';   % row for clean JSON
p.t_change_hr        = sched_info.t_change_hr(:)';
p.phi                = phi;
p.ct_psi1            = ct_psi1;
p.mu_cp              = mu_cp;
p.Bo                 = Bo;
p.h_ft               = h_m  * fu.m_to_ft;
p.rw_ft              = rw_m * fu.m_to_ft;
p.C_bbl_psi          = 0;                     % changed to zero
p.p_init_psia        = pi_psia;
p.kh_mDft            = kh_mDft;
p.re_ft              = re_ft;
p.nx                 = nx;
p.ny                 = ny;
p.dx_ft              = dx * fu.m_to_ft;
p.n_timesteps        = n_t;
p.t_start_hr         = sched_info.t_start_hr;
p.t_end_hr           = sched_info.t_end_hr;
p.t_wbs_end_hr       = t_wbs_hr;
p.t_pss_start_hr     = t_pss_hr;
p.dp_prime_psi       = dp_prime;
p.m_slope_psi_cycle  = m_slope;
p.well_cell_matlab   = wc;
p.sim_wall_time_s    = t_elapsed;

write_json(fullfile(OUT, 'params.json'), p);
fprintf('  params.json\n');

fprintf('\n==========================================================\n');
fprintf(' DONE: %s  (%.1f s wall time)\n', case_name, t_elapsed);
fprintf('==========================================================\n\n');

end