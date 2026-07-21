function run_single_case(k_mD, S_skin, schedule_id, out_root, mrst_root)
%RUN_SINGLE_CASE  Parameterised MRST simulation for reservoir ROM training.
%  (Instrumented version — records wall time, memory, and step timing.)
%
%  USAGE
%    run_single_case(k_mD, S_skin, schedule_id)
%    run_single_case(k_mD, S_skin, schedule_id, out_root)
%    run_single_case(k_mD, S_skin, schedule_id, out_root, mrst_root)
%
%  VARIABLE INPUTS
%    k_mD        : permeability [mD]          e.g. 5, 10, 20, 50, 100, 200, 500
%    S_skin      : skin factor  [-]           e.g. 0, 2, 5, 8, 10, 12, 15
%    schedule_id : rate schedule ID  (1-7)    see utils/make_rate_schedule.m
%
%  FIXED PARAMETERS (identical across all training runs):
%    phi = 0.18        ct  = 1.2e-5 psi-1    mu  = 1.2 cp
%    Bo  = 1.15        h   = 150 ft           rw  = 0.35 ft
%    C   = 0.003 bbl/psi                      pi  = 4500 psia
%    Grid: 50x50 cells, dx = dy = 30 m  ->  1500 m x 1500 m
%    350 log-spaced output steps, t = 0.001 to 5000 hr
%
%  PROFILING OUTPUT (written to params.json and printed to console):
%    sim_wall_time_s   : total simulateScheduleAD wall time [s]
%    setup_wall_time_s : grid/rock/fluid/model setup time [s]
%    export_wall_time_s: data export time [s]
%    mem_before_MB     : memory before simulation [MB]
%    mem_after_MB      : memory after simulation [MB]
%    mem_peak_delta_MB : estimated peak delta (after - before) [MB]
%    steps_completed   : actual timesteps completed
%    steps_requested   : requested timesteps
%    wall_per_step_ms  : mean wall time per completed step [ms]

%% ── 0. Defaults ───────────────────────────────────────────────────────
if nargin < 4 || isempty(out_root)
    out_root = fullfile(fileparts(mfilename('fullpath')), ...
                        '..', 'data', '01_training');
end
if nargin < 5 || isempty(mrst_root)
    mrst_root = 'E:\Datasets Module\MRST2025';
end

%% ── 0b. Initialise MRST ───────────────────────────────────────────────
if isempty(which('mrstPath'))
    addpath(mrst_root);
    run(fullfile(mrst_root, 'startup.m'));
end
mrstModule add ad-core ad-blackoil ad-props

SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR, 'utils'));
fu = si_to_field();

%% ── 0c. Output directory ──────────────────────────────────────────────
case_name = sprintf('k%03d_s%02d_sched%d', round(k_mD), round(S_skin), schedule_id);
OUT       = fullfile(out_root, case_name);
if ~exist(OUT, 'dir'), mkdir(OUT); end

fprintf('\n');
fprintf('==========================================================\n');
fprintf(' run_single_case : %s\n', case_name);
fprintf('  k = %.1f mD  |  S = %.0f  |  Schedule ID = %d\n', ...
        k_mD, S_skin, schedule_id);
fprintf('==========================================================\n\n');

%% ── *** PROFILING: overall timer starts here *** ─────────────────────
t_total_start = tic;

%% ── 1. Fixed reservoir parameters ────────────────────────────────────
phi       = 0.18;
ct_SI     = 1.740e-9;
mu_SI     = 1.2e-3;
Bo        = 1.15;
rho_SI    = 780;
h_m       = 45.72;
rw_m      = 0.1067;
pi_SI     = 310.3 * barsa;
C_bbl_psi = 0.003;
C_SI      = C_bbl_psi * 0.158987 / 6894.76;

%% ── 2. Variable parameters ────────────────────────────────────────────
k_SI       = k_mD   * milli * darcy;
S_skin_val = S_skin;

%% ── 3. Grid ───────────────────────────────────────────────────────────
t_setup_start = tic;

nx = 50;  ny = 50;
dx = 30;  dy = 30;
G  = cartGrid([nx, ny, 1], [nx*dx, ny*dy, h_m]);
G  = computeGeometry(G);

re_m  = (nx * dx) / 2;
re_ft = re_m * fu.m_to_ft;

fprintf('Grid    : %d x %d = %d cells | re = %.0f ft\n', ...
        nx, ny, G.cells.num, re_ft);

%% ── 4. Rock, fluid, model ────────────────────────────────────────────
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

setup_wall_time_s = toc(t_setup_start);
fprintf('Setup wall time : %.2f s\n', setup_wall_time_s);

%% ── 5. Rate schedule ──────────────────────────────────────────────────
sched_info = make_rate_schedule(schedule_id);
n_steps    = sched_info.n_steps;
dt_vec_s   = sched_info.dt_vec_s;
ctrl_vec   = sched_info.ctrl_vec;

fprintf('Schedule: [%d] %s\n', schedule_id, sched_info.label);
fprintf('  %d steps | t = %.4f to %.0f hr\n', ...
        n_steps, sched_info.t_start_hr, sched_info.t_end_hr);

%% ── 6. Build wells ────────────────────────────────────────────────────
wc = (round(ny/2) - 1)*nx + round(nx/2);
W_controls = cell(sched_info.n_controls, 1);
for ic = 1 : sched_info.n_controls
    q_rate    = sched_info.rates_STBd(ic);
    q_SI_ctrl = -(q_rate * Bo * 0.158987) / 86400;
    W_controls{ic} = addWell([], G, rock, wc, ...
        'Type',   'rate',        ...
        'Val',    q_SI_ctrl,     ...
        'Radius', rw_m,          ...
        'Skin',   S_skin_val,    ...
        'Name',   'PROD1',       ...
        'Comp_i', 1);
    W_controls{ic}.Cs = C_SI;
    fprintf('  Control %d: q = %d STB/day  (WI = %.4e m3/(Pa.s))\n', ...
            ic, q_rate, W_controls{ic}.WI);
end

%% ── 7. Assemble schedule ──────────────────────────────────────────────
schedule.step.val     = dt_vec_s;
schedule.step.control = ctrl_vec;
for ic = 1 : sched_info.n_controls
    schedule.control(ic).W = W_controls{ic};
end

%% ── 8. Analytical diagnostics ────────────────────────────────────────
kh_mDft   = k_mD * (h_m * fu.m_to_ft);
mu_cp     = mu_SI  * fu.Pas_to_cp;
ct_psi1   = ct_SI  * fu.Pam1_to_psi1;
pi_psia   = pi_SI  * fu.Pa_to_psia;
q_ref     = sched_info.rates_STBd(1);
t_wbs_hr  = 170000 * C_bbl_psi * mu_cp * exp(0.14*S_skin_val) / kh_mDft;
t_pss_hr  = (phi * mu_cp * ct_psi1 * re_ft^2) / (0.000264 * k_mD) * 0.1;
dp_prime  = 70.6  * q_ref * mu_cp * Bo / kh_mDft;
m_slope   = 162.6 * q_ref * mu_cp * Bo / kh_mDft;

fprintf('\n--- Analytical diagnostics ---\n');
fprintf('  WBS ends  ~ %.4f hr\n',  t_wbs_hr);
fprintf('  MTR ends  ~ %.1f hr\n',  t_pss_hr);
fprintf('  dp prime  = %.4f psi\n', dp_prime);
fprintf('  MTR slope = %.4f psi/cycle\n', m_slope);

%% ── 9. Memory snapshot BEFORE simulation ─────────────────────────────
% memstats() returns memory info on all platforms via whos/feature
try
    [~, mem_sys] = memory();
    mem_info_before = mem_sys.PhysicalMemory.Available / 1e6;
catch
    mem_info_before = NaN;  % memory() not supported on this platform
end
fprintf('\nMemory before simulation : %.1f MB (MATLAB workspace)\n', ...
        mem_info_before);

%% ── 10. Run simulation ───────────────────────────────────────────────
fprintf('\nRunning MRST simulation (%d steps)...\n', n_steps);
t_sim_start = tic;
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
sim_wall_time_s = toc(t_sim_start);

n_t = numel(states);
wall_per_step_ms = sim_wall_time_s / max(n_t, 1) * 1000;

fprintf('  Completed %d / %d steps\n', n_t, n_steps);
fprintf('  Simulation wall time : %.2f s  (%.0f ms/step avg)\n', ...
        sim_wall_time_s, wall_per_step_ms);

%% ── 11. Memory snapshot AFTER simulation ─────────────────────────────
try
    [~, mem_sys] = memory();
    mem_info_after = mem_sys.PhysicalMemory.Available / 1e6;
catch
    mem_info_after = NaN;
end
mem_peak_delta = mem_info_after - mem_info_before;
fprintf('Memory after  simulation : %.1f MB  (delta: %+.1f MB)\n', ...
        mem_info_after, mem_peak_delta);

if n_t < n_steps
    warning('run_single_case: only %d of %d steps completed for %s', ...
            n_t, n_steps, case_name);
end

%% ── 12. Extract well solution ────────────────────────────────────────
time_hr   = cumsum(dt_vec_s(1:n_t)) / 3600;
bhp_psia  = zeros(1, n_t);
rate_STBd = zeros(1, n_t);
for k_idx = 1 : n_t
    ws = states{k_idx}.wellSol(1);
    bhp_psia(k_idx)  =  ws.bhp     * fu.Pa_to_psia;
    rate_STBd(k_idx) =  abs(ws.qWs) * fu.m3s_to_STBd;
end
delta_p_psi = pi_psia - bhp_psia;

fprintf('\nWell solution summary:\n');
fprintf('  BHP range  : %.1f - %.1f psia\n', min(bhp_psia), max(bhp_psia));
fprintf('  Rate range : %.1f - %.1f STB/day\n', min(rate_STBd), max(rate_STBd));

%% ── 13. Export ───────────────────────────────────────────────────────
t_export_start = tic;

X_psia = zeros(G.cells.num, n_t);
for k_idx = 1 : n_t
    X_psia(:, k_idx) = states{k_idx}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUT, 'snapshots_psia.mat'), 'X_psia', '-v7.3');
save(fullfile(OUT, 'well_data.mat'), ...
     'bhp_psia', 'delta_p_psi', 'rate_STBd', 'time_hr', '-v7.3');

export_wall_time_s = toc(t_export_start);

%% ── 14. Params JSON (with profiling fields) ───────────────────────────
p = struct();
p.case_name            = case_name;
p.k_mD                 = k_mD;
p.S_skin               = S_skin_val;
p.schedule_id          = schedule_id;
p.schedule_label       = sched_info.label;
p.rates_STBd           = sched_info.rates_STBd(:)';
p.t_change_hr          = sched_info.t_change_hr(:)';
p.phi                  = phi;
p.ct_psi1              = ct_psi1;
p.mu_cp                = mu_cp;
p.Bo                   = Bo;
p.h_ft                 = h_m  * fu.m_to_ft;
p.rw_ft                = rw_m * fu.m_to_ft;
p.C_bbl_psi            = C_bbl_psi;
p.p_init_psia          = pi_psia;
p.kh_mDft              = kh_mDft;
p.re_ft                = re_ft;
p.nx                   = nx;
p.ny                   = ny;
p.dx_ft                = dx * fu.m_to_ft;
p.n_timesteps          = n_t;
p.t_start_hr           = sched_info.t_start_hr;
p.t_end_hr             = sched_info.t_end_hr;
p.t_wbs_end_hr         = t_wbs_hr;
p.t_pss_start_hr       = t_pss_hr;
p.dp_prime_psi         = dp_prime;
p.m_slope_psi_cycle    = m_slope;
p.well_cell_matlab     = wc;
% ── NEW PROFILING FIELDS ──────────────────────────────────────────────
p.sim_wall_time_s      = sim_wall_time_s;
p.setup_wall_time_s    = setup_wall_time_s;
p.export_wall_time_s   = export_wall_time_s;
p.total_wall_time_s    = toc(t_total_start);
p.wall_per_step_ms     = wall_per_step_ms;
p.steps_completed      = n_t;
p.steps_requested      = n_steps;
p.mem_before_MB        = mem_info_before;
p.mem_after_MB         = mem_info_after;
p.mem_peak_delta_MB    = mem_peak_delta;

write_json(fullfile(OUT, 'params.json'), p);

%% ── 15. Final summary ────────────────────────────────────────────────
total_wall_time_s = toc(t_total_start);
fprintf('\n==========================================================\n');
fprintf(' PROFILING SUMMARY : %s\n', case_name);
fprintf('==========================================================\n');
fprintf('  Grid + model setup  : %7.2f s\n', setup_wall_time_s);
fprintf('  MRST simulation     : %7.2f s   (%.0f ms/step)\n', ...
        sim_wall_time_s, wall_per_step_ms);
fprintf('  Data export         : %7.2f s\n', export_wall_time_s);
fprintf('  Total wall time     : %7.2f s   (%.2f min)\n', ...
        total_wall_time_s, total_wall_time_s/60);
fprintf('  Steps completed     : %d / %d\n', n_t, n_steps);
fprintf('  Memory delta (sim)  : %+.1f MB\n', mem_peak_delta);
fprintf('==========================================================\n\n');

end