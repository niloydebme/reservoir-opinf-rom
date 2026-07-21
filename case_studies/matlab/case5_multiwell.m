%==========================================================================
%  case6_multiwell.m  —  MRST 2025b
%
%  Three-well multi-rate production with STAGGERED STARTUP.
%  Wells start one after another and each changes rate independently.
%
%  PHYSICAL RATIONALE FOR STAGGERED STARTUP:
%    Real field development proceeds sequentially — W1 is drilled and
%    put on production first, W2 is brought online once W1 has established
%    a pressure transient (here at t=500 hr, during RF), and W3 is added
%    near the interference onset time.  Each well also changes rate at
%    a different time, creating a rich multi-event history.
%
%  STARTUP SEQUENCE:
%    W1 : t =  0.001 hr  q = 600 STB/day  (first producer)
%    W2 : t =    500 hr  q = 800 STB/day  (W1 in RF)
%    W3 : t =  1,500 hr  q = 1000 STB/day (near t_int = 1321 hr)
%
%  RATE CHANGE SEQUENCE (each well changes independently):
%    W1 changes at t =  3,000 hr : 600 → 400 STB/day  (during interference)
%    W2 changes at t =  6,000 hr : 800 → 500 STB/day  (during interference)
%    W3 changes at t = 10,000 hr : 1000 → 700 STB/day (approaching BDF)
%    All shut-in : t = 20,000 hr  (buildup)
%
%  PERIOD STRUCTURE (7 periods, 420 total log-spaced steps):
%    P1:  0.001 –   500 hr  ( 60 steps) W1=600  W2=0    W3=0
%    P2:    500 –  1500 hr  ( 60 steps) W1=600  W2=800  W3=0
%    P3:  1,500 –  3,000 hr  ( 60 steps) W1=600  W2=800  W3=1000
%    P4:  3,000 –  6,000 hr  ( 60 steps) W1=400  W2=800  W3=1000
%    P5:  6,000 – 10,000 hr  ( 60 steps) W1=400  W2=500  W3=1000
%    P6: 10,000 – 20,000 hr  ( 70 steps) W1=400  W2=500  W3=700
%    P7: 20,000 – 30,000 hr  ( 50 steps) all shut-in (buildup)
%
%  TIMESCALE HIERARCHY:
%    t_WBS =    9.3 hr   wellbore storage end
%    W2 on =  500   hr   W1 in radial flow, clear individual signature
%    W3 on = 1500   hr   near interference onset
%    t_int = 1321   hr   W1-W2 interference onset
%    W1 Dq = 3000   hr   rate change during interference
%    W2 Dq = 6000   hr   rate change during interference
%    W3 Dq = 10000  hr   rate change approaching BDF
%    t_BDF = 13212  hr   boundary domination onset
%    shutin= 20000  hr   full pressure buildup
%    t_end = 30000  hr
%
%  Reservoir 10x6 km, k=20 mD, wells at ix=40,50,60 spacing=1000m
%  Parameters identical to Cases 1, 4, and previous Case 6.
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
fprintf(' Case 6 — Three-Well Staggered Startup + Multi-Rate\n');
fprintf(' W1 starts t=0  |  W2 at t=500hr  |  W3 at t=1500hr\n');
fprintf(' Each well changes rate independently\n');
fprintf('==========================================================\n\n');

SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR, 'utils'));
fu = si_to_field();

OUT = fullfile(SCRIPT_DIR,'..','data','01_synthetic','case6_multiwell');
if ~exist(OUT,'dir'), mkdir(OUT); end
fprintf('Output  : %s\n\n', OUT);

%% ── 1. Fluid and well parameters (identical to Cases 1, 4) ──────────────
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

q_SI_func = @(q_STBd) -(q_STBd * Bo * 0.158987) / 86400;

%% ── 2. Grid: 100×100, dx=100m, dy=60m ───────────────────────────────────
nx = 100;  ny = 100;
dx = 100;  dy = 60;
G  = cartGrid([nx, ny, 1], [nx*dx, ny*dy, h_m]);
G  = computeGeometry(G);
re_x_m  = (nx*dx)/2;
re_x_ft = re_x_m * fu.m_to_ft;

fprintf('Grid    : %dx%d = %d cells\n', nx, ny, G.cells.num);
fprintf('          dx=%.0fm  Lx=%.0fm (%.1fkm)\n', dx, nx*dx, nx*dx/1000);
fprintf('          dy=%.0fm  Ly=%.0fm (%.1fkm)\n', dy, ny*dy, ny*dy/1000);

%% ── 3. Rock ──────────────────────────────────────────────────────────────
k_mD  = 20.0;
k_SI  = k_mD * milli * darcy;
rock  = makeRock(G, k_SI, phi);

%% ── 4. Fluid, model, initial state ──────────────────────────────────────
fluid  = initSimpleADIFluid('phases','W','mu',mu_SI,'rho',rho_SI, ...
                            'c',ct_SI,'pRef',pi_SI);
model  = GenericBlackOilModel(G, rock, fluid, ...
         'water',true,'oil',false,'gas',false);
state0 = initResSol(G, pi_SI, 1.0);
fprintf('Model   : %s\n', class(model));

%% ── 5. Well cells ────────────────────────────────────────────────────────
iy_well  = 50;
ix_wells = [40, 50, 60];
wc_wells = (iy_well-1)*nx + ix_wells;    % [4940, 4950, 4960]

fprintf('\nWell cells:\n');
for wi = 1:3
    fprintf('  W%d : cell=%d  x=%.0fm\n', wi, wc_wells(wi), ...
            (ix_wells(wi)-0.5)*dx);
end
fprintf('  Spacing: %.0fm = %.0fft\n', ...
        (ix_wells(2)-ix_wells(1))*dx, ...
        (ix_wells(2)-ix_wells(1))*dx*fu.m_to_ft);

%% ── 6. Rate schedule ─────────────────────────────────────────────────────
%
%  q_sched(pp, wi) = rate of well wi during period pp [STB/day]
%  0 = shut-in (well present but not producing)
%
%  Staggered startup:
%    P1: W1 only    (W2, W3 shut-in)
%    P2: W1 + W2   (W3 shut-in)
%    P3: W1+W2+W3  all producing (constant rates)
%    P4: W1 rate change (others constant)
%    P5: W2 rate change (others constant)
%    P6: W3 rate change (others constant)
%    P7: all shut-in (buildup)
%
q_sched = [ 800,    -1e-10,    -1e-10;    % P1
            800,  400,    -1e-10;    % P2
            800,  400,  300;    % P3
            600,  400,  300;    % P4  W1 changes
            600,  400,  300;    % P5  W2 changes
            600,  400,  300;    % P6  W3 changes
              -1e-10,    -1e-10,    -1e-10];   % P7  buildup

%  Period time boundaries [t_start, t_end, n_steps]
periods_hr = [0.001,   500,  60;
                500,  1500,  60;
               1500,  3000,  60;
               3000,  6000,  60;
               6000, 10000,  60;
              10000, 20000,  70;
              20000, 30000,  50];

n_periods  = size(periods_hr, 1);

%% ── 7. Build per-period well controls ────────────────────────────────────
W_ctrl = cell(n_periods, 1);

fprintf('\nRate schedule:\n');
fprintf('  %4s  %8s  %8s  %6s  %6s  %6s\n', ...
        'Per', 't_start', 't_end', 'W1', 'W2', 'W3');
fprintf('  %s\n', repmat('-',1,50));

for pp = 1:n_periods
    W_p = [];
    for wi = 1:3
        q_val = q_sched(pp, wi);
        W_p = addWell(W_p, G, rock, wc_wells(wi), ...
            'Type',   'rate', ...
            'Val',    q_SI_func(q_val), ...
            'Radius', rw_m, ...
            'Skin',   S_skin, ...
            'Cs',     C_SI, ...
            'Name',   sprintf('PROD%d', wi), ...
            'Comp_i', 1);
    end
    W_ctrl{pp} = W_p;
    fprintf('  %4d  %8.3f  %8.0f  %6.0f  %6.0f  %6.0f\n', ...
            pp, periods_hr(pp,1), periods_hr(pp,2), ...
            q_sched(pp,1), q_sched(pp,2), q_sched(pp,3));
end

%% ── 8. Build MRST schedule ───────────────────────────────────────────────
dt_all   = [];
ctrl_idx = [];

for pp = 1:n_periods
    t_s  = periods_hr(pp, 1);
    t_e  = periods_hr(pp, 2);
    n_st = periods_hr(pp, 3);
    t_nodes  = logspace(log10(t_s*3600), log10(t_e*3600), n_st+1)';
    dt_p     = diff(t_nodes);
    dt_all   = [dt_all; dt_p];
    ctrl_idx = [ctrl_idx, pp*ones(1, n_st)];
end

schedule.control = struct('W', W_ctrl);
schedule.step.val     = dt_all;
schedule.step.control = ctrl_idx';

fprintf('\nSchedule: %d total steps | t = 0.001 to %.0f hr\n', ...
        numel(dt_all), periods_hr(end,2));

%% ── 9. Analytical timescales ─────────────────────────────────────────────
mu_cp   = mu_SI  * fu.Pas_to_cp;
ct_psi1 = ct_SI  * fu.Pam1_to_psi1;
pi_psia = pi_SI  * fu.Pa_to_psia;
h_ft    = h_m    * fu.m_to_ft;
kh_mDft = k_mD   * h_ft;

t_wbs_hr = 170000 * C_bbl_psi * mu_cp * exp(0.14*S_skin) / kh_mDft;
t_bdf_hr = (phi*mu_cp*ct_psi1*re_x_ft^2)/(0.000264*k_mD)*0.1;
d_well_ft = (ix_wells(2)-ix_wells(1))*dx*fu.m_to_ft;
t_int_hr  = phi*mu_cp*ct_psi1*d_well_ft^2/(4*0.000264*k_mD);

% Bourdet flat for each period
q_tot = sum(q_sched, 2);
dp_flat = 70.6 * q_tot * mu_cp * Bo / kh_mDft;

fprintf('\nTimescale hierarchy:\n');
fprintf('  t_WBS           = %6.1f hr\n', t_wbs_hr);
fprintf('  W2 startup      = %6.0f hr\n', periods_hr(2,1));
fprintf('  W3 startup      = %6.0f hr\n', periods_hr(3,1));
fprintf('  t_int (W1-W2)   = %6.0f hr\n', t_int_hr);
fprintf('  W1 rate change  = %6.0f hr\n', periods_hr(4,1));
fprintf('  W2 rate change  = %6.0f hr\n', periods_hr(5,1));
fprintf('  W3 rate change  = %6.0f hr\n', periods_hr(6,1));
fprintf('  t_BDF           = %6.0f hr\n', t_bdf_hr);
fprintf('  All shut-in     = %6.0f hr\n', periods_hr(7,1));
fprintf('  t_end           = %6.0f hr\n', periods_hr(end,2));
fprintf('\n  t_BDF/t_int = %.0fx (%.2f decades of pure interference)\n', ...
        t_bdf_hr/t_int_hr, log10(t_bdf_hr/t_int_hr));

fprintf('\nBourdet flat per period (total rate):\n');
for pp = 1:n_periods-1
    fprintf('  P%d (qtot=%4.0f STB/d): dp''=%.4f psi\n', ...
            pp, q_tot(pp), dp_flat(pp));
end

%% ── 10. Simulate ─────────────────────────────────────────────────────────
fprintf('\nRunning simulation...\n');
t_wall = tic;
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
n_t = numel(states);
fprintf('  Completed %d / %d timesteps in %.1f s\n', ...
        n_t, numel(dt_all), toc(t_wall));

%% ── 11. Extract well data ────────────────────────────────────────────────
time_hr   = cumsum(dt_all(1:n_t)) * fu.s_to_hr;
bhp_psia  = zeros(3, n_t);
rate_STBd = zeros(3, n_t);

for k_idx = 1:n_t
    ws = states{k_idx}.wellSol;
    for wi = 1:3
        bhp_psia(wi,  k_idx) = ws(wi).bhp     * fu.Pa_to_psia;
        rate_STBd(wi, k_idx) = abs(ws(wi).qWs) * fu.m3s_to_STBd;
    end
end

delta_p_psi  = pi_psia - bhp_psia;
total_rate   = sum(rate_STBd, 1);
cum_prod_STB = cumtrapz(time_hr, total_rate);

fprintf('\nBHP ranges:\n');
for wi = 1:3
    fprintf('  W%d : %.1f – %.1f psia\n', wi, ...
            min(bhp_psia(wi,:)), max(bhp_psia(wi,:)));
end
fprintf('Total cumulative : %.0f STB\n', cum_prod_STB(end));

%% ── 12. Export ───────────────────────────────────────────────────────────
fprintf('\nExporting...\n');

X_psia = zeros(G.cells.num, n_t);
for k_idx = 1:n_t
    X_psia(:, k_idx) = states{k_idx}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUT,'snapshots_psia.mat'),'X_psia','-v7.3');
fprintf('  snapshots_psia.mat  [%d x %d]\n', G.cells.num, n_t);

save(fullfile(OUT,'well_data.mat'), ...
     'bhp_psia','delta_p_psi','rate_STBd','time_hr', ...
     'total_rate','cum_prod_STB','-v7.3');
fprintf('  well_data.mat       [3 x %d]\n', n_t);

perm_mD = rock.perm * fu.m2_to_mD;
poro    = rock.poro;
save(fullfile(OUT,'rock_field.mat'),'perm_mD','poro','-v7.3');
centroids_ft = G.cells.centroids * fu.m_to_ft;
save(fullfile(OUT,'grid_field.mat'),'centroids_ft','-v7.3');
fprintf('  rock_field.mat  grid_field.mat\n');

% params.json
p.description      = 'Case 6: 3-well staggered startup multirate, 10x6km, k=20mD';
p.p_init_psia      = pi_psia;
p.k_mD             = k_mD;
p.k_well_cell_mD   = k_mD;
p.phi              = phi;
p.ct_psi1          = ct_psi1;
p.mu_cp            = mu_cp;
p.Bo               = Bo;
p.h_ft             = h_ft;
p.rw_ft            = rw_m * fu.m_to_ft;
p.S_skin           = S_skin;
p.C_bbl_psi        = C_bbl_psi;
p.kh_mDft          = kh_mDft;
p.re_x_ft          = re_x_ft;
p.nx               = nx;
p.ny               = ny;
p.dx_ft            = dx * fu.m_to_ft;
p.dy_ft            = dy * fu.m_to_ft;
p.n_timesteps      = n_t;
p.wc_W1            = wc_wells(1);
p.wc_W2            = wc_wells(2);
p.wc_W3            = wc_wells(3);
p.d_well_ft        = d_well_ft;
p.t_wbs_end_hr     = t_wbs_hr;
p.t_int_hr         = t_int_hr;
p.t_pss_start_hr   = t_bdf_hr;
% Well startup times
p.t_W1_start_hr    = periods_hr(1,1);
p.t_W2_start_hr    = periods_hr(2,1);
p.t_W3_start_hr    = periods_hr(3,1);
% Rate change times
p.t_W1_change_hr   = periods_hr(4,1);
p.t_W2_change_hr   = periods_hr(5,1);
p.t_W3_change_hr   = periods_hr(6,1);
p.t_shutin_hr      = periods_hr(7,1);
% Rates per period
for pp = 1:n_periods
    p.(sprintf('q_P%d_W1',pp)) = q_sched(pp,1);
    p.(sprintf('q_P%d_W2',pp)) = q_sched(pp,2);
    p.(sprintf('q_P%d_W3',pp)) = q_sched(pp,3);
    p.(sprintf('q_tot_P%d',pp))= q_tot(pp);
    p.(sprintf('dp_prime_P%d_psi',pp)) = dp_flat(pp);
    p.(sprintf('t_period%d_end_hr',pp)) = periods_hr(pp,2);
end

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
fprintf(' Case 6 complete\n');
fprintf(' Startup  : W1@t=0.001  W2@t=500hr  W3@t=1500hr\n');
fprintf(' Changes  : W1@t=3000  W2@t=6000  W3@t=10000hr\n');
fprintf(' Shut-in  : all@t=20000hr  |  t_end=30000hr\n');
fprintf(' t_int=%.0fhr  t_BDF=%.0fhr  gap=%.0fx\n', ...
        t_int_hr, t_bdf_hr, t_bdf_hr/t_int_hr);
fprintf(' Output   : %s\n', OUT);
fprintf('==========================================================\n');