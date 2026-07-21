%=========================================================================
%  CASE 5 — Non-stationary events: rate change + skin damage + infill well
%=========================================================================

mrst_root = 'E:\Datasets Module\MRST2025';
if isempty(which('mrstPath'))
    addpath(mrst_root);
    run(fullfile(mrst_root, 'startup.m'));
end
mrstModule add ad-core ad-blackoil ad-props

clear; clc;
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR,'utils'));
fu = si_to_field();

fprintf('=== Case 5: non-stationary events ===\n');

OUTDIR = fullfile(SCRIPT_DIR,'..','data','01_synthetic','case5_nonstationary');
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

%--- Grid (SI) ----------------------------------------------------------
nx = 50; ny = 50; dz_m = 10; rw_m = 0.10;
G  = cartGrid([nx ny 1], [500, 500, dz_m]);
G  = computeGeometry(G);
fprintf('Grid: %d x %d = %d cells\n', nx, ny, G.cells.num);

%--- Rock (SI) ----------------------------------------------------------
k_SI  = 80 * milli * darcy;
phi   = 0.20;
rock  = makeRock(G, k_SI, phi);

%--- Fluid + model (SI) ------------------------------------------------
p_init_SI = 260 * barsa;
mu_SI     = 1e-3;
ct_SI     = 1e-9;
rho_SI    = 800;
fluid = initSimpleADIFluid('phases','W','mu',mu_SI,'rho',rho_SI, ...
                            'c',ct_SI,'pRef',p_init_SI);
model = GenericBlackOilModel(G, rock, fluid, ...
            'water',true,'oil',false,'gas',false);
state0 = initResSol(G, p_init_SI, 1.0);
fprintf('Model: %s\n', class(model));

%--- Well locations -----------------------------------------------------
wc1 = (round(ny/2)-1)*nx + round(nx/2);   % PROD1 central
wc2 = (round(ny/4)-1)*nx + round(nx/4);   % PROD2 infill
fprintf('PROD1 cell: %d  |  PROD2 cell: %d\n', wc1, wc2);

%--- Unit conversion factor for rates ----------------------------------
STBd_to_SI = 0.158987 / 86400;   % STB/day → m³/s

%=======================================================================
%  PHASE DEFINITIONS
%  KEY FIX: BOTH wells are present in every phase.
%  Inactive well (PROD2 in phases 1-3) has rate = 0 (shut-in).
%  This keeps the well array dimension constant across all controls.
%
%  Phase 1 (0-20 d)  : baseline    PROD1=-100 STB/d skin=0   PROD2=0
%  Phase 2 (20-40 d) : rate change PROD1=-200 STB/d skin=0   PROD2=0
%  Phase 3 (40-60 d) : skin damage PROD1=-200 STB/d skin=10  PROD2=0
%  Phase 4 (60-80 d) : infill      PROD1=-200 STB/d skin=10  PROD2=-80
%  Phase 5 (80-100 d): stable      PROD1=-200 STB/d skin=10  PROD2=-80
%=======================================================================

%  Table columns: [q_PROD1  skin_PROD1  q_PROD2  skin_PROD2]
%  Rates in STB/day (negative = producer, 0 = shut-in)
phase_table = [ ...
  -100,  0,   0,  0;   % Phase 1: baseline
  -200,  0,   0,  0;   % Phase 2: rate change
  -200, 10,   0,  0;   % Phase 3: skin damage
  -200, 10, -80,  0;   % Phase 4: infill well open
  -200, 10, -80,  0];  % Phase 5: stable

phase_labels   = {'baseline','rate_change','skin_damage', ...
                  'infill_well','stable'};
event_times_d  = [0, 20, 40, 60, 80];   % days
event_times_hr = event_times_d * 24;    % hours
n_spp          = 200;                   % steps per phase

schedule.control      = [];
schedule.step.val     = [];
schedule.step.control = [];

for ip = 1:size(phase_table,1)

    T_phase = 20 * day;
    dt_p    = logspace(log10(0.002*day), log10(T_phase), n_spp);
    dt_p    = diff([0; cumsum(dt_p(:))]);

    q1_SI   = phase_table(ip,1) * STBd_to_SI;   % PROD1 rate [m³/s]
    sk1     = phase_table(ip,2);                 % PROD1 skin
    q2_SI   = phase_table(ip,3) * STBd_to_SI;   % PROD2 rate (0 = shut-in)
    sk2     = phase_table(ip,4);                 % PROD2 skin

    % Both wells always present — PROD2 shut-in (rate=0) in phases 1-3
    W_phase = addWell([], G, rock, wc1, ...
                'Type','rate','Val',q1_SI, ...
                'Radius',rw_m,'Skin',sk1, ...
                'Name','PROD1','Comp_i',1);
    W_phase = addWell(W_phase, G, rock, wc2, ...
                'Type','rate','Val',q2_SI, ...
                'Radius',rw_m,'Skin',sk2, ...
                'Name','PROD2','Comp_i',1);

    schedule.control(end+1).W     = W_phase;
    schedule.step.val             = [schedule.step.val; dt_p];
    schedule.step.control         = [schedule.step.control; ...
                                      ip * ones(n_spp,1)];

    fprintf('Phase %d (%s): PROD1=%.0f STBd skin=%.0f | PROD2=%.0f STBd skin=%.0f\n', ...
            ip, phase_labels{ip}, ...
            phase_table(ip,1), phase_table(ip,2), ...
            phase_table(ip,3), phase_table(ip,4));
end
fprintf('Total timesteps: %d  |  Total time: %.0f hours\n', ...
        numel(schedule.step.val), sum(schedule.step.val)*fu.s_to_hr);

%--- Simulate ----------------------------------------------------------
fprintf('\nRunning simulation...\n');
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
wellSols = cellfun(@(s) s.wellSol, states, 'UniformOutput', false);
n_t = numel(states);
n_c = G.cells.num;
fprintf('Simulation complete. %d timesteps.\n', n_t);

%=======================================================================
%  EXPORT — PTA/RTA standard field units
%=======================================================================

%----- (A) Snapshot matrix [psia] --------------------------------------
X_psia = zeros(n_c, n_t);
for k = 1:n_t
    X_psia(:,k) = states{k}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUTDIR,'snapshots_psia.mat'), 'X_psia', '-v7.3');
fprintf('Saved snapshots_psia.mat  [%d x %d]\n', n_c, n_t);

%----- (B) Well data: BHP [psia], rate [STB/day], time [hours] ---------
time_hr    = cumsum(schedule.step.val) * fu.s_to_hr;
p_init_psia = p_init_SI * fu.Pa_to_psia;
n_w         = 2;   % always 2 wells now

bhp_psia   = zeros(n_w, n_t);
rate_STBd  = zeros(n_w, n_t);

for k = 1:n_t
    ws = wellSols{k};
    for w = 1:numel(ws)
        bhp_psia(w,k)  = ws(w).bhp  * fu.Pa_to_psia;
        rate_STBd(w,k) = ws(w).qWs  * fu.m3s_to_STBd;
    end
end

delta_p_psi = p_init_psia - bhp_psia;    % drawdown [psi]
well_names  = {'PROD1','PROD2'};

save(fullfile(OUTDIR,'well_data.mat'), ...
     'bhp_psia','delta_p_psi','rate_STBd','time_hr','well_names','-v7.3');
fprintf('Saved well_data.mat\n');

%----- (C) Event log in HOURS (key for streaming DMD in Python) --------
event_log.times_hr           = event_times_hr;
event_log.times_days         = event_times_d;
event_log.labels             = phase_labels;
event_log.steps_per_phase    = n_spp;
event_log.cumulative_steps   = n_spp * (1:5);
event_log.phase_table_STBd   = phase_table;
event_log.col_labels         = {'q_PROD1_STBd','skin_PROD1', ...
                                 'q_PROD2_STBd','skin_PROD2'};
event_log.time_unit          = 'hours';
event_log.pressure_unit      = 'psia';
event_log.rate_unit          = 'STB/day';
save(fullfile(OUTDIR,'event_log.mat'), 'event_log', '-v7.3');
fprintf('Saved event_log.mat\n');

%----- (D) Rock [mD, fraction] -----------------------------------------
perm_mD = rock.perm * fu.m2_to_mD;
poro    = rock.poro;
save(fullfile(OUTDIR,'rock_field.mat'), 'perm_mD','poro','-v7.3');

%----- (E) Grid geometry [ft] ------------------------------------------
centroids_ft = G.cells.centroids * fu.m_to_ft;
save(fullfile(OUTDIR,'grid_field.mat'), 'centroids_ft','-v7.3');

%----- (F) Units metadata ----------------------------------------------
units.pressure        = 'psia';
units.pressure_change = 'psi';
units.time            = 'hours';
units.rate_liquid     = 'STB/day  (negative=production, 0=shut-in)';
units.permeability    = 'mD';
units.thickness       = 'ft';
units.porosity        = 'fraction (dimensionless)';
units.compressibility = 'psi^-1';
units.viscosity       = 'cp';
units.radius          = 'ft';
units.area            = 'ft^2';
units.skin            = 'dimensionless';
save(fullfile(OUTDIR,'units.mat'), 'units','-v7.3');

%----- (G) Params JSON — all field units --------------------------------
h_ft  = dz_m * fu.m_to_ft;
p.k_mD              = k_SI     * fu.m2_to_mD;
p.phi               = phi;
p.ct_psi1           = ct_SI    * fu.Pam1_to_psi1;
p.mu_cp             = mu_SI    * fu.Pas_to_cp;
p.rho_lbft3         = rho_SI   * fu.kgm3_to_lbft3;
p.rw_ft             = rw_m     * fu.m_to_ft;
p.h_ft              = h_ft;
p.p_init_psia       = p_init_psia;
p.kh_mDft           = (k_SI * fu.m2_to_mD) * h_ft;
p.domain_ft         = 500      * fu.m_to_ft;
p.n_phases          = 5;
p.n_steps_per_phase = n_spp;
p.event_times_hr    = '0,480,960,1440,1920';
p.well1_cell        = wc1;
p.well2_cell        = wc2;
p.well2_active_from_phase = 4;
p.case_name         = 'nonstationary_events';
p.unit_system       = 'PTA_RTA_field_units';
p.fix_note          = 'both_wells_always_defined_PROD2_shutin_phases_1to3';

fid = fopen(fullfile(OUTDIR,'params.json'),'w');
fprintf(fid,'{\n');
fn = fieldnames(p);
for i = 1:numel(fn)
    v = p.(fn{i});
    if ischar(v), fprintf(fid,'  "%s": "%s"',fn{i},v);
    else,         fprintf(fid,'  "%s": %g',  fn{i},v); end
    if i < numel(fn), fprintf(fid,','); end
    fprintf(fid,'\n');
end
fprintf(fid,'}\n');
fclose(fid);

%--- Summary -----------------------------------------------------------
fprintf('\n=== Case 5 complete ===\n');
fprintf('Output folder  : %s\n', OUTDIR);
fprintf('Snapshots      : [%d cells x %d steps] in psia\n', n_c, n_t);
fprintf('p_init         : %.1f psia\n', p_init_psia);
fprintf('PROD1 BHP      : %.1f - %.1f psia\n', ...
        min(bhp_psia(1,:)), max(bhp_psia(1,:)));
fprintf('PROD2 BHP      : %.1f - %.1f psia\n', ...
        min(bhp_psia(2,:)), max(bhp_psia(2,:)));
fprintf('Time range     : %.5f - %.2f hours\n', min(time_hr), max(time_hr));
fprintf('k = %.1f mD  |  kh = %.1f mD-ft  |  ct = %.4e psi^-1\n', ...
        p.k_mD, p.kh_mDft, p.ct_psi1);
fprintf('\nEvent times (hours): 0 | 480 | 960 | 1440 | 1920\n');
fprintf('%-8s %-15s %10s %10s %10s %10s\n', ...
        'Phase','Label','PROD1 q','PROD1 S','PROD2 q','PROD2 S');
for ip = 1:5
    fprintf('%-8d %-15s %8.0f %10.0f %10.0f %10.0f\n', ...
            ip, phase_labels{ip}, phase_table(ip,1), phase_table(ip,2), ...
            phase_table(ip,3), phase_table(ip,4));
end