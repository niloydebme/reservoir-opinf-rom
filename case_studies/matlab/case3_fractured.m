%==========================================================================
%  case3_fractured.m  —  MRST 2025b
%
%  Hydraulically fractured well in a LARGE reservoir, single-phase,
%  constant-rate drawdown.  Non-uniform (logarithmic) grid in x so that
%  the near-fracture region is resolved accurately while the far-field
%  drainage area is large enough for a clear pseudo-radial plateau.
%
%  WHY non-uniform grid and t_end=50000hr:
%    The x-direction uses logspace cell widths: finest near the fracture
%    (capturing the LF-to-PSR transition accurately) and coarsest at the
%    domain boundary.  The y-direction is uniform (fracture is a column).
%    k_m=5 mD (4x tighter than Cases 1&2) + large domain (6000x6000m)
%    gives t_BDF=19025hr — well beyond the 5000hr window of Cases 1&2.
%    t_end=50000hr gives 0.42 log-decades of clear BDF signature.
%
%  Grid:
%    nx = 101 (odd, symmetric about fracture, logspace x-spacing)
%    ny = 120 (uniform y-spacing)
%    Lx = Ly = 6000 m  (square domain — large in BOTH directions)
%    xf = 400 m = 1312 ft  (fracture half-length, spans ny_frac rows)
%
%  Flow regime windows (analytical):
%    WBS   : 0.001  –    0.68 hr
%    LF    : 0.68   – 1429    hr  (3.32 log-decades, Bourdet slope 1/2)
%    PSR   : 1429   – 19025   hr  (1.12 log-decades, Bourdet flat)
%    BDF   : 19025  – 50000   hr  (0.42 log-decades, Bourdet rising)
%
%  Fluid/well parameters identical to Cases 1 & 2:
%    phi=0.18  ct=1.2e-5 psi-1  mu=1.2 cp  Bo=1.15  rho=780 kg/m3
%    h=150ft   rw=0.35ft        S=0  (fracture bypasses skin)
%    C=0.003 bbl/psi  q=800 STB/day  pi=4500 psia
%
%  Output: data/01_synthetic/case3_fractured/
%    snapshots_psia.mat, well_data.mat, rock_field.mat,
%    grid_field.mat, frac_cells.mat, params.json
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
fprintf(' Case 3 — Hydraulically Fractured Well (Large Reservoir)\n');
fprintf(' Non-uniform grid  |  6000x6000 m  |  xf=400 m  |  S=0\n');
fprintf(' WBS -> LF (3.32 dec) -> PSR (1.12 dec) -> BDF\n');
fprintf('==========================================================\n\n');

SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(SCRIPT_DIR, 'utils'));
fu = si_to_field();

OUT = fullfile(SCRIPT_DIR,'..','data','01_synthetic','case3_fractured');
if ~exist(OUT,'dir'), mkdir(OUT); end
fprintf('Output  : %s\n\n', OUT);

%% ── 1. Fluid and well parameters (identical to Cases 1 & 2) ─────────────
phi       = 0.18;
ct_SI     = 1.740e-9;              % Pa-1  (= 1.2e-5 psi-1)
mu_SI     = 1.2e-3;                % 1.2 cp
Bo        = 1.15;
rho_SI    = 780;                   % kg/m3
pi_SI     = 310.3 * barsa;         % ~4500 psia
h_m       = 45.72;                 % m  (= 150 ft)
rw_m      = 0.1067;                % m  (= 0.35 ft)
S_skin    = 0.0;                   % no skin
C_bbl_psi = 0.003;
C_SI      = C_bbl_psi * 0.158987 / 6894.76;   % 6.918e-8 m3/Pa
q_STBd    = 800;
q_SI      = -(q_STBd * Bo * 0.158987) / 86400;

%% ── 2. Non-uniform grid ──────────────────────────────────────────────────
%
%  x-direction (perpendicular to fracture):
%    Logspace spacing, finest at x=Lx/2 (fracture) coarsest at boundaries.
%    nx=101 (odd) so there is exactly one column at x=Lx/2.
%
%  y-direction (along fracture):
%    Uniform spacing, ny=120 rows.
%
%  Domain: Lx=Ly=6000 m  (square — large in BOTH directions)
%
nx   = 101;  ny   = 120;
Lx   = 6000; Ly   = 6000;
h_ax = h_m;

% x-node vector: symmetric logspace about x=Lx/2
half_nx   = (nx - 1) / 2;                         % 50 half-steps
x_half    = logspace(log10(0.25), log10(Lx/2), half_nx + 1);  % [0.25 .. 3000] m
x_full    = [Lx/2 - fliplr(x_half(2:end)), x_half + Lx/2];   % [0 .. 6000] m
if numel(x_full) ~= nx + 1
    x_full = linspace(0, Lx, nx + 1);             % fallback: uniform
    warning('tensorGrid x fallback to uniform');
end

% y-node vector: uniform
y_full = linspace(0, Ly, ny + 1);

G  = tensorGrid(x_full, y_full, [0, h_ax]);
G  = computeGeometry(G);

re_m   = Lx / 2;                    % drainage radius (x-direction)
re_ft  = re_m * fu.m_to_ft;
xf_m   = 400;                       % fracture half-length
xf_ft  = xf_m * fu.m_to_ft;

fprintf('Grid    : %dx%d = %d cells  (non-uniform x)\n', nx, ny, G.cells.num);
fprintf('          Lx=%.0fm (%.0fft)  Ly=%.0fm (%.0fft)\n', ...
        Lx, Lx*fu.m_to_ft, Ly, Ly*fu.m_to_ft);
fprintf('          re=%.0fm=%.0fft  xf=%.0fm=%.0fft\n', re_m,re_ft,xf_m,xf_ft);

% Minimum and maximum dx (for reference)
dx_min = min(diff(x_full));
dx_max = max(diff(x_full));
fprintf('          dx_min=%.3fm  dx_max=%.1fm\n', dx_min, dx_max);

%% ── 3. Permeability field — matrix + fracture ────────────────────────────
%
%  Fracture: all cells whose x-centroid ≈ Lx/2  AND  y-centroid within xf.
%  Coordinate-based search is required for non-uniform grids.
%
k_matrix_mD   = 5.0;
k_fracture_mD = 100000.0;
k_matrix_SI   = k_matrix_mD  * milli * darcy;
k_fracture_SI = k_fracture_mD * milli * darcy;

cx = G.cells.centroids(:,1);
cy = G.cells.centroids(:,2);

% Fracture cells: x ≈ Lx/2 (within half the finest dx) and y within xf
x_tol      = dx_min * 2;
frac_cells = find( abs(cx - Lx/2)       < x_tol & ...
                   abs(cy - Ly/2)       <= xf_m );

if isempty(frac_cells)
    error('No fracture cells found — check x_tol or xf_m definition');
end

k_field_SI = k_matrix_SI * ones(G.cells.num, 1);
k_field_SI(frac_cells) = k_fracture_SI;
rock = makeRock(G, k_field_SI, phi);

fprintf('\nFracture cells: %d  (|x-Lx/2|<%.3fm  |y-Ly/2|<=%.0fm)\n', ...
        numel(frac_cells), x_tol, xf_m);
fprintf('  k_matrix   = %.1f mD\n', k_matrix_mD);
fprintf('  k_fracture = %.0f mD  (k_f/k_m = %.0f)\n', ...
        k_fracture_mD, k_fracture_mD/k_matrix_mD);

%% ── 4. Fluid, model, initial conditions ─────────────────────────────────
fluid  = initSimpleADIFluid('phases','W','mu',mu_SI,'rho',rho_SI, ...
                            'c',ct_SI,'pRef',pi_SI);
model  = GenericBlackOilModel(G, rock, fluid, ...
         'water',true,'oil',false,'gas',false);
state0 = initResSol(G, pi_SI, 1.0);
fprintf('\nModel   : %s\n', class(model));

%% ── 5. Well — cell closest to (Lx/2, Ly/2), rate + WBS, S=0 ─────────────
[~, wc] = min( (cx - Lx/2).^2 + (cy - Ly/2).^2 );

W = addWell([], G, rock, wc,      ...
    'Type',   'rate',              ...
    'Val',    q_SI,                ...
    'Radius', rw_m,                ...
    'Skin',   S_skin,              ...
    'Cs',     C_SI,                ...
    'Name',   'FRAC_PROD',         ...
    'Comp_i', 1);

fprintf('Well    : cell %d  (x=%.2fm, y=%.2fm)\n', wc, cx(wc), cy(wc));
fprintf('          q=%d STB/day | S=%.0f | C=%.3f bbl/psi\n', ...
        q_STBd, S_skin, C_bbl_psi);

% Well cell size (for params.json / Peaceman reference)
% Find neighbouring cell in x (wc±1 in MATLAB ordering)
wc_ix = mod(wc-1, nx) + 1;   % 1-indexed x-position of well cell
if wc_ix < nx
    dx_wc = abs(cx(wc+1) - cx(wc)) * 2;
else
    dx_wc = abs(cx(wc) - cx(wc-1)) * 2;
end
dy_wc = Ly / ny;              % uniform y-spacing

fprintf('          dx_wc=%.4fm  dy_wc=%.4fm\n', dx_wc, dy_wc);

%% ── 6. Schedule : 400 log-spaced steps, 0.001 to 50000 hr ───────────────
t_start_hr = 0.001;
t_end_hr   = 50000;
n_steps    = 400;
t_nodes_s  = logspace(log10(t_start_hr*3600), ...
                       log10(t_end_hr  *3600), n_steps)';
dt_vec_s   = diff([0; t_nodes_s]);
schedule   = simpleSchedule(dt_vec_s, 'W', W);

fprintf('\nSchedule: %d steps | t = %.4f to %.0f hr\n', ...
        n_steps, t_start_hr, t_end_hr);

%% ── 7. Analytical regime boundaries ─────────────────────────────────────
mu_cp   = mu_SI  * fu.Pas_to_cp;
ct_psi1 = ct_SI  * fu.Pam1_to_psi1;
pi_psia = pi_SI  * fu.Pa_to_psia;
h_ft    = h_m    * fu.m_to_ft;
kh_m    = k_matrix_mD * h_ft;

t_wbs_hr = 170000 * C_bbl_psi* mu_cp * exp(0.14*S_skin) / kh_m;

m_LF = 4.064 * q_STBd * Bo * mu_cp / ...
       (xf_ft * h_ft * sqrt(k_matrix_mD * phi * mu_cp * ct_psi1));

t_psr_hr     = 1600 * phi * mu_cp * ct_psi1 * xf_ft^2 / k_matrix_mD;
dp_prime_psr = 70.6 * q_STBd * mu_cp * Bo / kh_m;
t_bdf_hr     = (phi * mu_cp * ct_psi1 * re_ft^2) / ...
               (0.000264 * k_matrix_mD) * 0.1;
m_slope_psr  = 162.6 * q_STBd * mu_cp * Bo / kh_m;

fprintf('\nFlow regime boundaries:\n');
fprintf('  WBS ends    : t ~ %.3f hr\n',  t_wbs_hr);
fprintf('  LF window   : %.3f – %.0f hr  (%.2f decades)  slope=1/2\n', ...
        t_wbs_hr, t_psr_hr, log10(t_psr_hr/t_wbs_hr));
fprintf('  PSR onset   : t ~ %.0f hr   dp''=%.4f psi\n', t_psr_hr, dp_prime_psr);
fprintf('  BDF onset   : t ~ %.0f hr\n', t_bdf_hr);
fprintf('  PSR window  : %.2f log-decades\n', log10(t_bdf_hr/t_psr_hr));
fprintf('  BDF window  : %.2f log-decades\n', log10(t_end_hr/t_bdf_hr));
fprintf('  m_LF        : %.4f psi / hr^0.5\n', m_LF);

%% ── 8. Simulate ──────────────────────────────────────────────────────────
fprintf('\nRunning simulation...\n');
t_wall = tic;
[~, states, ~] = simulateScheduleAD(state0, model, schedule);
n_t = numel(states);
fprintf('  Completed %d / %d timesteps in %.1f s\n', ...
        n_t, n_steps, toc(t_wall));

%% ── 9. Extract well solution ─────────────────────────────────────────────
time_hr   = cumsum(dt_vec_s(1:n_t)) * fu.s_to_hr;
bhp_psia  = zeros(1, n_t);
rate_STBd = zeros(1, n_t);
for k_idx = 1:n_t
    ws = states{k_idx}.wellSol(1);
    bhp_psia(k_idx)  = ws.bhp      * fu.Pa_to_psia;
    rate_STBd(k_idx) = abs(ws.qWs) * fu.m3s_to_STBd;
end
delta_p_psi = pi_psia - bhp_psia;

fprintf('\nBHP     : %.1f – %.1f psia\n',  min(bhp_psia),    max(bhp_psia));
fprintf('Delta-p : %.4f – %.2f psi\n',   min(delta_p_psi), max(delta_p_psi));
n_wbs=sum(time_hr<t_wbs_hr); n_lf=sum(time_hr>=t_wbs_hr&time_hr<t_psr_hr);
n_psr=sum(time_hr>=t_psr_hr&time_hr<t_bdf_hr); n_bdf=sum(time_hr>=t_bdf_hr);
fprintf('WBS(%d)  LF(%d)  PSR(%d)  BDF(%d)\n', n_wbs,n_lf,n_psr,n_bdf);

%% ── 10. Export ───────────────────────────────────────────────────────────
fprintf('\nExporting...\n');

X_psia = zeros(G.cells.num, n_t);
for k_idx = 1:n_t
    X_psia(:, k_idx) = states{k_idx}.pressure * fu.Pa_to_psia;
end
save(fullfile(OUT,'snapshots_psia.mat'),'X_psia','-v7.3');
fprintf('  snapshots_psia.mat  [%d x %d] psia\n', G.cells.num, n_t);

save(fullfile(OUT,'well_data.mat'), ...
     'bhp_psia','delta_p_psi','rate_STBd','time_hr','-v7.3');
fprintf('  well_data.mat       [%d timesteps]\n', n_t);

perm_mD = rock.perm * fu.m2_to_mD;
poro    = rock.poro;
k_well_cell_mD = perm_mD(wc);
save(fullfile(OUT,'rock_field.mat'),'perm_mD','poro','-v7.3');

% Grid centroids in ft (both x and y columns)
centroids_ft = G.cells.centroids(:,1:2) * fu.m_to_ft;
save(fullfile(OUT,'grid_field.mat'),'centroids_ft','-v7.3');

% Fracture cell indices (MATLAB 1-indexed)
save(fullfile(OUT,'frac_cells.mat'),'frac_cells','-v7.3');
fprintf('  frac_cells.mat      [%d cells]\n', numel(frac_cells));

% params.json — all keys consistent with Cases 1 & 2
p.description        = 'Case 3: fractured well, non-uniform grid, WBS+LF+PSR+BDF';
p.p_init_psia        = pi_psia;
p.k_mD               = k_matrix_mD;
p.k_matrix_mD        = k_matrix_mD;
p.k_fracture_mD      = k_fracture_mD;
p.k_well_cell_mD     = k_well_cell_mD;
p.dx_wc_ft           = dx_wc * fu.m_to_ft;    % well-cell x-width [ft]
p.dy_wc_ft           = dy_wc * fu.m_to_ft;    % well-cell y-width [ft]
p.phi                = phi;
p.ct_psi1            = ct_psi1;
p.mu_cp              = mu_cp;
p.Bo                 = Bo;
p.h_ft               = h_ft;
p.rw_ft              = rw_m * fu.m_to_ft;
p.S_skin             = S_skin;
p.C_bbl_psi          = C_bbl_psi;
p.q_STBd             = q_STBd;
p.kh_mDft            = kh_m;
p.re_ft              = re_ft;
p.xf_ft              = xf_ft;
p.xf_m               = xf_m;
p.nx                 = nx;
p.ny                 = ny;
p.Lx_ft              = Lx * fu.m_to_ft;
p.Ly_ft              = Ly * fu.m_to_ft;
p.n_frac_cells       = numel(frac_cells);
p.n_timesteps        = n_t;
p.t_start_hr         = t_nodes_s(1)   / 3600;
p.t_end_hr           = t_nodes_s(end) / 3600;
p.t_wbs_end_hr       = t_wbs_hr;
p.t_lf_end_hr        = t_psr_hr;
p.t_psr_start_hr     = t_psr_hr;
p.t_pss_start_hr     = t_bdf_hr;
p.dp_prime_psi       = dp_prime_psr;
p.m_LF_psi_hr05      = m_LF;
p.m_slope_psi_cycle  = m_slope_psr;
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
fprintf(' Case 3 complete\n');
fprintf(' Grid    : %dx%d = %d cells (non-uniform x)\n', nx, ny, G.cells.num);
fprintf(' Domain  : %.0fx%.0f m  (%.0fx%.0f ft)\n', Lx,Ly,Lx*fu.m_to_ft,Ly*fu.m_to_ft);
fprintf(' Fracture: k_m=%.0f mD  k_f=%.0f mD  xf=%.0f ft  S=%.0f\n', ...
        k_matrix_mD, k_fracture_mD, xf_ft, S_skin);
fprintf(' WBS : 0.001 – %.3f hr  (%d steps)\n', t_wbs_hr, n_wbs);
fprintf(' LF  : %.3f – %.0f hr   (%d steps)\n', t_wbs_hr, t_psr_hr, n_lf);
fprintf(' PSR : %.0f – %.0f hr   (%d steps)\n', t_psr_hr, t_bdf_hr, n_psr);
fprintf(' BDF : %.0f – %.0f hr   (%d steps)\n', t_bdf_hr, t_end_hr, n_bdf);
fprintf(' Output: %s\n', OUT);
fprintf('==========================================================\n');