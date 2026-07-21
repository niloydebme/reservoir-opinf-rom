function info = make_rate_schedule(schedule_id, t_end_hr, n_steps)
%MAKE_RATE_SCHEDULE  Time discretisation and rate-period info for a given ID.
%
%  info = make_rate_schedule(schedule_id)
%  info = make_rate_schedule(schedule_id, t_end_hr, n_steps)
%
%  Schedule definitions (rates in STB/day, change times in hours):
%
%  ── Training schedules (IDs 1–7, unchanged) ──────────────────────────
%    1  Constant  800  full duration
%    2  Constant  400  full duration
%    3  Constant 1200  full duration
%    4  1000 → 500     change at t = 1000 hr
%    5   400 → 900     change at t = 1000 hr
%    6   800 → 0       change at t = 2500 hr  (buildup / shut-in)
%    7  600 → 1000 → 300  changes at t = 1667 hr and t = 3333 hr
%
%  ── Validation schedules (IDs 8–9, never used in training) ───────────
%    8   600 → 900 → 0     changes at t = 800 hr and t = 3000 hr
%        (ramp-UP then build-up — ordering never seen during training)
%    9   900 → 600 → 900 → 400   changes at t = 400, 1200, 3500 hr
%        (4-period non-monotonic — training max is 3 periods;
%         change times do not coincide with any training schedule)
%
%  OUTPUT struct fields:
%    .t_nodes_s   [n_steps × 1]       cumulative elapsed times [s]
%    .dt_vec_s    [n_steps × 1]       per-step timestep sizes  [s]
%    .n_controls  scalar              number of distinct rate periods
%    .rates_STBd  [n_controls × 1]   surface production rate per period [STB/day]
%    .t_change_hr [n_controls-1 × 1] elapsed times at which rate changes [hr]
%    .ctrl_vec    [n_steps × 1]       1-based MRST control index per timestep
%    .label       string              human-readable tag used in folder naming
%    .t_start_hr  scalar              first output time [hr]
%    .t_end_hr    scalar              final output time [hr]
%    .n_steps     scalar              number of output times

if nargin < 2 || isempty(t_end_hr), t_end_hr = 5000; end
if nargin < 3 || isempty(n_steps),  n_steps  = 350;  end

t_start_hr = 0.001;

% Log-spaced cumulative times (same strategy as case1_homogeneous.m)
t_nodes_s  = logspace(log10(t_start_hr * 3600), ...
                       log10(t_end_hr   * 3600), n_steps)';
dt_vec_s   = diff([0; t_nodes_s]);
t_hr_nodes = t_nodes_s / 3600;   % used below for control assignment

switch schedule_id

    %% ── Training schedules (unchanged) ──────────────────────────────
    case 1
        rates     = 800;
        t_changes = [];
        label     = 'q800c';
    case 2
        rates     = 400;
        t_changes = [];
        label     = 'q400c';
    case 3
        rates     = 1200;
        t_changes = [];
        label     = 'q1200c';
    case 4
        rates     = [1000; 500];
        t_changes = 1000;
        label     = 'q1000to500';
    case 5
        rates     = [400; 900];
        t_changes = 1000;
        label     = 'q400to900';
    case 6
        rates     = [800; 1e-7];
        t_changes = 2500;
        label     = 'q800buildup';
    case 7
        rates     = [600; 1000; 300];
        t_changes = [1667; 3333];
        label     = 'q600_1000_300';

    %% ── Validation schedules (new) ───────────────────────────────────
    case 8
        % Ramp-UP then build-up.
        % Training has step-DOWN→BU (sched 6) and low→high→low (sched 7).
        % This ordering — low → high → shut-in — was never seen during
        % training.  Change times (800, 3000 hr) differ from all training
        % schedules.
        rates     = [600; 900; 1e-7];
        t_changes = [800; 3000];
        label     = 'q600_900_bu';

    case 9
        % 4 control periods with non-monotonic (up-down-up-down) rates.
        % Training max is 3 periods.  Change times (400, 1200, 3500 hr)
        % do not coincide with any training schedule change time.
        rates     = [900; 600; 900; 400];
        t_changes = [400; 1200; 3500];
        label     = 'q900_600_900_400';

    otherwise
        error('make_rate_schedule: unknown schedule_id %d. Valid IDs: 1–9.', ...
              schedule_id);
end

% Build MRST control index vector:
% ctrl_vec(j) = k means timestep j uses schedule.control(k).W
n_controls = numel(rates);
ctrl_vec   = ones(n_steps, 1);
for ic = 1 : numel(t_changes)
    ctrl_vec(t_hr_nodes > t_changes(ic)) = ic + 1;
end

info.t_nodes_s   = t_nodes_s;
info.dt_vec_s    = dt_vec_s;
info.n_controls  = n_controls;
info.rates_STBd  = rates(:);
info.t_change_hr = t_changes(:);
info.ctrl_vec    = ctrl_vec;
info.label       = label;
info.t_start_hr  = t_start_hr;
info.t_end_hr    = t_end_hr;
info.n_steps     = n_steps;

end
