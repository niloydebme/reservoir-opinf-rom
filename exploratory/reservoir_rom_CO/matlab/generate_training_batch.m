%GENERATE_TRAINING_BATCH  Run all (k, S, schedule) training simulations.
%
%  Parameter grid:
%    k     = [5, 10, 20, 50, 100, 200, 500]  mD   (7 log-spaced values)
%    S     = [0, 2, 5, 8, 10, 12, 15]              (7 linear values)
%    sched = 1 : 7                                  (7 rate schedules)
%
%  Total : 7 × 7 × 7 = 343 simulations  ≈  17 hr at 3 min/sim
%
%  RESUMABLE: already-completed cases (found in manifest.csv) are skipped
%  automatically, so you can stop and restart without re-running anything.
%
%  Manifest file:  data/01_training/manifest.csv
%    Columns: case_name, k_mD, S_skin, schedule_id, schedule_label,
%             n_timesteps, sim_wall_s, status

clear; clc;

%% ── Configuration ────────────────────────────────────────────────────
MRST_ROOT  = 'E:\Datasets Module\MRST2025';
SCRIPT_DIR = fileparts(mfilename('fullpath'));
OUT_ROOT   = fullfile(SCRIPT_DIR, '..', 'data', '01_training_CO');

if ~exist(OUT_ROOT, 'dir'), mkdir(OUT_ROOT); end

addpath(fullfile(SCRIPT_DIR, 'utils'));

%% ── Parameter grid ───────────────────────────────────────────────────
k_grid    = [5, 10, 20, 50, 100, 200, 500];   % mD
S_grid    = [0,  2,  5,  8,  10,  12,  15];   % dimensionless
sched_ids = 1 : 7;

n_k     = numel(k_grid);
n_S     = numel(S_grid);
n_sched = numel(sched_ids);
n_total = n_k * n_S * n_sched;

fprintf('Training batch\n');
fprintf('  k grid    : [%s] mD  (%d values)\n', ...
        num2str(k_grid, '%.0f '), n_k);
fprintf('  S grid    : [%s]     (%d values)\n', ...
        num2str(S_grid, '%.0f '), n_S);
fprintf('  Schedules : %d to %d  (%d types)\n', ...
        sched_ids(1), sched_ids(end), n_sched);
fprintf('  Total     : %d simulations\n\n', n_total);
fprintf('  Output    : %s\n\n', OUT_ROOT);

%% ── Initialise MRST once ─────────────────────────────────────────────
if isempty(which('mrstPath'))
    addpath(MRST_ROOT);
    run(fullfile(MRST_ROOT, 'startup.m'));
end
mrstModule add ad-core ad-blackoil ad-props

%% ── Manifest: create header if file does not exist ───────────────────
manifest_path = fullfile(OUT_ROOT, 'manifest.csv');
if ~exist(manifest_path, 'file')
    fid = fopen(manifest_path, 'w');
    fprintf(fid, 'case_name,k_mD,S_skin,schedule_id,schedule_label,n_timesteps,sim_wall_s,status\n');
    fclose(fid);
    completed = {};
else
    % Read already-completed case names
    T = readtable(manifest_path, 'Delimiter', ',');
    if ~isempty(T) && ismember('status', T.Properties.VariableNames)
        ok_rows  = strcmp(T.status, 'ok');
        completed = T.case_name(ok_rows);
    else
        completed = {};
    end
end
fprintf('Already completed: %d / %d cases\n\n', numel(completed), n_total);

%% ── Main loop ────────────────────────────────────────────────────────
counter   = 0;
n_skipped = 0;
n_errors  = 0;
t_batch   = tic;

for ik = 1 : n_k
    for iS = 1 : n_S
        for isched = 1 : n_sched

            k_val = k_grid(ik);
            S_val = S_grid(iS);
            s_id  = sched_ids(isched);

            counter   = counter + 1;
            case_name = sprintf('k%03d_s%02d_sched%d', ...
                                round(k_val), round(S_val), s_id);

            % ── Skip if already done ──────────────────────────────────
            if any(strcmp(completed, case_name))
                n_skipped = n_skipped + 1;
                fprintf('[%3d/%3d] SKIP  %s\n', counter, n_total, case_name);
                continue
            end

            fprintf('[%3d/%3d] START %s\n', counter, n_total, case_name);

            % ── Run simulation ────────────────────────────────────────
            t_case  = tic;
            status  = 'ok';
            n_t_out = NaN;

            try
                run_single_case(k_val, S_val, s_id, OUT_ROOT, MRST_ROOT);

                % Read back n_timesteps as a cross-check
                wd      = load(fullfile(OUT_ROOT, case_name, 'well_data.mat'), 'time_hr');
                n_t_out = numel(wd.time_hr);

            catch ME
                status   = 'error';
                n_errors = n_errors + 1;
                warning('BATCH:failed', '[%s] FAILED: %s', case_name, ME.message);
            end

            sim_s = toc(t_case);

            % ── Append one line to manifest ───────────────────────────
            si    = make_rate_schedule(s_id);
            fid   = fopen(manifest_path, 'a');
            fprintf(fid, '%s,%.1f,%.0f,%d,%s,%d,%.1f,%s\n', ...
                    case_name, k_val, S_val, s_id, si.label, ...
                    n_t_out, sim_s, status);
            fclose(fid);

            elapsed_batch = toc(t_batch);
            n_done_now    = counter - n_skipped - n_errors;
            if n_done_now > 0
                avg_s     = elapsed_batch / n_done_now;
                remaining = (n_total - counter) * avg_s;
                fprintf('       DONE %.1f s | avg %.1f s/sim | ETA %.0f min\n', ...
                        sim_s, avg_s, remaining/60);
            end

        end
    end
end

%% ── Summary ──────────────────────────────────────────────────────────
fprintf('\n==========================================================\n');
fprintf(' Batch complete\n');
fprintf('  Total   : %d  |  New: %d  |  Skipped: %d  |  Errors: %d\n', ...
        counter, counter - n_skipped - n_errors, n_skipped, n_errors);
fprintf('  Wall time: %.1f min\n', toc(t_batch)/60);
fprintf('  Manifest : %s\n', manifest_path);
fprintf('==========================================================\n');