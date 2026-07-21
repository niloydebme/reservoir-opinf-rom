%GENERATE_VALIDATION_BATCH  Run validation simulations with two NEW rate
%  schedules (IDs 8 and 9) at unseen (k, S) parameter values.
%
%  HOW IT DIFFERS FROM generate_training_batch.m
%  ──────────────────────────────────────────────
%  1. Output goes to  data/02_validation/  (not 01_training)
%  2. k and S grids contain values OUTSIDE the training grid
%  3. Schedule IDs 8 and 9 are used — novel rate histories never seen
%     during operator inference
%  4. Everything else is identical: calls run_single_case.m unchanged
%
%  VALIDATION PARAMETER GRID
%  ──────────────────────────
%  k = [15, 30, 75, 150, 250] mD   — all off training grid {10,20,50,100,200,300}
%  S = [0, 3, 5]                   — S=0 on training grid (sanity check)
%                                    S=3 between training {2,4}
%                                    S=5 between training {4,6}

% New validation for inv char
%k_grid    = [ 25, 45, 85, 120, 170, 220, 270];   % mD  -- all off training grid
%S_grid    = [0, 1.5, 2.5, 4.5, 5.5, 6.5];                % --  -- mix of seen and unseen
%sched_ids = [2, 3, 5, 9];                   % new validation schedules only

%
%  NEW SCHEDULES  (defined in utils/make_rate_schedule.m, cases 8 and 9)
%  ID 8  q600_1100_bu       600 -> 1100 -> 0  STB/d   (changes at 800, 3000 hr)
%  ID 9  q900_300_1000_200  900 -> 300 -> 1000 -> 200  (changes at 400, 1200, 3500 hr)
%
%  Total : 5 x 3 x 2 = 30 simulations  ~1.5 hr at 3 min/sim
%
%  RESUMABLE: completed cases in val_manifest.csv are skipped automatically.
%
%  Manifest:  data/02_validation/val_manifest.csv

clear; clc;

%% -- Configuration ---------------------------------------------------
MRST_ROOT  = 'E:\Datasets Module\MRST2025';
SCRIPT_DIR = fileparts(mfilename('fullpath'));
OUT_ROOT   = fullfile(SCRIPT_DIR, '..', 'data', '02_validation');

if ~exist(OUT_ROOT, 'dir'), mkdir(OUT_ROOT); end

addpath(fullfile(SCRIPT_DIR, 'utils'));

%% -- Parameter grid --------------------------------------------------
k_grid    = [85, 130, 170, 220, 280];   % mD  -- all off training grid
S_grid    = [3, 5, 7];                % --  -- mix of seen and unseen
sched_ids = [4, 6, 9];                   % new validation schedules only

n_k     = numel(k_grid);
n_S     = numel(S_grid);
n_sched = numel(sched_ids);
n_total = n_k * n_S * n_sched;

fprintf('Validation batch\n');
fprintf('  k grid    : [%s] mD  (%d values)\n', num2str(k_grid, '%.0f '), n_k);
fprintf('  S grid    : [%s]     (%d values)\n', num2str(S_grid, '%.0f '), n_S);
fprintf('  Schedules : IDs [%s]  (%d types)\n', num2str(sched_ids), n_sched);
fprintf('  Total     : %d simulations\n\n', n_total);
fprintf('  Output    : %s\n\n', OUT_ROOT);

%% -- Initialise MRST once --------------------------------------------
if isempty(which('mrstPath'))
    addpath(MRST_ROOT);
    run(fullfile(MRST_ROOT, 'startup.m'));
end
mrstModule add ad-core ad-blackoil ad-props

%% -- Manifest --------------------------------------------------------
manifest_path = fullfile(OUT_ROOT, 'val_manifest.csv');
if ~exist(manifest_path, 'file')
    fid = fopen(manifest_path, 'w');
    fprintf(fid, 'case_name,k_mD,S_skin,schedule_id,schedule_label,n_timesteps,sim_wall_s,status\n');
    fclose(fid);
    completed = {};
else
    T = readtable(manifest_path, 'Delimiter', ',');
    if ~isempty(T) && ismember('status', T.Properties.VariableNames)
        ok_rows   = strcmp(T.status, 'ok');
        completed = T.case_name(ok_rows);
    else
        completed = {};
    end
end
fprintf('Already completed: %d / %d cases\n\n', numel(completed), n_total);

%% -- Main loop -------------------------------------------------------
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

            % -- Skip if already done ---------------------------------
            if any(strcmp(completed, case_name))
                n_skipped = n_skipped + 1;
                fprintf('[%3d/%3d] SKIP  %s\n', counter, n_total, case_name);
                continue
            end

            fprintf('[%3d/%3d] START %s\n', counter, n_total, case_name);

            % -- Run simulation (run_single_case.m unchanged) ---------
            t_case  = tic;
            status  = 'ok';
            n_t_out = NaN;

            try
                run_single_case(k_val, S_val, s_id, OUT_ROOT, MRST_ROOT);

                wd      = load(fullfile(OUT_ROOT, case_name, 'well_data.mat'), 'time_hr');
                n_t_out = numel(wd.time_hr);

            catch ME
                status   = 'error';
                n_errors = n_errors + 1;
                warning('BATCH:failed', '[%s] FAILED: %s', case_name, ME.message);
            end

            sim_s = toc(t_case);

            % -- Append to manifest -----------------------------------
            si  = make_rate_schedule(s_id);
            fid = fopen(manifest_path, 'a');
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

%% -- Summary ---------------------------------------------------------
fprintf('\n==========================================================\n');
fprintf(' Validation batch complete\n');
fprintf('  Total   : %d  |  New: %d  |  Skipped: %d  |  Errors: %d\n', ...
        counter, counter - n_skipped - n_errors, n_skipped, n_errors);
fprintf('  Wall time: %.1f min\n', toc(t_batch)/60);
fprintf('  Manifest : %s\n', manifest_path);
fprintf('==========================================================\n');
