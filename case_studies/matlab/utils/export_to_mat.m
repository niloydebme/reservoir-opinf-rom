function export_to_mat(outdir, states, wellSols, schedule, rock, G, params)
%EXPORT_TO_MAT  Package MRST simulation results .mat files
%
%  Arguments
%  ---------
%  outdir    : string, output directory (created if missing)
%  states    : cell array from simulateScheduleAD
%  wellSols  : cell array of well solutions
%  schedule  : MRST schedule struct
%  rock      : rock struct (perm, poro)
%  G         : grid struct
%  params    : struct of scalar parameters (phi, ct, mu, q, rw, etc.)

    if ~exist(outdir, 'dir'), mkdir(outdir); end
    n_t  = numel(states);         % number of timesteps
    n_c  = G.cells.num;           % number of grid cells

    %--- 1. Snapshot matrix X : pressure [Pa] at each cell × timestep ---
    X = zeros(n_c, n_t);
    for k = 1:n_t
        X(:, k) = states{k}.pressure;   % [Pa]
    end
    save(fullfile(outdir, 'snapshots.mat'), 'X', '-v7.3');
    fprintf('Saved snapshots.mat  [%d cells × %d timesteps]\n', n_c, n_t);

    %--- 2. BHP time series at each well ---
    n_w = numel(wellSols{1});
    BHP = zeros(n_w, n_t);
    Q   = zeros(n_w, n_t);
    for k = 1:n_t
        for w = 1:n_w
            BHP(w, k) = wellSols{k}(w).bhp;   % [Pa]
            Q(w, k)   = wellSols{k}(w).qOs;   % [m³/s] surface oil rate
        end
    end
    time_days = cumsum(schedule.step.val) / 86400;  % convert s → days

    save(fullfile(outdir, 'bhp_rate.mat'), 'BHP', 'Q', 'time_days', '-v7.3');
    fprintf('Saved bhp_rate.mat   [%d wells × %d timesteps]\n', n_w, n_t);

    %--- 3. Rock properties ---
    perm = rock.perm;    % [m²]  shape: [n_c × 1] or [n_c × 3]
    poro = rock.poro;    % [-]   shape: [n_c × 1]
    save(fullfile(outdir, 'rock.mat'), 'perm', 'poro', '-v7.3');

    %--- 4. Grid geometry (cell centroids only — compact) ---
    centroids = G.cells.centroids;  % [n_c × 3]  in metres
    save(fullfile(outdir, 'grid.mat'), 'centroids', '-v7.3');

    %--- 5. Transmissibility matrix A ---
    %  A = T * D  where T = TPFA transmissibility, D = divergence operator
    %  Useful for verifying your SPD constraint in Python
    T = computeTrans(G, rock);
    save(fullfile(outdir, 'transmissibility.mat'), 'T', '-v7.3');

    %--- 6. Scalar parameters as JSON ---
    fid = fopen(fullfile(outdir, 'params.json'), 'w');
    fprintf(fid, '{\n');
    fields = fieldnames(params);
    for i = 1:numel(fields)
        val = params.(fields{i});
        if ischar(val)
            fprintf(fid, '  "%s": "%s"', fields{i}, val);
        else
            fprintf(fid, '  "%s": %g', fields{i}, val);
        end
        if i < numel(fields), fprintf(fid, ','); end
        fprintf(fid, '\n');
    end
    fprintf(fid, '}\n');
    fclose(fid);
    fprintf('Saved params.json\n');
    fprintf('Done. Output in: %s\n', outdir);
end