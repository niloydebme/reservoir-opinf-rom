function write_json(filepath, s)
%WRITE_JSON  Write a flat struct to a JSON file.
%  Handles: scalars, strings, and 1-D numeric arrays.
%
%  write_json(filepath, s)

fid = fopen(filepath, 'w');
if fid == -1
    error('write_json: cannot open file for writing: %s', filepath);
end

fn = fieldnames(s);
fprintf(fid, '{\n');
for i = 1 : numel(fn)
    v = s.(fn{i});
    fprintf(fid, '  "%s": ', fn{i});
    if ischar(v)
        fprintf(fid, '"%s"', v);
    elseif isnumeric(v) && isscalar(v)
        if isnan(v)
            fprintf(fid, 'null');
        else
            fprintf(fid, '%.10g', v);
        end
    elseif isnumeric(v) && numel(v) > 1
        % 1-D numeric array  →  JSON array
        v = v(:)';   % ensure row
        fprintf(fid, '[');
        fprintf(fid, '%.10g, ', v(1 : end-1));
        fprintf(fid, '%.10g]', v(end));
    elseif isempty(v)
        fprintf(fid, '[]');
    else
        fprintf(fid, 'null');
    end
    if i < numel(fn), fprintf(fid, ','); end
    fprintf(fid, '\n');
end
fprintf(fid, '}\n');
fclose(fid);

end