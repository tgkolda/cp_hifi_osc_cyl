function [XYTV, T, info] = load_osc_cyl(field, csvdir, verbose)
%LOAD_OSC_CYL  Load ParaView-exported oscillating-cylinder slice for CP-HIFI.
%
%   [XYTV, T, info] = LOAD_OSC_CYL()
%   [XYTV, T, info] = LOAD_OSC_CYL(field)
%   [XYTV, T, info] = LOAD_OSC_CYL(field, csvdir)
%   [XYTV, T, info] = LOAD_OSC_CYL(field, csvdir, verbose)
%
%   Reads the per-timestep CSVs written by amrex_to_csv.py (one file per
%   plotfile, all fields as columns) and returns coordinate-format data,
%   ready for XYTV_TO_TENSOR_UNALIGNED.
%
%   Each CSV has a header row of column names, e.g.
%       Points:0, Points:1, Points:2, density, vfrac, vorticity,
%       y_velocity, z_velocity
%   For the default x=0 slice, Points:0 (x) is constant; the two varying
%   coordinates are Points:1 (y) and Points:2 (z). This loader detects the
%   constant coordinate automatically, so it also works for a y= or z= slice.
%
%   Speed: the files are large (~879k rows each). We read only the columns
%   we need (the three Points + the requested field(s)) with a hard-coded
%   textscan format, which is ~3.5x faster than readmatrix, and print
%   per-file progress because a full load is ~tens of seconds.
%
%   Inputs:
%     field   : field column(s) to return as value(s). Either a single name
%               (char) or a cell array of names for a single-pass multi-field
%               load (e.g. {'z_velocity','vfrac'} reads the file set once).
%               Default 'z_velocity'. Any exported field works: 'y_velocity',
%               'density', 'vorticity', 'vfrac'.
%     csvdir  : directory holding osc_cyl__*.csv. Defaults to
%               <this file's folder>/csv.
%     verbose : print per-file progress (default true).
%
%   Returns:
%     XYTV : N x (3+F) matrix. Columns are [c1, c2, t, v1, ..., vF], where
%            c1,c2 are the two non-constant slice coordinates (for the x=0
%            slice: y, z), t is the integer plt index, and v1..vF are the
%            requested field(s) in the order given. For a single field this
%            is the usual N x 4 [c1, c2, t, v].
%     T    : 1 x Nt sorted vector of plt indices (e.g. 0, 100, 200, ...).
%     info : struct describing the load: .fields (cellstr), .valueCols (the
%            column indices of v1..vF in XYTV), .coordNames (the two kept
%            Points:* columns), .constName (the dropped constant axis),
%            .nDupesDropped, .files (cellstr, in time order).
%
%   Example (single field):
%     XYTV = load_osc_cyl('z_velocity');
%     TU   = xytv_to_tensor_unaligned(XYTV);
%
%   Example (one-pass field + mask, as in run_cp_hifi.m):
%     [D, T, info] = load_osc_cyl({'z_velocity','vfrac'});
%     keep = D(:, info.valueCols(2)) > 0.5;   % open-flow cells (vfrac>0.5)
%     XYTV = D(keep, [1 2 3 info.valueCols(1)]);
%
%   Notes:
%     - Cells inside the cylinder have vfrac == 0 (embedded boundary); mask
%       them as above.
%     - For CP-HIFI only the relative spacing of t matters, and the plt
%       stride is uniform, so the integer plt index is a fine time axis.
%
%   See README.md in this directory for the full export recipe.

    if nargin < 1 || isempty(field)
        field = 'z_velocity';
    end
    if nargin < 2 || isempty(csvdir)
        csvdir = fullfile(fileparts(mfilename('fullpath')), 'csv');
    end
    if nargin < 3 || isempty(verbose)
        verbose = true;
    end
    fields = cellstr(field);          % accept a char or a cellstr of names
    nf = numel(fields);

    % Accept plain CSVs and gzip-compressed CSVs (osc_cyl__<n>.csv.gz, written
    % by amrex_to_csv.py --gzip). A '*.csv' glob does not match '*.csv.gz', so
    % the two lists are disjoint and concatenating cannot double-count.
    files = [dir(fullfile(csvdir, 'osc_cyl__*.csv')); ...
             dir(fullfile(csvdir, 'osc_cyl__*.csv.gz'))];
    if isempty(files)
        error('load_osc_cyl:noFiles', ...
              'No osc_cyl__*.csv[.gz] files found in %s', csvdir);
    end

    % Sort files by the numeric plt index in the filename (a plain sort
    % would put osc_cyl__1000.csv before osc_cyl__200.csv). The optional
    % .gz suffix is ignored when reading the index.
    idx = zeros(numel(files), 1);
    for k = 1:numel(files)
        tok = regexp(files(k).name, 'osc_cyl__(\d+)\.csv(\.gz)?$', 'tokens', 'once');
        idx(k) = str2double(tok{1});
    end
    [T, order] = sort(idx(:).');
    files = files(order);

    % Map column names -> indices from the header of the first file. The
    % header has quoted names like "Points:0"; strip quotes and whitespace.
    [hpath, hclean] = resolve_csv(fullfile(csvdir, files(1).name)); %#ok<NASGU>
    hdr = read_header(hpath);
    clear hclean;                          % delete the temp now if it was a .gz
    ncols = numel(hdr);
    cP = [find_col(hdr, 'Points:0'), ...
          find_col(hdr, 'Points:1'), ...
          find_col(hdr, 'Points:2')];
    if any(isnan(cP))
        error('load_osc_cyl:noPoints', ...
              'Missing a Points:* column. Header: %s', strjoin(hdr, ', '));
    end
    cV = zeros(1, nf);
    for j = 1:nf
        cV(j) = find_col(hdr, fields{j});
        if isnan(cV(j))
            error('load_osc_cyl:noField', ...
                  'Field "%s" not found. Available columns: %s', ...
                  fields{j}, strjoin(hdr, ', '));
        end
    end

    % We only ever need the three Points columns plus the field column(s).
    % Read just those (textscan skips the rest with %*f). keepIdx is sorted,
    % so build a map from original column index -> position in the read matrix.
    keepIdx = unique([cP, cV]);            % sorted, distinct
    colpos  = zeros(1, ncols);
    colpos(keepIdx) = 1:numel(keepIdx);
    pPoints = colpos(cP);                  % positions of Points:0/1/2 in M
    pV      = colpos(cV);                  % positions of the field(s) in M

    if verbose
        fprintf('load_osc_cyl: %d files, field(s) {%s}, from %s\n', ...
                numel(files), strjoin(fields, ', '), csvdir);
    end

    blocks = cell(numel(files), 1);
    keepLocal = [];                        % which two Points cols vary (set on file 1)
    t0 = tic;
    for k = 1:numel(files)
        tf = tic;
        [rpath, rclean] = resolve_csv(fullfile(csvdir, files(k).name)); %#ok<NASGU>
        M = read_numeric_csv(rpath, ncols, keepIdx);
        clear rclean;                      % delete the temp now if it was a .gz
        P = M(:, pPoints);                 % the three coordinate columns
        V = M(:, pV);                      % the field value column(s), N x nf

        % On the first file, decide which Points column is the (constant)
        % slice-normal axis: the narrowest-ranging one. The other two are the
        % in-plane coordinates we keep.
        if isempty(keepLocal)
            ranges = max(P, [], 1) - min(P, [], 1);
            [~, constLocal] = min(ranges);
            keepLocal = setdiff(1:3, constLocal);

            % De-dup is now done once, at export time, by amrex_to_csv.py: the
            % default x=0 slice lies on a cell face, so ParaView emits two
            % coincident cell-centers per (c1,c2); the exporter drops the copy.
            % We therefore no longer de-dup every file here. Guard on the FIRST
            % file that the CSVs really are clean -- coincident (c1,c2) left in
            % would be DOUBLED by the tensor builder (which sums duplicate
            % subscripts). A failure means stale, pre-de-dup CSVs: re-export.
            c12 = [P(:, keepLocal(1)), P(:, keepLocal(2))];
            if size(unique(c12, 'rows'), 1) ~= size(c12, 1)
                error('load_osc_cyl:staleDuplicates', ...
                    ['%s holds coincident (c1,c2) points. These CSVs predate ' ...
                     'the de-dup in amrex_to_csv.py; the tensor builder would ' ...
                     'double their values. Re-export the slices (and delete ' ...
                     'any cached D.mat) before loading.'], files(k).name);
            end
        end

        blocks{k} = [P(:, keepLocal(1)), P(:, keepLocal(2)), ...
                     repmat(T(k), size(M, 1), 1), V];

        if verbose
            fprintf('  [%2d/%2d] %-20s %8d rows, %.2fs\n', ...
                    k, numel(files), files(k).name, size(blocks{k}, 1), toc(tf));
        end
    end
    XYTV = vertcat(blocks{:});

    nDup = 0;                              % retained for back-compat; CSVs are pre-de-duped
    if verbose
        fprintf('load_osc_cyl: done -- %d rows, %.1fs total\n', ...
                size(XYTV, 1), toc(t0));
    end

    pointNames = {'Points:0', 'Points:1', 'Points:2'};
    info = struct( ...
        'fields',        {fields}, ...
        'valueCols',     3 + (1:nf), ...
        'coordNames',    {pointNames(keepLocal)}, ...
        'constName',     pointNames{setdiff(1:3, keepLocal)}, ...
        'nDupesDropped', nDup, ...
        'files',         {{files.name}});
end


function M = read_numeric_csv(csvfile, ncols, keepIdx)
%READ_NUMERIC_CSV  Fast read of selected columns from a 1-header numeric CSV.
%   Returns an N x numel(keepIdx) matrix whose columns are the requested
%   columns in ascending index order. Uses a hard-coded textscan format
%   (%f for kept columns, %*f to skip the rest), which is markedly faster
%   than readmatrix for these large, uniformly-numeric files.
    spec = repmat({'%*f'}, 1, ncols);
    spec(keepIdx) = {'%f'};
    fmt = strjoin(spec, '');
    fid = fopen(csvfile, 'r');
    if fid < 0
        error('load_osc_cyl:open', 'Cannot open %s', csvfile);
    end
    cleanup = onCleanup(@() fclose(fid));
    C = textscan(fid, fmt, 'Delimiter', ',', 'HeaderLines', 1, ...
                 'CollectOutput', true);
    M = C{1};
end


function [p, c] = resolve_csv(fullpath)
%RESOLVE_CSV  Plain-CSV path for a .csv or .csv.gz file.
%   For a gzip-compressed CSV, decompress it to a fresh temp directory and
%   return that path plus an onCleanup handle that removes the temp when it
%   goes out of scope (so only one timestep is ever decompressed on disk at a
%   time). For a plain CSV, return the path unchanged with an empty cleanup.
    if endsWith(fullpath, '.gz')
        td = tempname;
        names = gunzip(fullpath, td);      % writes <td>/<basename without .gz>
        p = names{1};
        c = onCleanup(@() rmtemp(td));
    else
        p = fullpath;
        c = [];
    end
end


function rmtemp(td)
%RMTEMP  Remove a temp directory (and its decompressed CSV) if it exists.
    if isfolder(td)
        rmdir(td, 's');
    end
end


function names = read_header(csvfile)
%READ_HEADER  Return the CSV column names as a cellstr (quotes stripped).
    fid = fopen(csvfile, 'r');
    if fid < 0
        error('load_osc_cyl:open', 'Cannot open %s', csvfile);
    end
    line = fgetl(fid);
    fclose(fid);
    parts = strsplit(line, ',');
    names = strtrim(erase(parts, '"'));
end


function c = find_col(names, target)
%FIND_COL  Index of `target` in the cellstr `names`, or NaN if absent.
    c = find(strcmp(names, target), 1);
    if isempty(c)
        c = NaN;
    end
end
