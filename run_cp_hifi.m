%% RUN CP HIFI on the oscillating-cylinder dataset
cd('C:\Users\TammyKolda\Documents\GIT Repos\CP-HIFI\lbnl_data');

%% --- Check paths ---
% Ensure tensor_toolbox is in path
if ~exist('tensor','file')
	warning('Tensor Toolbox not found in path. Attempting to add it now...');
    addpath(genpath('../../tensor_toolbox'));
end

% Add CP-ALS-HIFI to path
if ~exist('cp_als_hifi','file')
	warning('CP-ALS-HIFI not found in path. Adding it now...');
	addpath(genpath('../cp_hifi_code'));
end

%% --- Load the oscillating-cylinder slice ---
% The original data is 3D plus temporal.
% The x-values  go from -0.15 to 0.15, and we only consider x=0.
% The csv files were created with amrex_to_csv.py, using a step of 100 
% (so 25 timesteps total, from t=0 to t=2400) and exporting the y-z plane at x=0.
% We can go back and use a finer time step, down to t-delta=10 rather than t-delta=100.
% which exported the x=0 plane as a (y,z) slice over time. 
% So the tensor we're building is
%  mode 1 = y (roughly -5 to 5)
%  mode 2 = z (roughly 0 to 20)
%  mode 3 = t (time steps 0, 100, ..., 2400) - 25 timesteps total
%
% Pick a field from: `density`, `y_velocity`, `z_velocity`, `vorticity`.
% The `vfrac` field is also included which marks the embedded boundary 
% (cylinder) with vfrac=0; we use this to mask out the solid region since 
% it carries no meaningful flow data and would otherwise inject a spurious zero block into the tensor.

% Choose field
field = 'density';

% Read CVS files
% D columns: [y, z, t, <field>, vfrac].
if ~exist('D.mat','file')
    D = load_osc_cyl({field, 'vfrac'});
    % -v7.3 (HDF5+deflate): required once D tops 2 GB -- all 234 timesteps is
    % ~132M rows -- and gives a compressed cache MATLAB reloads directly.
    save('D.mat',"D","field",'-v7.3')
else
    load D.mat
end

%% --- Sanity-check the assembled data ---
% The full export is 234 de-duped slices of 565,888 points each, so D should
% have a uniform point count per timestep. A non-uniform count means a stale
% or duplicated CSV slipped in and doubled a timestep (the tensor builder sums
% coincident subscripts); a different total means a partial export. Soft
% warnings, not errors, so a deliberately reduced CSV set still runs.
if 0
expNt  = 234;
expPts = 565888;
[tvals_all, ~, gT] = unique(D(:,3));
nt      = numel(tvals_all);
ptsPerT = accumarray(gT, 1);            % rows per timestep
fprintf('data check: %d timesteps, %d total rows, %d-%d pts/timestep\n', ...
        nt, size(D,1), min(ptsPerT), max(ptsPerT));
if nt ~= expNt
    warning('run_cp_hifi:timesteps', ...
        'expected %d timesteps but found %d', expNt, nt);
end
if min(ptsPerT) ~= max(ptsPerT)
    warning('run_cp_hifi:nonuniform', ...
        ['non-uniform points-per-timestep (%d..%d): a duplicated or stale ' ...
         'CSV likely doubled a timestep. Re-export and delete D.mat.'], ...
        min(ptsPerT), max(ptsPerT));
elseif ptsPerT(1) ~= expPts
    warning('run_cp_hifi:ptcount', ...
        'uniform %d pts/timestep but expected %d', ptsPerT(1), expPts);
end
clear tvals_all gT ptsPerT nt expNt expPts
end
%% --- Build the tensor ---
% Spatial window: restrict to y in [-3,3] and z in [5,10].
ylim_keep = [-0.5 0.5];
zlim_keep = [6.5 8];

% Scale factor applied to the field values (mode 4) when building the tensor.
val_scale = 10;

% Switch: when true, drop solid-region cells (vfrac<=0.5); when false, keep them.
mask_vfrac = false;

% Mask cells outside the spatial window, and (optionally) the solid region.
keep = D(:, 1) >= ylim_keep(1) & D(:, 1) <= ylim_keep(2) ... % y window
     & D(:, 2) >= zlim_keep(1) & D(:, 2) <= zlim_keep(2);    % z window
if mask_vfrac
    keep = keep & D(:, 5) > 0.5;   % open-flow cells (vfrac>0.5)
    mask_desc = 'vfrac+window';
else
    mask_desc = 'window only';
end
YZTV = D(keep, [1 2 3 4]);        % [y z t <field>]
YZTV(:, 4) = val_scale * YZTV(:, 4);   % scale the field values
YZTV(:, 4) = YZTV(:, 4)  - median(YZTV(:, 4));
fprintf('mask: kept %d of %d cells (%.1f%% removed by %s)\n', ...
        sum(keep), numel(keep), 100*mean(~keep), mask_desc);
fprintf('value scale: field multiplied by %g\n', val_scale);

% Preserve the full (pre-subsample) observations for the comparison figure,
% so the "Orig" row shows all the data rather than the thinned tensor.
YZTV_full = YZTV;

% Random subsample: each YZTV row becomes one tensor nonzero, so cap the row
% count to keep nnz tractable. Set to Inf to disable. Seeded for reproducibility.
target_nnz = 1e6;
if size(YZTV,1) > target_nnz
    rng(0);
    keepIdx = randsample(size(YZTV,1), target_nnz, false);
    fprintf('subsample: %d -> %d rows (%.1f%% kept)\n', ...
            size(YZTV,1), target_nnz, 100*target_nnz/size(YZTV,1));
    YZTV = YZTV(keepIdx, :);
    clear keepIdx
else
    fprintf('subsample: %d rows <= target %d, keeping all\n', ...
            size(YZTV,1), target_nnz);
end

fprintf('Building tensor...\n');
tBuild = tic;
TU = xytv_to_tensor_unaligned(YZTV);
fprintf('  tensor built in %.1fs: size %s, nnz %d\n', ...
        toc(tBuild), mat2str(size(TU)), nnz(TU));

% clear keep D YZTV tBuild;
%% --- Fiber plots: shared setup fixed third mode to tidx ---

% fibers to draw per plot
nFibers = 200; 

% Set a target time in the range [0,2300] and then 
% let tidx be the index of the closest time point.
target_time = 1000;
[~, tidx] = min(abs(TU.xvals{3} - target_time)); 
clear target_time
fiber_time = TU.xvals{3}(tidx);

% Filter to just subscripts for time tidx
idx = (TU.subs(:,3) == tidx);
subs = TU.subs(idx, :);  
vals = TU.vals(idx);
ylabel_str = strrep(field, '_', '\_');

%% plot some mode-1 fibers (vary y; fix z and t)
figure(1); clf;

% Per-mode Gaussian overlays. Each mode's width is the RKHS kernel sigma used
% for that mode below (hifi_info{m}.kfunc = kernfunc_gaussian(gauss_width<m>)).
% Mode 1 (y):
gauss_width1  = 0.02;    % sigma
gauss_height1 = -0.6;      % peak amplitude
gauss_center1 = 0;      % center in y

all_zidx = unique(subs(:,2));
sampled_zidx = sort(randsample(all_zidx,nFibers,false));

hold on;
for j = 1:nFibers
    zidx = sampled_zidx(j);
    idx = (subs(:,2) == zidx);
    xx = subs(idx,1);
    yy = vals(idx,1);
    [xxs,ord] = sort(xx);
    xxv = TU.xvals{1}(xxs);
    plot(xxv, yy(ord), '.-');               % column vectors -> one clean line
end

% Overlay the Gaussian last so it sits on top of all the fibers.
yg = linspace(min(TU.xvals{1}), max(TU.xvals{1}), 400);
gg = gauss_height1 * exp(-((yg - gauss_center1).^2) ./ (2*gauss_width1.^2));
plot(yg, gg, 'k-', 'LineWidth', 2);
hold off;
title_str = sprintf('%d Mode-1 Fibers for time t=%g', nFibers, fiber_time);
title(title_str);
xlabel('y');
ylabel(ylabel_str);
grid on;

%% plot some mode-2 fibers (vary z; fix y and t)
figure(2); clf;
% Mode 2 (z):
gauss_width2  = 0.03;    % sigma
gauss_height2 = -0.7;      % peak amplitude
gauss_center2 = 7.25;   % center in z (mid-window)

all_yidx = unique(subs(:,1));
sampled_yidx = sort(randsample(all_yidx,nFibers,false));

hold on;
for j = 1:nFibers
    yidx = sampled_yidx(j);
    idx = (subs(:,1) == yidx);
    xx = subs(idx,2);
    yy = vals(idx,1);
    [xxs,ord] = sort(xx);
    xxv = TU.xvals{2}(xxs);
    plot(xxv, yy(ord), '.-');               % column vectors -> one clean line
end

% Overlay the Gaussian last so it sits on top of all the fibers.
zg = linspace(min(TU.xvals{2}), max(TU.xvals{2}), 400);
gg = gauss_height2 * exp(-((zg - gauss_center2).^2) ./ (2*gauss_width2.^2));
plot(zg, gg, 'k-', 'LineWidth', 2);
hold off;
title_str = sprintf('%d Mode-2 Fibers for time t=%g', nFibers, fiber_time);
title(title_str);
xlabel('z');
ylabel(ylabel_str);
grid on;


%% plot some mode-3 fibers (vary t; fix y and z)
% NOTE: this plot varies t, so it must use the FULL tensor data, not the
% time-filtered subs/vals built above (those hold only t = tidx).
figure(3); clf;
% Mode 3 (t):
gauss_width3  = 50;     % sigma (in t units)
gauss_height3 = -0.7;    % peak amplitude
gauss_center3 = 1000;    % center in t

zvals = TU.xvals{2};
yvals = TU.xvals{1};
tvals = TU.xvals{3};
zidx  = round(numel(zvals)/2);                % a fixed z-plane (mid-domain)
zz    = zvals(zidx);

atz = TU.subs(:,2) == zidx;                   % all points on this z-plane (all y, all t)
yi = TU.subs(atz,1);
ti = TU.subs(atz,3);
vv = TU.vals(atz);
yUsed = unique(yi);                           % y-fibers present on this z-plane
pick  = unique(yUsed(round(linspace(1, numel(yUsed), min(nFibers, numel(yUsed))))));

hold on;
for j = 1:numel(pick)
    m = (yi == pick(j));
    [ts, ord] = sort(ti(m));
        xx = tvals(ts);
        yy = vv(m);
    plot(xx(:), yy(ord), '.-');
end

% Overlay the Gaussian last so it sits on top of all the fibers.
tg = linspace(min(tvals), max(tvals), 400);
gg = gauss_height3 * exp(-((tg - gauss_center3).^2) ./ (2*gauss_width3.^2));
plot(tg, gg, 'k-', 'LineWidth', 2);
hold off;
title(sprintf('mode-3 fibers at z = %.4g  (%d of %d y-fibers)', ...
              zz, numel(pick), numel(yUsed)));
xlabel('t');
ylabel(ylabel_str);
grid on;
fprintf('mode-3: %d fibers at z=%.4g, value range [%.4g %.4g]\n', ...
        numel(pick), zz, min(vv), max(vv));

%% Downsample TU


%% Setup for CP-HIFI
% Kernel/regularization carried over from the vortex run as a starting
% point; these likely need retuning for this dataset's coordinate scales.
hifi_info{1}.inf         = true;
hifi_info{1}.kfunc       = kernfunc_gaussian(gauss_width1);  % matches the figure-1 overlay
hifi_info{1}.lambda      = 1e-3;
hifi_info{1}.rho         = 1e-6;
hifi_info{2}.inf         = true;
hifi_info{2}.kfunc       = kernfunc_gaussian(gauss_width2);  % matches the figure-2 overlay
hifi_info{2}.lambda      = 1e-3;
hifi_info{2}.rho         = 1e-6;
hifi_info{3}.inf         = true;
hifi_info{3}.kfunc       = kernfunc_gaussian(gauss_width3);  % matches the figure-3 overlay
hifi_info{3}.lambda      = 1e-3;
hifi_info{3}.rho         = 1e-6;

% Run cp-hifi
R = 100;
genargs = {'printitn',1,'tol',1e-4,'maxiters',10,...
	'solver',{'pcg','pcg','pcg'}};
fprintf('Running CP-HIFI: R=%d on size %s, nnz %d ...\n', ...
        R, mat2str(size(TU)), nnz(TU));
tHifi = tic;
M = cp_als_hifi(TU, R, hifi_info, genargs{:});
fprintf('CP-HIFI finished in %.1fs\n', toc(tHifi));

%%
figure(5); clf;
viz(extract(M,1:5),'Figure',5);

%% Display the resulting reconstruction
% Resample the HIFI modes (1 = y, 2 = z) onto a uniform fine grid so
% imagesc renders correctly. Mode 3 (t) is now also a HIFI mode, but we
% leave it on its native grid here since the display shows discrete timesteps.
fprintf('Reconstructing for display (resample + full)...\n');
tRec = tic;
Mu = resample_mode(M,  1, 100, [], false);
Mu = resample_mode(Mu, 2, 100, [] ,false);

Xhat = full(Mu);
fprintf('  reconstruction ready in %.1fs\n', toc(tRec));
yv = Xhat.xvals{1};
zv = Xhat.xvals{2};
tv = Xhat.xvals{3};

% Pick a half-dozen timesteps spread across the available time range.
nShow = 6;
tsel  = unique(round(linspace(1, numel(tv), nShow)));

% Shared axis limits; color limits driven by the *data* range across all
% shown timesteps so the originals always render correctly. If the
% reconstruction overshoots, it will visibly saturate -- a feature, not a bug.
mSel = ismember(YZTV_full(:,3), tv(tsel));
cl = [min(YZTV_full(mSel,4)), max(YZTV_full(mSel,4))];
yl = [min(yv) max(yv)];
zl = [min(zv) max(zv)];

% Layout: 3 rows x nShow cols. Row 1 = full data, row 2 = sampled data (what
% the tensor was fit on), row 3 = reconstruction; each column is one timestep.
nRows = 3;
figure(4); clf;
set(gcf, 'Units', 'normalized', 'Position', [0.04 0.12 0.92 0.78]);  % wide rectangle, fits screen
for k = 1:numel(tsel)
    ti = tsel(k);
    tt = tv(ti);
    rec = squeeze(Xhat.data(:,:,ti));   % [Ny x Nz]

    % Full data: all observed points at this timestep.
    mf = YZTV_full(:,3) == tt;
    % Sampled data: the thinned points actually fed to the tensor.
    ms = YZTV(:,3) == tt;

    % Row 1: full data.
    subplot(nRows, numel(tsel), k);
    scatter(YZTV_full(mf,1), YZTV_full(mf,2), 18, YZTV_full(mf,4), 'filled', 's');
    axis equal;
    xlim(yl);
    ylim(zl);
    clim(cl);
    title(sprintf('Full t = %.4g', tt));
    if k == 1, ylabel('z'); end
    xlabel('y');

    % Row 2: sampled data.
    subplot(nRows, numel(tsel), numel(tsel) + k);
    scatter(YZTV(ms,1), YZTV(ms,2), 18, YZTV(ms,4), 'filled', 's');
    axis equal;
    xlim(yl);
    ylim(zl);
    clim(cl);
    title(sprintf('Sampled t = %.4g', tt));
    if k == 1, ylabel('z'); end
    xlabel('y');

    % Row 3: reconstruction.
    subplot(nRows, numel(tsel), 2*numel(tsel) + k);
    imagesc(yv, zv, rec.');
    axis xy equal;
    xlim(yl);
    ylim(zl);
    clim(cl);
    title(sprintf('Recon t = %.4g', tt));
    if k == 1, ylabel('z'); end
    xlabel('y');
end

% One shared colorbar for the whole figure.
colorbar('Position', [0.93 0.11 0.015 0.815]);

%% Side-by-side movie: full | sampled | reconstruction, animated over time
% Three panels in lockstep across every timestep. We draw each panel once,
% keep its graphics handle, then update only the data each frame (CData /
% XData / YData) so playback stays smooth instead of re-plotting every frame.
fps = 10;                       % playback rate
nT  = numel(tv);

figure(6); clf;
set(gcf, 'Units', 'normalized', 'Position', [0.04 0.30 0.92 0.40]);

% Initialize each panel from the first timestep, capture handles.
ti = 1; tt = tv(ti);
mf = YZTV_full(:,3) == tt;
ms = YZTV(:,3)      == tt;

ax1 = subplot(1,3,1);
hFull = scatter(YZTV_full(mf,1), YZTV_full(mf,2), 18, YZTV_full(mf,4), 'filled', 's');
axis equal; xlim(yl); ylim(zl); clim(cl); xlabel('y'); ylabel('z'); title('Full');

ax2 = subplot(1,3,2);
hSamp = scatter(YZTV(ms,1), YZTV(ms,2), 18, YZTV(ms,4), 'filled', 's');
axis equal; xlim(yl); ylim(zl); clim(cl); xlabel('y'); title('Sampled');

ax3 = subplot(1,3,3);
hRec = imagesc(yv, zv, squeeze(Xhat.data(:,:,ti)).');
axis xy equal; xlim(yl); ylim(zl); clim(cl); xlabel('y'); title('Reconstruction');

colorbar('Position', [0.93 0.11 0.015 0.815]);
hSup = sgtitle(sprintf('t = %.4g   (frame %d/%d)', tt, ti, nT));

% Animate. Loops once through all timesteps; wrap in `while ishandle(...)` if
% you want it to repeat. Updates data in place for smooth playback.
for ti = 1:nT
    tt = tv(ti);
    mf = YZTV_full(:,3) == tt;
    ms = YZTV(:,3)      == tt;

    set(hFull, 'XData', YZTV_full(mf,1), 'YData', YZTV_full(mf,2), 'CData', YZTV_full(mf,4));
    set(hSamp, 'XData', YZTV(ms,1),      'YData', YZTV(ms,2),      'CData', YZTV(ms,4));
    set(hRec,  'CData', squeeze(Xhat.data(:,:,ti)).');
    set(hSup,  'String', sprintf('t = %.4g   (frame %d/%d)', tt, ti, nT));

    drawnow;
    pause(1/fps);
end

%% Storage comparison

fprintf('Original data: %d nonzeros\n', nnz(TU));
fprintf('Total storage in MB for scarce tensor: %.2f\n', ndims(TU)*nnz(TU)*16/1e6);

fprintf('CP-HIFI model: %d components\n', R);
fprintf('Total storage in MB for CP-HIFI: %.2f\n', R*sum(size(M))*16/1e6);
