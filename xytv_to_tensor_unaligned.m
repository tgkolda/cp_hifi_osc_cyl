function TU = xytv_to_tensor_unaligned(XYTV)
%XYTV_TO_TENSOR_UNALIGNED  Pack [x y t v] rows into a tensor_unaligned.
%
%   TU = XYTV_TO_TENSOR_UNALIGNED(XYTV)
%
%   XYTV is an N x 4 matrix whose columns are [x, y, t, v]. Returns a
%   tensor_unaligned TU of size [Nx Ny Nt] whose x-values (the real-
%   valued coordinates stored alongside the integer subs) are the sorted
%   unique values along each mode.
%
%   Assumes cp_hifi_code and Tensor Toolbox are already on the path.
%
%   Example:
%     XYTV = load_lbnl_yvel;
%     TU   = xytv_to_tensor_unaligned(XYTV);
%     size(TU)     % [Nx Ny Nt]
%     nnz(TU)      % number of observed entries

    x = XYTV(:,1);  y = XYTV(:,2);  t = XYTV(:,3);  v = XYTV(:,4);

    [xv, ~, ix] = unique(x);
    [yv, ~, iy] = unique(y);
    [tv, ~, it] = unique(t);

    subs  = [ix, iy, it];
    sz    = [numel(xv), numel(yv), numel(tv)];
    xvals = {xv(:).', yv(:).', tv(:).'};   % row vectors per class contract

    TU = tensor_unaligned(subs, v, sz, xvals);
end
