function summary = extractDendrites_new(varargin)
%EXTRACTDENDRITES_NEW Deprecated compatibility wrapper for extractDendrites.
%
%   The epoch-aware HDF5-backed voltage extraction pipeline now lives in
%   extractDendrites.m. This wrapper preserves older notebooks/scripts that still
%   call extractDendrites_new while avoiding two independent implementations.
%
%   Prefer:
%       summary = extractDendrites(sessionDir, params);

warning('extractDendrites_new:Deprecated', ...
    ['extractDendrites_new is deprecated. Use extractDendrites instead; ' ...
     'this wrapper will call extractDendrites with the same inputs.']);
summary = extractDendrites(varargin{:});
end
