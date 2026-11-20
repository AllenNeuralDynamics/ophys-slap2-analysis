function [trialTable, dirs] = loadTrialTableLegacy(path)
%LOADTRIALTABLELEGACY Load a trialTable.mat or GIAnT-MATLAB trial_table.h5
%   as the legacy (camelCase) struct used by summarize_Voltage/verifyFiles.
%
%   [trialTable, dirs] = loadTrialTableLegacy(path)
%
%   dirs.dataDr : folder holding the raw .dat files and reference TIFFs
%   dirs.mocoDr : folder holding registered TIFFs and *_ALIGNMENTDATA files
%   dirs.saveDr : folder to write ExperimentSummary into
%
%   For a .mat trial table all three are the trial table folder (legacy
%   layout: everything lives in the neuron folder). For an .h5 trial table
%   written by GIAnT-MATLAB, raw data live in trialTable.datadr, motion
%   correction outputs in <savedr>/motion_correction, and results go to
%   trialTable.savedr. Filenames in fnRegDS/fnAdata/fnRaw stay relative;
%   callers resolve them against dirs.mocoDr / dirs.dataDr.
%
%   Field mapping (h5 -> legacy):
%     slap2_info.first_line / last_line  -> firstLine / lastLine
%     slap2_info.ref_stack.Path#         -> refStack{#}
%     motion_correction.fn_reg_ds        -> fnRegDS
%     motion_correction.fn_adata         -> fnAdata
%     source_extraction.fn_raw (or motion_correction.fn_raw, or filename) -> fnRaw
%     true_trial_ix                      -> trueTrialIx

[dr, ~, ext] = fileparts(char(path));
if isempty(dr)
    dr = pwd;
end

if strcmpi(ext, '.h5')
    s = readStructH5(path);
    trialTable = normalizeH5(s);
    dirs.dataDr = pickDir(s, 'datadr', dr);
    dirs.saveDr = pickDir(s, 'savedr', dr);
    dirs.mocoDr = fullfile(dirs.saveDr, 'motion_correction');
else
    S = load(path, 'trialTable');
    trialTable = S.trialTable;
    dirs.dataDr = dr;
    dirs.mocoDr = dr;
    dirs.saveDr = dr;
end
end

function d = pickDir(s, fieldName, fallback)
d = fallback;
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    v = s.(fieldName);
    if iscell(v), v = v{1}; end
    d = char(v);
end
end

function t = normalizeH5(s)
t = s;
filename = s.filename;
if isstring(filename), filename = cellstr(filename); end
if ~iscell(filename), filename = {filename}; end
t.filename = filename;
[nDMDs, nTrials] = size(filename);

if isfield(s, 'slap2_info')
    si = s.slap2_info;
    if isfield(si, 'first_line'), t.firstLine = double(si.first_line); end
    if isfield(si, 'last_line'),  t.lastLine  = double(si.last_line);  end
    if isfield(si, 'ref_stack') && isstruct(si.ref_stack)
        rs = si.ref_stack;
        t.refStack = cell(1, nDMDs);
        for dmdIx = 1:nDMDs
            keys = {sprintf('Path%d', dmdIx), sprintf('DMD%d', dmdIx)};
            for k = 1:numel(keys)
                if isfield(rs, keys{k})
                    t.refStack{dmdIx} = rs.(keys{k});
                    break
                end
            end
        end
    end
end

if isfield(s, 'true_trial_ix')
    t.trueTrialIx = s.true_trial_ix;
else
    t.trueTrialIx = 1:nTrials;
end

if isfield(s, 'motion_correction')
    mc = s.motion_correction;
    if isfield(mc, 'fn_reg_ds'), t.fnRegDS = toCell(mc.fn_reg_ds, [nDMDs nTrials]); end
    if isfield(mc, 'fn_adata'),  t.fnAdata = toCell(mc.fn_adata, [nDMDs nTrials]);  end
    if isfield(mc, 'registration_failed'), t.registrationFailed = logical(mc.registration_failed); end
end
if ~isfield(t, 'fnAdata')
    error('loadTrialTableLegacy:NotRegistered', ...
        'trial_table.h5 has no motion_correction/fn_adata. Run registration first.');
end
if ~isfield(t, 'fnRegDS')
    t.fnRegDS = t.fnAdata; % only used for an existence check in verifyFiles
end

% raw source files for voltage extraction: SLAP2 reads the original .dat files
if isfield(s, 'source_extraction') && isfield(s.source_extraction, 'fn_raw')
    t.fnRaw = toCell(s.source_extraction.fn_raw, [nDMDs nTrials]);
elseif isfield(s, 'motion_correction') && isfield(s.motion_correction, 'fn_raw')
    t.fnRaw = toCell(s.motion_correction.fn_raw, [nDMDs nTrials]);
else
    t.fnRaw = t.filename;
end
end

function c = toCell(v, sz)
if isstring(v)
    c = cellstr(v);
elseif ischar(v)
    c = {v};
else
    c = v;
end
if numel(c) == prod(sz)
    c = reshape(c, sz);
end
end
