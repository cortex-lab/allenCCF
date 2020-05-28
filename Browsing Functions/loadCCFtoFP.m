
function CCFtoFPtable = loadCCFtoFP(fn)
% load the lookup table from CCF to FP slices

[~, fnBase] = fileparts(fn);

fid = fopen(fn, 'r');


titles = textscan(fid, '%s%s%s%s%s%s', 1, 'delimiter', ',');
titles = cellfun(@(x)x{1}, titles, 'uni', false);
data = textscan(fid, '%d%s%s%d%f%f', 'delimiter', ',');


CCFtoFPtable = table(data{:}, 'VariableNames', titles);

fclose(fid);