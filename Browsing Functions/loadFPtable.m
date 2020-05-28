
function FPtable = loadFPtable(fn)
% load the lookup table from CCF to FP slices

[~, fnBase] = fileparts(fn);

fid = fopen(fn, 'r');

titles = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s', 1, 'delimiter', ',');
titles = cellfun(@(x)x{1}, titles, 'uni', false);
data = textscan(fid, '%s%s%d%d%d%d%d%d%s%s%s', 'delimiter', ',');

FPtable = table(data{:}, 'VariableNames', titles);


fclose(fid);