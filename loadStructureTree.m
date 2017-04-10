
function structureTreeTable = loadStructureTree(fn)
% note: use an edited version that doesn't have commas in any fields... 

% fn = 'structure_tree_safe.csv';

fid = fopen(fn, 'r');

titles = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s', 1, 'delimiter', ',');
titles = cellfun(@(x)x{1}, titles, 'uni', false);
titles{1} = 'index'; % this is blank in the file

data = textscan(fid, '%d%s%d%s%d%s%d%d%d%d%d%s%s%d%d%s%d%s%s%d%d', 'delimiter', ',');

structureTreeTable = table(data{:}, 'VariableNames', titles);

fclose(fid);