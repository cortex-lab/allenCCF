
function [treePairsInd, treePairs] = makeSTtree(st)
% make structure tree into a list of parent/child relationships. This makes
% it very quick to find all indices that are a child of another, etc. 
%
% Use "loadStructureTree" to get the input argument
%
% output has three columns: parent, child, number of levels between
% e.g. [4 6 2] means that 4 is a grandparent of 6 (2 levels above)
%
% treePairs uses the "id"s of the structures, treePairsInd uses the
% indices, i.e. the row numbers.
%
% This makes several operations trivial, e.g.: 
% - Find all direct children of a certain structure:
%  >> ch = tp(tp(:,1)==targetID & tp(:,3)==1,2);
% - Find siblings of a certain structure:
%  >> parentID = tp(tp(:,2)==targetID & tp(:,3)==1,1);
%  >> sib = % as above, for children
% - Find all structures anywhere below a target:
%  >> below = tp(tp(:,1)==targetID,2);

% st = loadStructureTree('J:/allen/structure_tree_safe.csv');

treePairs = [];
treePairsInd = [];
ids = st.id;
for s = 1:size(st,1)
    spath = st.structure_id_path{s};
    nums = regexp(spath, '/*/', 'split');
    nums = cellfun(@str2double, nums(2:end-1));
    
    inds = arrayfun(@(x)find(ids==x), nums);
    
    for x = 1:numel(nums)-1
        
        % treePairs has:
        % - parentID, childID, numLevelsBetween (1=direct parent/child
        % relationship, 2=grandchild, etc)
        treePairs(end+1,:) = [nums(x) nums(end) numel(nums)-x];
        
        treePairsInd(end+1,:) = [inds(x) inds(end) numel(nums)-x];
        
    end
end