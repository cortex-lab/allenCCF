

function isA = isAreaOrContains(testAcr, targetAcr, st)
% function isA = isAreaOrContains(testAcr, targetAcr, st)
%
% Determines whether an area with acronym given in testAcr is the same as
% targetAcr or is below targetAcr in the hierarchy. E.g.:
% isAreaOrContains('MG', 'TH') = true (medial geniculate is in the
% thalamus)
% isAreaOrContains('MG', 'MG') = true
% isAreaOrContains('TH', 'MG') = false
%
% testAcr can be a cell array, in which case isA is a vector of the same
% size

if ischar(testAcr)
    testAcr = {testAcr};
end
assert(iscellstr(testAcr),...
    'testAcr must be a string or cell array of strings, the acronyms of query areas');

assert(ischar(targetAcr), 'targetAcr must be a single string acronym')

targetInd = strcmp(st.acronym, targetAcr);
assert(~isempty(targetInd), 'targetAcr %s not found in the structure tree', targetAcr);
targetID = st(targetInd,:).id;
targetStr = sprintf('/%d/', targetID);

assert(all(cellfun(@(x)any(strcmp(st.acronym,x)), testAcr)), 'one or more test acronyms not found');

isA = cellfun(@(x)contains(st(strcmp(st.acronym, x),:).structure_id_path, targetStr),testAcr);
