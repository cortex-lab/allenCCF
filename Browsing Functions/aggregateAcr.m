
function acrOut = aggregateAcr(acrIn, varargin)
% function acrOut = aggregateAcr(acrIn, varargin)
% 
% Translate voxel acronyms into a useful list to split analysis by. 
% 
% e.g.:
% >> acrIn = {'VISp5', 'VISp6a', 'MB', 'VISpm1', 'MB', 'LP', 'VISp5', 'LGd', 'VISp'};
% >> acrOut = aggregateAcr(acrIn); acrOut
% acrOut =
%   1×9 cell array
%     {'VISp'}    {'VISp'}    {[NaN]}    {'VISpm'}    {[NaN]}    {'LP'}    {'VISp'}    {'LGd'}    {'VISp'}

% Optional input argument aggType can be 1 for a short list of high-level
% regions, or 2 [default] for a long list of sensible regions (manually
% selected by N.A.S.). E.g. 
% >> acrOut = aggregateAcr(acrIn, 1); acrOut
% acrOut =
%   1×9 cell array  
%     {'Isocortex'} {'Isocortex'} {'MB'} {'Isocortex'} {'MB'} {'TH'} {'Isocortex'} {'TH'} {'Isocortex'}

% This function will return NaN for the acronyms that don't have a parent
% in the list or that aren't valid acronyms. To find these in the output,
% use: 
% >> cellfun(@(x)any(isnan(x)), acrOut)

% If you are calling this function many times, you can speed up execution
% by providing the tree-pair representation of the structure tree, which is
% a bit slow to calculate. Provide this as the third input argument. You
% can make it like this:
% >> [~, tp] = makeSTtree(st);

st = loadStructureTree();

if nargin>2
    tp = varargin{2};
else
    % see doc of this function to know what tp is
    [~, tp] = makeSTtree(st);
end

if nargin>1
    aggType = varargin{1}; 
    if isempty(aggType); aggType = 2; end
else
    aggType = 2; 
end

aggList = getAggList(aggType); 

aggTP = tp(ismember(tp(:,1), aggList),:); 
% aggTP now contains only children of desired list
% since these lists are chosen such that every entry of structure tree only
% has one (at most) parent on this list, the second column now contains
% just one entry for each of our desired areas.

% we also want exact matches so let's just add those in here. 
aggTP = [aggTP; horzcat(aggList', aggList', zeros(size(aggList))')];


if ~iscell(acrIn) % can be either a single acronym string, or a cell array
    acrIn = {acrIn};
end

allAcr = st.acronym; 
allId = st.id;

notAnAcr = ~ismember(acrIn, allAcr);
if any(notAnAcr)
    firstBadAcr = acrIn{find(notAnAcr,1)};
    warning(...
        sprintf('One or more of your acronyms didn''t parse. One that didn''t is ''%s''.', ...
        firstBadAcr));
end

acrInClean = acrIn(~notAnAcr); % only process the good ones, the rest we set to NaN later

idIn = cellfun(@(x)allId(strcmp(allAcr,x)), acrInClean);

tpLoc = arrayfun(@(x)find(aggTP(:,2)==x), idIn, 'uni', false); 

noMatch = cellfun(@isempty, tpLoc);
multipleMatch = cellfun(@(x)numel(x)>1, tpLoc);

if any(multipleMatch)
    % -- todo -- report one that is bad
    warning('Some acronyms appear to be children of two or more parents in the list. This probably means your aggList is badly constructed. Will return the first match.');
    tpLoc = cellfun(@(x)x(1), tpLoc, 'uni', false); 
end

% not throwing a warning here because this is probably common/normal, e.g.
% voxels that are "MB" or "root" or other things. 
% These are returned as NaN.
tpLocClean = tpLoc(~noMatch);

tpLocClean = cell2mat(tpLocClean);

aggParent = aggTP(tpLocClean,1);

aggOut = arrayfun(@(x)allAcr{allId==x}, aggParent, 'uni', false);

% now put in NaN's for the ones that didn't work
acrOut1 = cell(size(tpLoc));
acrOut1(noMatch) = {NaN};
acrOut1(~noMatch) = aggOut;

acrOut = cell(size(acrIn)); 
acrOut(notAnAcr) = {NaN};
acrOut(~notAnAcr) = acrOut1;
    
function aggList = getAggList(aggType)

switch aggType
    case 1 % 10 top-level divisions: {'Isocortex', 'HPF', 'OLF', 'CTXsp', 'CNU', 'TH', 'HY', 'MB', 'HB', 'CB'}
        aggList = [315 1089 698 703 623 549 1097 313 1065 512];
        
    case 2 % 316 manually selected regions
        aggList = [184 985 993 353 329 337 345 369 361 182305689 ...
                    378 1057 677 1011 1002 1027 1018 402 394 ...
                    409 385 425 533 312782574 312782628 39 48 972 44 723 ...
                    731 738 746 104 111 119 894 879 886 312782546 417 ...
                    541 922 895 507 151 159 597 605 814 961 619 639 647 ...
                    788 566 382 423 463 726 982 19 918 926 843 1037 1084 ...
                    502 484682470 589508447 484682508 583 952 966 131 295 ...
                    319 780 672 56 998 754 250 258 266 310 333 23 292 536 ...
                    1105 403 1022 1031 342 298 564 596 581 351 629 685 718 ...
                    725 733 741 563807435 406 609 1044 475 170 218 1020 1029 ...
                    325 560581551 255 127 64 1120 1113 155 59 362 366 1077 ...
                    149 15 181 560581559 189 599 907 575 930 560581563 262 ...
                    27 563807439 178 321 483 186 390 38 30 118 ...
                    223 72 263 272 830 452 523 1109 126 133 347 286 338 ...
                    576073699 689 88 210 491 525 557 515 980 1004 63 693 ...
                    946 194 226 364 576073704 173 470 614 797 302 4 580 271 ...
                    874 381 749 607344830 246 128 294 795 215 531 ...
                    628 634 706 1061 549009203 616 214 35 549009211 975 115 ...
                    606826663 757 231 66 75 58 374 1052 12 100 197 591 872 ...
                    612 7 867 398 280 880 599626927 898 931 1093 318 534 574 ...
                    621 549009215 549009219 549009223 549009227 679 147 162 ...
                    604 146 238 350 358 207 96 101 711 1039 903 642 651 429 ...
                    437 445 589508451 653 661 135 839 1048 372 83 136 106 ...
                    203 235 307 395 852 859 938 177 169 995 1069 209 202 225 ...
                    217 765 773 781 206 230 222 912 976 984 1091 936 944 951 ...
                    957 968 1007 1056 1064 1025 1033 1041 1049 989 91 846 ...
                    589508455];

end
