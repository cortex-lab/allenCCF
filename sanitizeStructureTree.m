

function sanitizeStructureTree(inFilename, outFilename)
% function sanitizeStructureTree(inFilename, outFilename)
%
% make the structure tree that you get from:
% 
% into something that you can load into a table with loadStructureTree.m
%
% This means taking out all the commas within fields. 

fid = fopen(inFilename, 'r');
ofid = fopen(outFilename, 'w');

while ~feof(fid)
    q = textscan(fid, '%s', 1, 'delimiter', '\n');
    
    thisLine = q{1}{1};
    
    doubleQuote = strfind(thisLine, '"');
    while ~isempty(doubleQuote)
        commas = strfind(thisLine, ',');
        thisLine(commas(commas>doubleQuote(1) & commas<doubleQuote(2))) = '';
        doubleQuote = doubleQuote(3:end);
    end
    
    fprintf(ofid, '%s\n', thisLine); 
end

fclose(fid);        
fclose(ofid);