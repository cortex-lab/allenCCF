
function plotNeuronOnSliceFromCoord(ccfCoord, acr, ax, av, allenst, regColor)
% function plotNeuronOnSliceFromCoord(ccfCoord, acr, ax, av, allenst)
% 
%

if nargin<6
    regColor = 0.75*[1 1 1];
end

% indices to highlight
thisid = ['/' num2str(allenst.id(strcmp(allenst.acronym, acr))) '/'];
indAndCh = find(cellfun(@(x)contains(x, thisid), allenst.structure_id_path));

% plotting the slice with highlight
thisSlice = squeeze(av(round(ccfCoord(1)/10),:,:));
im = sliceOutlineWithRegionVec(thisSlice, indAndCh, regColor, ax);

% plotting the neuron
h = plot(ax, ccfCoord(3)/10, ccfCoord(2)/10, 'ro', 'MarkerFaceColor', 'r');
