
function plotBrainOutlinesByAxis(axH,av)

for v = 1:3
    
    axes(axH(v)); box off; hold on
    
    sliceIm = squeeze(max(av,[],v)>1);
    
    c = contourc(double(sliceIm), [0.5 0.5]);
    coordsReg = makeSmoothCoords(c);
    for cidx = 1:numel(coordsReg)
        h = plot(axH(v), coordsReg(cidx).x,coordsReg(cidx).y, 'k', 'LineWidth', 2.0); hold on;
        
    end
    
    set(axH(v), 'YDir', 'reverse'); 
    axis equal

end

view([90 -90]); % sagittal should be flipped
set(axH(v), 'XDir', 'reverse'); 
