
function [f,h, whichSlice] = sagittalSlices(tv, ccfCoords, f)

if nargin<3
    f = figure;
    p = get(f, 'Position');
    set(f, 'Position',   [p(1) p(2)      1597         700]);
    set(f, 'Color', 'k');
    slicesExist = false;
else 
    slicesExist = true;
end

slicePoints = 0.25:0.25:3.0; %mm
bregma = allenCCFbregma(); %voxels
sliceCCF = -slicePoints*1000+bregma(3)*10;

cc = ccfCoords;

r = randn(size(cc,1),2)*45;

[~,whichSlice] = min(abs(cc(:,3)-sliceCCF),[],2);

h = [];

apCoords = (0:size(tv,1)-1)*10;
dvCoords = (0:size(tv,2)-1)*10;

inclAP = apCoords>1000 & apCoords<11000; 
inclDV = dvCoords<6000;

for s = 1:numel(sliceCCF)
    subtightplot(3,4,s);
    if ~slicesExist
        imagesc(apCoords(inclAP), dvCoords(inclDV), ...
            tv(inclAP,inclDV,round(sliceCCF(s)/10))');
        colormap gray
        axis image
        caxis([0 400]);
        axis off
    end
    
%     title(sprintf('ML = %.2f', slicePoints(s)));
    
    hold on;
    
    inclN = whichSlice==s & cc(:,3)>0;
            
    h(s) = scatter(cc(inclN,1)+r(inclN,1),cc(inclN,2)+r(inclN,2),4,'r','filled');
end