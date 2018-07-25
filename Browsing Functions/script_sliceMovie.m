
%% Create the plots
szAP = size(tv,1);
szDV = size(tv,2);
szLR = size(tv,3);

fh = figure; set(fh, 'Position', [-1896         -72        1479        1024]);

subtightplot(2,2,1); 
bregma = allenCCFbregma();
isBrain = av>1; % >0 for original av, >1 for by_index
gridIn3D(double(isBrain), 0.5, 50, bregma);
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
view([-30    25]);

sliceCoronal = fill3([0 0 0 0], [0 szLR szLR 0], [0 0 szDV szDV], 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
sliceHoriz = fill3([0 szAP szAP 0], [0 0 szLR szLR], [0 0 0 0], 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
sliceSag = fill3([0 szAP szAP 0], [0 0 0 0], [0 0 szDV szDV], 'g', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);

subtightplot(2,2,2); 
imCoronal = plotTVslice(squeeze(tv(1,:,:)));

subtightplot(2,2,3); 
imHoriz = plotTVslice(squeeze(tv(:,1,:)));

subtightplot(2,2,4); 
imSag = plotTVslice(squeeze(tv(:,:,1)));

%% Page through as a movie

makeVideo = true;

if makeVideo
    WriterObj = VideoWriter('allenCCFsliceMovie');
    WriterObj.FrameRate=30;
    open(WriterObj);
end

nFr = 440;
apInds = round(linspace(1, szAP, nFr));
dvInds = round(linspace(1, szDV, nFr));
lrInds = round(linspace(1, szLR, nFr));

for f = 1:nFr
    sliceCoronal.Vertices(:,1) = apInds(f);
    sliceHoriz.Vertices(:,3) = dvInds(f);
    sliceSag.Vertices(:,2) = lrInds(f);
    
    set(imCoronal, 'CData', squeeze(tv(apInds(f), :,:)));
    set(imHoriz, 'CData', squeeze(tv(:,dvInds(f),:)));
    set(imSag, 'CData', squeeze(tv(:,:,lrInds(f)))');
    
    drawnow;
    if makeVideo
        frame = getframe(fh);
        writeVideo(WriterObj,frame);
    end
    
end

if makeVideo
    close(WriterObj);
end
