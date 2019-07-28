

function [im, h] = sliceOutlineWithRegion(avSlice, regInd, regColor, ax)

if nargin<2
    regInd = [];
    regColor = [];
end
if nargin<4
    ax = [];
end

diffilt = -1/8*ones(3,3);
diffilt(2,2) = 1;
thisd = conv2(avSlice, diffilt, 'same');

thisbw = abs(sign(thisd));
se = strel('square', 2);
im = imerode(thisbw, se); % im is size of avSlice, 0 for bckg and 1 for borders
im(end,:) = 0; im(:,end) = 0; % these were ones 

if ~isempty(regInd) % specified a region to highlight
    highlightPts = avSlice==regInd & im==0;
    im = uint8(~im*255);
    im1 = im; im1(highlightPts)=uint8(regColor(1)*255);
    im2 = im; im2(highlightPts)=uint8(regColor(2)*255);
    im3 = im; im3(highlightPts)=uint8(regColor(3)*255);
    im = cat(3, im1, im2, im3);    
    
else
    im = repmat(uint8(~im*255), [1,1,3]);
end

if ~isempty(ax)
    h = image(ax, im);
    axis image
    axis off
end