function coords = makeSmoothCoords(c)

coords = [];
ii = 1;coordsInd = 1; 
buff = 10;
while ii<size(c,2)
    n = c(2,ii);
    if n>100
        x = c(1,ii+1:ii+n);
        y = c(2,ii+1:ii+n);
        x = [x(end-buff:end) x x(1:buff)]; % buffer makes ends meet
        y = [y(end-buff:end) y y(1:buff)];
        x = smooth(x,25,'loess');
        y = smooth(y,25,'loess');
        x = [x(buff+1:end-buff-1);x(buff+1)];
        y = [y(buff+1:end-buff-1);y(buff+1)];
        coords(coordsInd).x = x;
        coords(coordsInd).y = y;
        coordsInd = coordsInd+1;                                        
        
    end
    ii = ii+n+1;
end