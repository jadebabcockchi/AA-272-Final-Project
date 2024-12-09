function tropo = tropomodel(zd, el)
    sz = size(el,1);
    tropo = ones(sz,1) * NaN;
    for i = 1:sz
        if abs(el(i)) >= 0
            tropo(i) = zd/sind(el(i));
        end
    end
    %tropo = 1./sind(el); 
end