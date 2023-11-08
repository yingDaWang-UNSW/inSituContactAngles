function [] = writeTiffStackYDW(fname,domain)

domain=uint8(domain);
for m = 1:size(domain,3)
    disp(['Writing Slice ', num2str(m), ' of ', num2str(size(domain,3))])
    tiff = domain(:,:,m);
    if m==1
        imwrite(tiff,fname)
    else
        imwrite(tiff,fname,'WriteMode', 'append')
    end
end