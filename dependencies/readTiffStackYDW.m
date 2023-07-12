function [data] = readTiffStackYDW(fname,varargin)


info = imfinfo(fname);
num_images = numel(info);
verbose=1;
if nargin>1
    try
        verbose=varargin{1};
    catch
    end
    
    try
        nz=varargin{2};
        if nz<num_images
            num_images=nz;
        end
    catch
    end
    
end
disp(['Reading ', fname])

data=zeros(info(1).Height, info(1).Width, num_images,'uint8');
for k = 1:num_images
    if verbose
        disp(['Reading Slice ', num2str(k), ' of ', num2str(num_images)])
    end
    A = imread(fname, k, 'Info', info);
    A = uint8(A);
    data(:,:,k)=A;
end
end