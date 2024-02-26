function [data] = interpolateContactLoopstoSolid(domain,solid,resamplingFactor)


data=imresize3(domain,1/resamplingFactor,'nearest');
[nzX, nzY, nzZ] = ind2sub(size(data), find(data > 1));

[zX, zY, zZ] = ind2sub(size(data), find(data <= 1));

nonzero_points = [nzX, nzY, nzZ];
nonzero_values = data(data >1);

F = scatteredInterpolant(nonzero_points, double(nonzero_values), 'natural');

zero_points = [zX, zY, zZ];
interpolated_values=zeros(size(zero_points,1),1);
for i=1:size(zero_points,1)
    if mod(i,1000)==0
        disp(['Interpolating Points: ', num2str(i),' of ',num2str(size(zero_points,1))])
    end
    interpolated_values(i) = F(zero_points(i,:));
end

data(data <= 1) = interpolated_values;
data=imresize3(data,resamplingFactor);
data(~solid)=0;