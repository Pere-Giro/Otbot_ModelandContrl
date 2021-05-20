function [ Array ] = from_file( nomfitxer )
%FROM_FILE Summary of this function goes here
%   Detailed explanation goes here

fid=fopen(nomfitxer);

[Array(1:14,1),c]=fscanf(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f',[14 1]);

fclose(fid);

end

