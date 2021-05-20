function [ Array ] = to_file( nomfitxer,Geom_Vector, savePath )
%TO_FILE Summary of this function goes here
%   Detailed explanation goes here

fid=fopen([savePath strcat(nomfitxer,'.txt')],'w');
for i=1:max(size(Geom_Vector))
    fprintf(fid,num2str(Geom_Vector(i),16));
    fprintf(fid,'\t\t');
end
fprintf(fid,'\n');

fclose(fid);

Array=Geom_Vector;

end

