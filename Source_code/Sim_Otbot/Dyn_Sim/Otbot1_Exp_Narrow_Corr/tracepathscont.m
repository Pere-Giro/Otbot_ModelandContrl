function [] = tracepathscont(Gpvec,color)
%TRACEPATHS Summary of this function goes here
%   This function is designed to plots the trajectory path of a given point
%   as a continuos line

hold on
plot(Gpvec(1,:),Gpvec(2,:),strcat('-.',color))
hold off

end

