function [ P ] = readpoints( filename )
%READPOINTS Summary of this function goes here
%   Detailed explanation goes here

file = fopen(filename, 'r');
N = fscanf(file, '%d', 1);
P = fscanf(file, '%f%f%f', [3,N]);
fclose(file);

end

