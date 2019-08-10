function result = improve_path (ncities,path,Dist)
% improve_path.m
% Author: Mike Matton 
% Date: October 2008
%
% This function improves a single tsp path in path representation by removing 
% local loops up to pathlength 3. 
%
% Input parameters:
%    ncities         - The number of cities
%    path            - The path to improve (path representation)
%    Dist            - Matrix with distances between cities
% Output parameter:
%    result          - The improved path


maxlen =  min(3, ncities / 2) ;

for len = 2:maxlen
    for start = 1:ncities
        stop = mod((start + len - 1 -1 ),ncities) + 1;
        gain = Dist(path(1,mod(start -1 +ncities-1,ncities)+1),path(1,start)) ...
               + Dist(path(1,stop),path(1,mod(stop+1 -1 ,ncities)+1)) ...
               - Dist(path(1,mod(start -1 +ncities-1,ncities)+1),path(1,stop)) ...
               - Dist(path(1,start),path(1,mod(stop+1 -1 ,ncities)+1));
        if ( gain > 0.0 ) 
            path = SwapSubpath(ncities, path, start, len);
        end
    end
end
result = path;

