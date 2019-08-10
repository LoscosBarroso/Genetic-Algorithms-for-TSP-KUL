% tsp_ImprovePopulation.m
% Author: Mike Matton
% 
% This function improves a tsp population by removing local loops from
% each individual.
%
% Syntax: improvedPopulation = tsp_ImprovePopulation(popsize, ncities, pop, improve, dists)
%
% Input parameters:
%   popsize           - The population size
%   ncities           - the number of cities
%   pop               - the current population (path representation)
%   improve           - Improve the population (0 = no improvement, <>0 = improvement)
%   dists             - distance matrix with distances between the cities
%
% Output parameter:
%   improvedPopulation  - the new population after loop removal (if improve
%                          <> 0, else the unchanged population).

function newpop = tsp_ImprovePopulation(popsize, ncities, pop, dists)
    
    for i=1:popsize  
        pop(i,:) = improve_path(ncities, pop(i,:),dists);
    end
    newpop = pop;
end

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
end

function result = SwapSubpath(ncities, path, start, length)
    i = start;
    j = mod(start-1 + length - 1, ncities)+1;

    while ( length > 1 ) 
       temp = path(i);
       path(i) = path(j);
       path(j) = temp;
       length = length - 2;
       if ( ++i >= ncities ) 
           i = 0;
       end
       if ( --j <  0 ) 
           j = ncities - 1;
       end
    end

    result = path;
end