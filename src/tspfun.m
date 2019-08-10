%
% ObjVal = tspfun(Phen, Dist)
% Implementation of the TSP fitness function
%	Phen contains the phenocode of the matrix coded in path
%	representation
%	Dist is the matrix with precalculated distances between each pair of cities
%	ObjVal is a vector with the fitness values for each candidate tour (=each row of Phen)
%

function ObjVal = tspfun(Phen, Dist)
    dims = size(Phen);
    ObjVal = zeros(dims(1),1);
    for k = 1:dims(1)
        ObjVal(k)=Dist(Phen(k,1),Phen(k,dims(2)));
        for t=2:dims(2)
            ObjVal(k)=ObjVal(k)+Dist(Phen(k,t),Phen(k,t-1));
        end
    end   
end

% End of function

