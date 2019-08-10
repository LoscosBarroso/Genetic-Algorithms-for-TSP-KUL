% low level function for TSP mutation
% reciprocal exchange : two random cities in a tour are swapped
% Representation is an integer specifying which encoding is used
%	1 : adjacency representation
%	2 : path representation
%

function NewChrom = reciprocal_exchange(OldChrom,Representation);

NewChrom=OldChrom;

if Representation==1 
	NewChrom=adj2path(NewChrom);
end

% swap two random cities in the tour
rndi=zeros(1,2);
while rndi(1)==rndi(2)
	rndi=rand_int(1,2,[1 size(NewChrom,2)]);
end

buffer=NewChrom(rndi(1));
NewChrom(rndi(1))=NewChrom(rndi(2));
NewChrom(rndi(2))=buffer;


if Representation==1
	NewChrom=path2adj(NewChrom);
end


% End of function

