%
% adj2path(Adj)
% function to convert between adjacency and path representation for TSP
% Adj, Path are row vectosr
%

function Path = adj2path(Adj);

	walking_index=1;
	Path=zeros(size(Adj));
	Path(1)=1;
	for t=2:size(Adj,2)
		Path(t)=Adj(walking_index);
		walking_index=Path(t);
	end


% End of function

