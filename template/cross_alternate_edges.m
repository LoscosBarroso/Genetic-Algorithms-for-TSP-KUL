% low level function for calculating an offspring
% given 2 parent in the Parents - agrument
% Parents is a matrix with 2 rows, each row
% represent the genocode of the parent
%
% Tim.Volodine@cs.kuleuven.ac.be
% Tim.Pillards@cs.kuleuven.ac.be
%

function Offspring=cross_alternate_edges(Parents);
	cols=size(Parents,2);
	Offspring=zeros(1,cols);
    [sx,InverseParents1]=sort(Parents(1,:));
    [sx,InverseParents2]=sort(Parents(2,:));
    InverseParents=[InverseParents1;InverseParents2];
    % InverseParents is the same tour but in reversed direction
    % to easily find both edges in a city
    AllParents=[Parents;InverseParents];
    % AllParents contains both directed edges for every city
    % the last two rows contain the inversed tours
	
	start_index=rand_int(1,1,[1 cols]);
%	start_index=1;
	walking_index=start_index; % select a random seed edge
	parentNr=1;
	visited_list=zeros(1,cols);
	visited_list(walking_index)=1;
	cities=1;
	while cities<cols
%		fprintf('%i',cities');
        direction=rand_int(1,1,[0,1]); % randomly chose the direction
		new_city=AllParents(parentNr+2*direction,walking_index);
        % if direction=1, use inverse direction (rows 3 and 4), hence the
        % +2*direction
        
        % if city already visited, try other direction
        if visited_list(new_city)==1
            direction=1-direction; % switch direction: 0->1, 1->0
            new_city=AllParents(parentNr+2*direction,walking_index);
        end
        
        %if other direction also visited, visit a random city
		if visited_list(new_city)==1		% & cities<cols-1
			% resolve the loop
			allowed_list=zeros(1,cols-cities);
			z=1;
			for t=1:size(visited_list,2)
				if visited_list(t)==0
					allowed_list(z)=t;
					z=z+1;
				end
			end
			if z-1~=cols-cities
				warning('something went wrong');
			end
			new_city=allowed_list(rand_int(1,1,[1 size(allowed_list,2)]));
		end
		
	  	Offspring(walking_index)=new_city;
		walking_index=new_city;
		cities=cities+1;
		visited_list(walking_index)=1;
		parentNr=3-parentNr;  % switch the parents (1->2 or 2->1)
	end
	Offspring(walking_index)=start_index;
% end function
