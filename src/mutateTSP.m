% MUTATETSP.M       (MUTATion for TSP high-level function)
%
% This function takes a matrix OldChrom containing the 
% representation of the individuals in the current population,
% mutates the individuals and returns the resulting population.
%
% Syntax:  NewChrom = mutate(MUT_F, OldChrom, MutOpt)
%
% Input parameter:
%    MUT_F     - String containing the name of the mutation function
%    OldChrom  - Matrix containing the chromosomes of the old
%                population. Each line corresponds to one individual.
%    MutOpt    - mutation rate
% Output parameter:
%    NewChrom  - Matrix containing the chromosomes of the population
%                after mutation in the same format as OldChrom.


function NewChrom = mutateTSP(MUT_F, OldChrom, MutOpt)

    % Check parameter consistency
       if nargin < 2,  error('Not enough input parameters'); end

    [rows,~]=size(OldChrom);
    NewChrom=OldChrom;

    for r=1:rows
        if rand<MutOpt
            NewChrom(r,:) = feval(MUT_F, OldChrom(r,:));
        end
    end
end


% low level function for TSP mutation
% inversion : a segment of the tour is inverted
% Representation is path representation
function NewChrom = inversion(OldChrom)

    NewChrom=OldChrom;
    % select two positions in the tour
    k = size(NewChrom,2);
    rndi=zeros(1,2);

    while rndi(1)==rndi(2)
        rndi=rand_int(1,2,[1 k]);
    end
    rndi = sort(rndi);
    d = rndi(2) - rndi(1)+1;
    NewChrom(rndi(1):rndi(2)) = NewChrom(rndi(2):-1:rndi(1));
    NewChrom = circshift(NewChrom,1-rndi(1));
    NewChrom(d+1:k) = circshift(NewChrom(d+1:k), randi([0 k-d]));
end

function NewChrom = simple_inversion(OldChrom)

    NewChrom=OldChrom;

    % select two positions in the tour

    rndi=zeros(1,2);

    while rndi(1)==rndi(2)
        rndi=rand_int(1,2,[1 size(NewChrom,2)]);
    end
    rndi = sort(rndi);
    NewChrom(rndi(1):rndi(2)) = NewChrom(rndi(2):-1:rndi(1));
end

function NewChrom = scramble(OldChrom)
    % select two positions in the tour

    rndi=zeros(1,2);
    rndi(2) = randi([1 size(OldChrom,2)]);
    rndi(1) = randi([2 size(OldChrom,2)]);
    NewChrom = circshift(OldChrom,rndi(2));
    temp = NewChrom(1:rndi(1));
    i = randperm(rndi(1));
    NewChrom(1:rndi(1)) = temp(i);
end

% low level function for TSP mutation
% swap mutation : two random cities in a tour are swapped

function NewChrom = swap_mutation(OldChrom)

    NewChrom=OldChrom;

    % swap two random cities in the tour
    rndi=zeros(1,2);
    while rndi(1)==rndi(2)
        rndi=rand_int(1,2,[1 size(NewChrom,2)]);
    end
    NewChrom([rndi(1) rndi(2)])=NewChrom([rndi(2) rndi(1)]);
end

function NewChrom = insert(OldChrom)
    % swap two random cities in the tour
    rndi=zeros(1,2);
    while rndi(1)==rndi(2)
        rndi=rand_int(1,2,[1 size(OldChrom,2)]);
    end
    NewChrom = circshift(OldChrom,1-rndi(1));
    i = rem(1+size(OldChrom,2)+rndi(2)-rndi(1),size(OldChrom,2));
    NewChrom(2:i)= circshift(NewChrom(2:i),1);
end