% crossTSP.M       (RECOMBINation high-level function)
%
% This function performs recombination between pairs of individuals
% and returns the new individuals after mating. The function handles
% multiple populations and calls the low-level recombination function
% for the actual recombination process.
%
% Syntax:  NewChrom = crossTSP(REC_F, OldChrom, RecOpt, SUBPOP)
%
% Input parameters:
%    REC_F     - String containing the name of the recombination or
%                crossover function
%    Chrom     - Matrix containing the chromosomes of the old
%                population. Each line corresponds to one individual
%    RecOpt    - (optional) Scalar containing the probability of 
%                recombination/crossover occurring between pairs
%                of individuals.
%                if omitted or NaN, 1 is assumed
%    SUBPOP    - (optional) Number of subpopulations
%                if omitted or NaN, 1 subpopulation is assumed
%
% Output parameter:
%    NewChrom  - Matrix containing the chromosomes of the population
%                after recombination in the same format as OldChrom.

%  Author:    Hartmut Pohlheim
%  History:   18.03.94     file created 24.11.18 file adapted and renamed 


function NewChrom = crossTSP(REC_F, Chrom, RecOpt, SUBPOP)
% Check parameter consistency
   if nargin < 2, error('Not enough input parameter'); end
   % Identify the population size (Nind)
   [Nind,Nvar] = size(Chrom);
   if nargin < 4, SUBPOP = 1; end
   if nargin > 3,
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar'); end
   end
   if (Nind/SUBPOP) ~= fix(Nind/SUBPOP), error('Chrom and SUBPOP disagree'); end
   Nind = Nind/SUBPOP;  % Compute number of individuals per subpopulation
   if nargin < 3, RecOpt = 0.7; end
   if nargin > 2,
      if isempty(RecOpt), RecOpt = 0.7;
      elseif isnan(RecOpt), RecOpt = 0.7;
      elseif length(RecOpt) ~= 1, error('RecOpt must be a scalar');
      elseif (RecOpt < 0 | RecOpt > 1), error('RecOpt must be a scalar in [0, 1]'); end
   end
% Select individuals of one subpopulation and call low level function
   NewChrom = [];
   for irun = 1:SUBPOP,
      ChromSub = Chrom((irun-1)*Nind+1:irun*Nind,:);  
      NewChromSub = feval(REC_F, ChromSub, RecOpt);
      NewChrom=[NewChrom; NewChromSub];
   end
end 

function NewChrom = order(OldChrom, XOVR)
    if nargin < 2, XOVR = NaN; end
    [rows,cols]=size(OldChrom);
       maxrows=rows;
       if rem(rows,2)~=0
           maxrows=maxrows-1;
       end
       for row=1:2:maxrows
        % crossover of the two chromosomes
        % results in 2 offsprings
        if rand<XOVR			% recombine with a given probability
            k = randi([0 cols-1]);
            OldChrom(row,:) = circshift(OldChrom(row,:),k);
            OldChrom(row+1,:) = circshift(OldChrom(row+1,:),k);
            k = randi([1 floor(cols/2)]);
            NewChrom(row,:) = cross_order([OldChrom(row,:);OldChrom(row+1,:)],k);
            NewChrom(row+1,:) = cross_order([OldChrom(row+1,:);OldChrom(row,:)],k);
        else
            NewChrom(row,:)=OldChrom(row,:);
            NewChrom(row+1,:)=OldChrom(row+1,:);
        end
       end
       if rem(rows,2)~=0
           NewChrom(rows,:)=OldChrom(rows,:);
       end
       
        function Offspring=cross_order(Parents,k)
            cols=size(Parents,2);
            Offspring=zeros(1,cols);
            Offspring(1:k) = Parents(1,1:k);
            c = k+1;
            while k < cols
                k = k+1;
                while(ismember(Parents(2,c), Offspring(1:k))) 
                    c = rem(c,cols)+1;
                end
                Offspring(k) = Parents(2,c);
                c = rem(c,cols)+1;
            end
        end
end

function NewChrom = edge_recombination(OldChrom, XOVR)
    if nargin < 2, XOVR = NaN; end
    [rows,cols]=size(OldChrom);
       maxrows=rows;
       if rem(rows,2)~=0
           maxrows=maxrows-1;
       end
       for row=1:2:maxrows
        % crossover of the two chromosomes
        % results in 2 offsprings
        if rand<XOVR			% recombine with a given probability
            NewChrom(row,:) =cross_edge_recombination([OldChrom(row,:);OldChrom(row+1,:)]);
            NewChrom(row+1,:)=cross_edge_recombination([OldChrom(row+1,:);OldChrom(row,:)]);
        else
            NewChrom(row,:)=OldChrom(row,:);
            NewChrom(row+1,:)=OldChrom(row+1,:);
        end
       end
       if rem(rows,2)~=0
           NewChrom(rows,:)=OldChrom(rows,:);
       end
        function Offspring=cross_edge_recombination(Parents)
            cols=size(Parents,2);
            Offspring=zeros(1,cols);
            banned = zeros(1,cols);
            edges = zeros(cols,4);
            edges(Parents(1,1),[1 2]) = Parents(1,[cols 2]);
            edges(Parents(2,1),[3 4]) = Parents(2,[cols 2]);
            for i = 2:cols-1
                edges(Parents(1,i),[1 2]) = Parents(1,[i-1 i+1]);
                edges(Parents(2,i),[3 4]) = Parents(2,[i-1 i+1]);
            end
            edges(Parents(1,cols),[1 2]) = Parents(1,[cols-1 1]);
            edges(Parents(2,cols),[3 4]) = Parents(2,[cols-1 1]);
            entry = randi([1 cols]);
            Offspring(1) = entry;
            banned(entry) = 1;
            [entry, edges] = next(edges,entry,cols,banned);
            for i = 2:cols-1
                Offspring(i) = entry;
                banned(entry) = 1;
                [entry, edges] = next(edges,entry,cols,banned);
            end
            [~, entry] = min(banned);
            Offspring(cols) = entry;
        end
        
    function [next, edges] = next(edges, entry, cols, banned)
        count = zeros(1,cols);
        for i = 1:cols
            if banned(i)
                count(i) = -1;
            else
                rep = zeros(1,4);
                k=1;
                for j = 1:4
                    if edges(i,j) == -1
                    elseif edges(i,j) == entry
                        edges(i,j) = -1;
                    else
                        if ~ismember(edges(i,j),rep)
                            rep(k) = edges(i,j);
                            k = k+1;
                            count(i) = count(i)+1;
                        end
                   end    
                end
            end   
        end
        countel=zeros(1,cols);
        for j = 1:4
            if (edges(entry,j) ~= -1)
                countel(edges(entry,j)) = countel(edges(entry,j)) + 1;
            end
        end
        [repelem, next] = max(countel);
        if repelem == 0
            [~, next] = min(banned);
        end
        if repelem == 1
            len = zeros(1,4);
            for j = 1:4
                if (edges(entry,j) == -1)
                	len(j) = 5;
                else
                    len(j) = count(edges(entry,j)); 
                end
            end
            [~, next] = min(len);
            next = edges(entry,next);
        end    
    end
end

