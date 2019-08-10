% SELECTTSP.M          (universal SELECTion)
%
% This function performs universal selection. The function handles
% multiple populations and calls the low level selection function
% for the actual selection process.

%
% Syntax:  SelCh = select(SEL_F, Chrom, FitnV, GGAP, SUBPOP)
%
% Input parameters:
%    SEL_F     - Name of the selection function
%    Chrom     - Matrix containing the individuals (parents) of the current
%                population. Each row corresponds to one individual.
%    ObjV      - Column vector containing the raw cost/optimality values of the
%                individuals in the population.
%    GGAP      - (optional) Rate of individuals to be selected
%                if omitted 1.0 is assumed
%    SUBPOP    - (optional) Number of subpopulations
%                if omitted 1 subpopulation is assumed
%
% Output parameters:
%    SelCh     - Matrix containing the selected individuals.

% Author:     Hartmut Pohlheim               Daniel Loscos
% History:    10.03.94     file created      28.11.18    file adapted

function SelCh = selectTSP(SEL_F, Chrom, ObjV, GGAP, SUBPOP);

    % Check parameter consistency
       if nargin < 3, error('Not enough input parameter'); end

       % Identify the population size (Nind)
       [NindCh,Nvar] = size(Chrom);
       [NindF,VarF] = size(ObjV);
       if NindCh ~= NindF, error('Chrom and FitnV disagree'); end
       if VarF ~= 1, error('FitnV must be a column vector'); end

       if nargin < 5, SUBPOP = 1; end
       if nargin > 4,
          if isempty(SUBPOP), SUBPOP = 1;
          elseif isnan(SUBPOP), SUBPOP = 1;
          elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar'); end
       end

       if (NindCh/SUBPOP) ~= fix(NindCh/SUBPOP), error('Chrom and SUBPOP disagree'); end
       Nind = NindCh/SUBPOP;  % Compute number of individuals per subpopulation

       if nargin < 4, GGAP = 1; end
       if nargin > 3,
          if isempty(GGAP), GGAP = 1;
          elseif isnan(GGAP), GGAP = 1;
          elseif length(GGAP) ~= 1, error('GGAP must be a scalar');
          elseif (GGAP < 0), error('GGAP must be a scalar bigger than 0'); end
       end

    % Compute number of new individuals (to select)
       NSel=max(floor(Nind*GGAP+.5),2);

    % Select individuals from population
       SelCh = [];
       for irun = 1:SUBPOP,
          ObjVSub = ObjV((irun-1)*Nind+1:irun*Nind);
          ChrIx=feval(SEL_F, ObjVSub, NSel)+(irun-1)*Nind;
          SelCh=[SelCh; Chrom(ChrIx,:)];
       end
end

function NewChrIx = ranksus(ObjV,Nsel);
   
%assign fitness values to entire population
   FitnV = ranking(ObjV);

% Identify the population size (Nind)
   [Nind,ans] = size(FitnV);

% Perform stochastic universal sampling
   cumfit = cumsum(FitnV);
   trials = cumfit(Nind) / Nsel * (rand + (0:Nsel-1)');
   Mf = cumfit(:, ones(1, Nsel));
   Mt = trials(:, ones(1, Nind))';
   [NewChrIx, ans] = find(Mt < Mf & [ zeros(1, Nsel); Mf(1:Nind-1, :) ] <= Mt);

% Shuffle new population
   [ans, shuf] = sort(rand(Nsel, 1));
   NewChrIx = NewChrIx(shuf);
end

function NewChrIx = fitness_prop(ObjV,Nsel);
   
  %assign fitness values to entire population
   FitnV = fitness(ObjV);

% Identify the population size (Nind)
   [Nind,ans] = size(FitnV);

% Perform a random sampling
   cumfit = cumsum(FitnV);
   trials = sort(cumfit(Nind) * rand(Nsel,1));
   Mf = cumfit(:, ones(1, Nsel));
   Mt = trials(:, ones(1, Nind))';
   [NewChrIx, ans] = find(Mt < Mf & [ zeros(1, Nsel); Mf(1:Nind-1, :) ] <= Mt);

% Shuffle new population
   [ans, shuf] = sort(rand(Nsel, 1));
   NewChrIx = NewChrIx(shuf);
end

function NewChrIx = tournament3(ObjV,Nsel);
   
  %assign fitness values to entire population
   FitnV = fitness(ObjV);

% Identify the population size (Nind)
   [Nind,ans] = size(FitnV);

% Perform a random sampling
   cumfit = cumsum(FitnV);
   t = cumfit(Nind) * rand(Nsel,1);
   t = max(t,cumfit(Nind) * rand(Nsel,1));
   trials = sort(max(t,cumfit(Nind) * rand(Nsel,1)));
   Mf = cumfit(:, ones(1, Nsel));
   Mt = trials(:, ones(1, Nind))';
   [NewChrIx, ans] = find(Mt < Mf & [ zeros(1, Nsel); Mf(1:Nind-1, :) ] <= Mt);

% Shuffle new population
   [ans, shuf] = sort(rand(Nsel, 1));
   NewChrIx = NewChrIx(shuf);
end

function NewChrIx = tournament5(ObjV,Nsel);
   
  %assign fitness values to entire population
   FitnV = fitness(ObjV);

% Identify the population size (Nind)
   [Nind,ans] = size(FitnV);

% Perform a random sampling
   cumfit = cumsum(FitnV);
   t = cumfit(Nind) * rand(Nsel,1);
   for i = 1:3
       t = max(t,cumfit(Nind) * rand(Nsel,1));
   end
   trials = sort(max(t,cumfit(Nind) * rand(Nsel,1)));
   Mf = cumfit(:, ones(1, Nsel));
   Mt = trials(:, ones(1, Nind))';
   [NewChrIx, ans] = find(Mt < Mf & [ zeros(1, Nsel); Mf(1:Nind-1, :) ] <= Mt);

% Shuffle new population
   [ans, shuf] = sort(rand(Nsel, 1));
   NewChrIx = NewChrIx(shuf);
end

%returns the fitness values
    function FitV = fitness(ObjV)
        maxcost = max(ObjV);
        FitV = ObjV;
        n = size(ObjV);
        for i = 1:n
           FitV(i) = maxcost - ObjV(i); 
        end
    end