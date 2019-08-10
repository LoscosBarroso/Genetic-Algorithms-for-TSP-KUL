testN = 20;
ResultsRank = 1:1:testN;
ResultsFitnessP = 1:1:testN;
ResultsT3 = 1:1:testN;
ResultsT5 = 1:1:testN;
% load the data sets
datasetslist = dir('datasets/');datasetslist = dir('datasets/');
datasets=cell( size(datasetslist,1)-2,1);datasets=cell( size(datasetslist,1)-2 ,1);
for i=1:size(datasets,1)
    datasets{i} = datasetslist(i+2).name;
end

NIND=200;		% Number of individuals
MAXGEN=80000;		% Maximum no. of generations
MAXTIME = 10;   % Maximum run time
NVAR=26;		% No. of variables
PRECI=1;		% Precision of variables
ELITIST=0.1;    % percentage of the elite population
GGAP=1-ELITIST;		% Generation gap
STOP_PERCENTAGE=.95;    % percentage of equal fitness individuals for stopping
PR_CROSS=.85;     % probability of crossover
PR_MUT=.1;       % probability of mutation
LOCALLOOP=1;      % local loop removal {1 on, 0 off}
EARLY=0;      % early termination criteria {0 'off',1 'No Improvement',2 'Low Diversity',3 'Both'}
STALL_PERCENTAGE=0.20; %percentage of MAXGEN after which no fitness improvement is considered an stall
RESET=0;      % reset on early termination {1 on, 0 off}
CROSSOVER = 'order';  % default crossover operator {'order','edge_recombination'}
MUTATION = 'scramble'; %default mutation operation {'inversion','simple_inversion','scramble','swap_mutation','insert'}
SELECTION = 'ranksus'; %default seleccion operator {'ranksus' (stochastic universal selection),'fitness_prop','tournament3','tournament5'}


dataset = 'xqf131.tsp';
data = load(['datasets/' dataset]);
x=data(:,1);y=data(:,2);
NVAR=size(data,1);
MUTATION = 'simple_inversion';
for i=1:1:testN
    ResultsRank(i) = run_ga(x, y, NIND, MAXGEN, MAXTIME, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, 'ranksus', CROSSOVER, MUTATION, LOCALLOOP, EARLY, STALL_PERCENTAGE, RESET);
    ResultsFitnessP(i) = run_ga(x, y, NIND, MAXGEN, MAXTIME, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, 'fitness_prop', CROSSOVER, MUTATION, LOCALLOOP, EARLY, STALL_PERCENTAGE, RESET);
    ResultsT3(i) = run_ga(x, y, NIND, MAXGEN, MAXTIME, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, 'tournament3', CROSSOVER, MUTATION, LOCALLOOP, EARLY, STALL_PERCENTAGE, RESET);
    ResultsT5(i) = run_ga(x, y, NIND, MAXGEN, MAXTIME, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, 'tournament5', CROSSOVER, MUTATION, LOCALLOOP, EARLY, STALL_PERCENTAGE, RESET);
    [ResultsRank(i) ResultsFitnessP(i) ResultsT3(i) ResultsT5(i)]
    i
end
['$' num2str(min(ResultsRank)) '$ & $' num2str(min(ResultsFitnessP)) '$ & $' num2str(min(ResultsT3)) '$ & $' num2str(min(ResultsT5)) '$']

['$' num2str(mean(ResultsRank)) ' \pm ' num2str(std(ResultsRank)) '$ & $' num2str(mean(ResultsFitnessP)) ' \pm ' num2str(std(ResultsFitnessP)) '$ & $' num2str(mean(ResultsT3)) ' \pm ' num2str(std(ResultsT3)) '$ & $' num2str(mean(ResultsT5)) ' \pm ' num2str(std(ResultsT5)) '$']

H = zeros(4,4);
P = zeros(4,4);
[H(1,1), P(1,1)] = ttest(ResultsRank,ResultsRank);
[H(1,2), P(1,2)] = ttest(ResultsRank,ResultsFitnessP);
[H(1,3), P(1,3)] = ttest(ResultsRank,ResultsT3);
[H(1,4), P(1,4)] = ttest(ResultsRank,ResultsT5);
[H(2,1), P(2,1)] = ttest(ResultsFitnessP,ResultsRank);
[H(2,2), P(2,2)] = ttest(ResultsFitnessP,ResultsFitnessP);
[H(2,3), P(2,3)] = ttest(ResultsFitnessP,ResultsT3);
[H(2,4), P(2,4)] = ttest(ResultsFitnessP,ResultsT5);
[H(3,1), P(3,1)] = ttest(ResultsT3,ResultsRank);
[H(3,2), P(3,2)] = ttest(ResultsT3,ResultsFitnessP);
[H(3,3), P(3,3)] = ttest(ResultsT3,ResultsT3);
[H(3,4), P(3,4)] = ttest(ResultsT3,ResultsT5);
[H(4,1), P(4,1)] = ttest(ResultsT5,ResultsRank);
[H(4,2), P(4,2)] = ttest(ResultsT5,ResultsFitnessP);
[H(4,3), P(4,3)] = ttest(ResultsT5,ResultsT3);
[H(4,4), P(4,4)] = ttest(ResultsT5,ResultsT5);



function realmin = run_ga(x, y, NIND, MAXGEN, MAXTIME, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, SELECTION, CROSSOVER, MUTATION, LOCALLOOP, EARLY, STALL_PERCENTAGE, RESET, ah1)

        stall = ceil(MAXGEN*STALL_PERCENTAGE);
        GGAP = 1 - ELITIST; %Percentage of non-elite individuals.
        mean_fits=zeros(1,MAXGEN+1);
        worst=zeros(1,MAXGEN+1);
        Dist=zeros(NVAR,NVAR);
        for i=1:size(x,1)
            for j=1:size(y,1)
                Dist(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            end
        end
        % initialize population
        Chrom=zeros(NIND,NVAR);
        for row=1:NIND
            Chrom(row,:)=randperm(NVAR);
        end
        gen=0;
        last=0;
        % number of individuals of equal fitness needed to stop
        stopN=ceil(STOP_PERCENTAGE*NIND);
        % evaluate initial population
        ObjV = tspfun(Chrom,Dist);
        realmin = ObjV(1);
        bestpath= Chrom(1,:);
        best=zeros(1,MAXGEN);
        GO = true;
        intime = true;
        % generational loop
        tic
        while gen<MAXGEN && (GO || RESET) && intime
            if(~GO)
                % initialize population
                GO=true;
                Chrom=zeros(NIND,NVAR);
                for row=1:NIND
                   Chrom(row,:)=randperm(NVAR);
                end
                % number of individuals of equal fitness needed to stop
                stopN=ceil(STOP_PERCENTAGE*NIND);
                % evaluate initial population
                ObjV = tspfun(Chrom,Dist);
            end
            sObjV=sort(ObjV);
          	best(gen+1)= min(ObjV);
        	minimum=best(gen+1);
            mean_fits(gen+1)=mean(ObjV);
            worst(gen+1)=max(ObjV);
            for t=1:size(ObjV,1)
                if (ObjV(t)==minimum)
                    break;
                end
            end
            if(minimum < realmin)
                bestpath= Chrom(t,:);
                realmin = minimum;
            end
            
            
            %early termination strategies
            if (EARLY ~=0)
                if (EARLY ~= 1 && sObjV(stopN)-sObjV(1) <= 1e-15)
                    GO = false;
                end    
                if(EARLY~= 2 && gen>(stall+last) && best(gen-stall)==best(gen))
                    GO = false;
                    last = gen;
                end
            end
        	%select individuals for breeding
        	SelCh = selectTSP('tournament5', Chrom, ObjV, GGAP); %selection via stochastic univesal sampling
        	%recombine individuals (crossover)
            SelCh = crossTSP(CROSSOVER,SelCh,PR_CROSS);
            SelCh = mutateTSP(MUTATION,SelCh,PR_MUT);
            %evaluate offspring, call objective function
        	ObjVSel = tspfun(SelCh,Dist);
            %reinsert offspring into population
        	[Chrom, ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
            if(LOCALLOOP)
                Chrom = tsp_ImprovePopulation(NIND, NVAR, Chrom, Dist);
            end
        	%increment generation counter
        	gen=gen+1;
            if(MAXTIME > 0)
               intime = toc < MAXTIME;
            end
        end
        %visualizeTSP(x, y, bestpath, realmin, gen, best, mean_fits, worst, ObjV, NIND);
end


function visualizeTSP(X,Y, Path, TotalDist, gen, best, mean_fits, worst, ObjV, NIND)
        figure();
        plot([X(Path(length(Path))) X(Path(1))],[Y(Path(length(Path))) Y(Path(1))], 'ko-','MarkerFaceColor','Black');
        hold on;
        for k = 1:(length(Path)-1)
            plot([X(Path(k)) X(Path(k+1))],[Y(Path(k)) Y(Path(k+1))], 'ko-','MarkerFaceColor','Black');
        end
    	title(['Best round length: ' num2str(TotalDist)]);
        hold off;
        
        figure();
        plot([0:gen-1],best(1:gen),'r-', [0:gen-1],mean_fits(1:gen),'b-', [0:gen-1],worst(1:gen),'g-');
        xlabel('Generation');
        ylabel('Distance (Min. - Gem. - Max.)');
        hold on
        a = toc;
        title(['Time elapsed (s): ' num2str(a(1))]);
        hold off
        drawnow;
    end