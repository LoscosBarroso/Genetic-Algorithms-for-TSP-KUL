function run_ga(x, y, NIND, MAXGEN, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, CROSSOVER, LOCALLOOP, EARLY, STALL_PERCENTAGE, RESET, ah1, ah2, ah3)
% usage: run_ga(x, y, 
%               NIND, MAXGEN, NVAR, 
%               ELITIST, STOP_PERCENTAGE, 
%               PR_CROSS, PR_MUT, CROSSOVER, 
%               ah1, ah2, ah3)
%
%
% x, y: coordinates of the cities
% NIND: number of individuals
% MAXGEN: maximal number of generations
% ELITIST: percentage of elite population
% STOP_PERCENTAGE: percentage of equal fitness (stop criterium)
% PR_CROSS: probability for crossover
% PR_MUT: probability for mutation
% CROSSOVER: the crossover operator
% calculate distance matrix between each pair of cities
% ah1, ah2, ah3: axes handles to visualise tsp
{NIND MAXGEN NVAR ELITIST STOP_PERCENTAGE PR_CROSS PR_MUT CROSSOVER LOCALLOOP EARLY STALL_PERCENTAGE RESET}

        stall = ceil(MAXGEN*STALL_PERCENTAGE);
        GGAP = 1 - ELITIST;
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
        	Chrom(row,:)=path2adj(randperm(NVAR));
            %Chrom(row,:)=randperm(NVAR);
        end
        gen=0;
        last=0;
        % number of individuals of equal fitness needed to stop
        stopN=ceil(STOP_PERCENTAGE*NIND);
        % evaluate initial population
        ObjV = tspfun(Chrom,Dist);
        realmin = ObjV(1);
        best=zeros(1,MAXGEN);
        GO = true;
        % generational loop
        while gen<MAXGEN && (GO || RESET)
            if(~GO)
                % initialize population
                GO=true;
                Chrom=zeros(NIND,NVAR);
                for row=1:NIND
                   Chrom(row,:)=path2adj(randperm(NVAR));
                   %Chrom(row,:)=randperm(NVAR);
                end
                % number of individuals of equal fitness needed to stop
                stopN=ceil(STOP_PERCENTAGE*NIND);
                % evaluate initial population
                ObjV = tspfun(Chrom,Dist);
            end
            sObjV=sort(ObjV);
          	best(gen+1)=min(ObjV);
        	minimum=best(gen+1);
            mean_fits(gen+1)=mean(ObjV);
            worst(gen+1)=max(ObjV);
            for t=1:size(ObjV,1)
                if (ObjV(t)==minimum)
                    break;
                end
            end
            if(minimum < realmin)
                bestpath=adj2path(Chrom(t,:));
                realmin = minimum;
            end
            visualizeTSP(x,y,bestpath, realmin, ah1, gen, best, mean_fits, worst, ah2, ObjV, NIND, ah3);
            
            if (EARLY ~=0)
                if (EARLY ~= 1 && sObjV(stopN)-sObjV(1) <= 1e-15)
                    GO = false;
                end    
                if(EARLY~= 2 && gen>(stall+last) && best(gen-stall)==best(gen))
                    GO = false;
                    last = gen;
                    gen
                end
            end
        	%assign fitness values to entire population
        	FitnV=ranking(ObjV);
        	%select individuals for breeding
        	SelCh=select('sus', Chrom, FitnV, GGAP);
        	%recombine individuals (crossover)
            SelCh = recombin(CROSSOVER,SelCh,PR_CROSS);
            SelCh=mutateTSP('inversion',SelCh,PR_MUT);
            %evaluate offspring, call objective function
        	ObjVSel = tspfun(SelCh,Dist);
            %reinsert offspring into population
        	[Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
            Chrom = tsp_ImprovePopulation(NIND, NVAR, Chrom,LOCALLOOP,Dist);
        	%increment generation counter
        	gen=gen+1;            
        end
end
