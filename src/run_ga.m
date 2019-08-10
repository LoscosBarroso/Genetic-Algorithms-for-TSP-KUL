function run_ga(x, y, NIND, MAXGEN, MAXTIME, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, SELECTION, CROSSOVER, MUTATION, LOCALLOOP, EARLY, STALL_PERCENTAGE, RESET, ah1, ah2, ah3)
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
{NIND MAXGEN MAXTIME NVAR ELITIST STOP_PERCENTAGE PR_CROSS PR_MUT SELECTION CROSSOVER MUTATION LOCALLOOP EARLY STALL_PERCENTAGE RESET}

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
            visualizeTSP(x, y, bestpath, realmin, ah1, gen, best, mean_fits, worst, ah2, ObjV, NIND, ah3);
            
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
end

    function visualizeTSP(X,Y, Path, TotalDist, figNr, gen, best, mean_fits, worst, figNr2, ObjV, NIND, ah3)
        axes(figNr);
        plot([X(Path(length(Path))) X(Path(1))],[Y(Path(length(Path))) Y(Path(1))], 'ko-','MarkerFaceColor','Black');
        hold on;
        for k = 1:(length(Path)-1)
            plot([X(Path(k)) X(Path(k+1))],[Y(Path(k)) Y(Path(k+1))], 'ko-','MarkerFaceColor','Black');
        end
    	title(['Best round length: ' num2str(TotalDist)]);
        hold off;
        axes(figNr2);
        
        plot([0:gen],best(1:gen+1),'r-', [0:gen],mean_fits(1:gen+1),'b-', [0:gen],worst(1:gen+1),'g-');
        xlabel('Generation');
        ylabel('Distance (Min. - Gem. - Max.)');
        hold on
        a = toc;
        title(['Time elapsed (s): ' num2str(a(1))]);
        hold off
        axes(ah3);
        bins = max([1 ceil((max(ObjV) - min(ObjV))/0.3)]);
        limits = get(ah3,'Xlim');
        limit_b = limits(2);
        hist(ObjV, bins);
        xlabel('Distance');
        ylabel('Number');
        limits = get(ah3,'Xlim');
        limit_a = limits(2);

        set(ah3,'Xlim',[0 max([limit_a limit_b])]);
        set(ah3,'Ylim',[0 NIND]);
        drawnow;
    end
    

