    function visualizeTSP(X,Y, Path, TotalDist, figNr, gen, best, mean_fits, worst, figNr2, ObjV, NIND, ah3)
        axes(figNr);
        plot(X(Path),Y(Path), 'ko-','MarkerFaceColor','Black');
        hold on;
        plot([X(Path(length(Path))) X(Path(1))],[Y(Path(length(Path))) Y(Path(1))], 'ko-','MarkerFaceColor','Black');
    	title(['Beste rondrit lengte: ' num2str(TotalDist)]);
        hold off;
        axes(figNr2);
        plot([0:gen],best(1:gen+1),'r-', [0:gen],mean_fits(1:gen+1),'b-', [0:gen],worst(1:gen+1),'g-');
        xlabel('Generation');
        ylabel('Distance (Min. - Gem. - Max.)');       
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