function [x y] = input_cities(ncities)
        % get the input cities
        fg1 = figure(1);clf;
        %subplot(2,2,2);
        axis([0 1 0 1]);
        title(ncities);
        hold on;
        x=zeros(ncities,1);y=zeros(ncities,1);
        for v=1:ncities
            [xi,yi]=ginput(1);
            x(v)=xi;
            y(v)=yi;
            plot(xi,yi, 'ko','MarkerFaceColor','Black');
            title(ncities-v);
        end
        hold off;
        set(fg1, 'Visible', 'off');
end