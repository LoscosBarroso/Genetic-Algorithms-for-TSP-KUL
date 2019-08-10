function tspgui()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NIND=50;		% Number of individuals
MAXGEN=100;		% Maximum no. of generations
MAXTIME = -1;   % Maximum run time
NVAR=26;		% No. of variables
PRECI=1;		% Precision of variables
ELITIST=0.05;    % percentage of the elite population
GGAP=1-ELITIST;		% Generation gap
STOP_PERCENTAGE=.95;    % percentage of equal fitness individuals for stopping
PR_CROSS=.95;     % probability of crossover
PR_MUT=.05;       % probability of mutation
LOCALLOOP=0;      % local loop removal
EARLY=0;      % early termination criteria
STALL_PERCENTAGE=0.20; %percentage of MAXGEN after which no fitness improvement is considered an stall
RESET=0;      % reset on early termination
CROSSOVER = 'order';  % default crossover operator
MUTATION = 'simple_inversion'; %default mutation operation
SELECTION = 'ranksus'; %default seleccion operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read an existing population
% 1 -- to use the input file specified by the filename
% 0 -- click the cities yourself, which will be saved in the file called
%USE_FILE=0;
%FILENAME='data/cities.xy';
%if (USE_FILE==0)
%    % get the input cities
%    fg1 = figure(1);clf;
%    %subplot(2,2,2);
%    axis([0 1 0 1]);
%    title(NVAR);
%    hold on;
%    x=zeros(NVAR,1);y=zeros(NVAR,1);
%    for v=1:NVAR
%        [xi,yi]=ginput(1);
%        x(v)=xi;
%        y(v)=yi;
%        plot(xi,yi, 'ko','MarkerFaceColor','Black');
%        title(NVAR-v);
%    end
%    hold off;
%    set(fg1, 'Visible', 'off');
%    dlmwrite(FILENAME,[x y],'\t');
%else
%    XY=dlmread(FILENAME,'\t');
%    x=XY(:,1);
%    y=XY(:,2);
%end

% load the data sets
datasetslist = dir('datasets/');datasetslist = dir('datasets/');
datasets=cell( size(datasetslist,1)-2,1);datasets=cell( size(datasetslist,1)-2 ,1);
for i=1:size(datasets,1)
    datasets{i} = datasetslist(i+2).name;
end

% start with first dataset
data = load(['datasets/' datasets{1}]);
x=data(:,1);y=data(:,2);
NVAR=size(data,1);

datasets

% initialise the user interface
%fh = figure('Visible','off','Name','TSP Tool','Position',[0,0,1024,768]);
fh = figure('Visible','off','Name','TSP Tool','Position',[0,0,1500,950]);
%ah1 = axes('Parent',fh,'Position',[.1 .55 .4 .4]);
ah1 = axes('Parent',fh,'Position',[.05 .53 .42 .42]);
plot(x,y,'ko')
%ah2 = axes('Parent',fh,'Position',[.55 .55 .4 .4]);
ah2 = axes('Parent',fh,'Position',[.53 .53 .42 .42]);
axes(ah2);
xlabel('Generation');
ylabel('Distance (Min. - Gem. - Max.)');
%ah3 = axes('Parent',fh,'Position',[.1 .1 .4 .4]);
ah3 = axes('Parent',fh,'Position',[.05 .05 .42 .42]);
axes(ah3);
title('Histogram');
xlabel('Distance');
ylabel('Number');

%ph = uipanel('Parent',fh,'Title','Settings','Position',[.55 .05 .45 .45]);
ph = uipanel('Parent',fh,'Title','Settings','Position',[.53 .05 .42 .42]);
datasetpopuptxt = uicontrol(ph,'Style','text','String','Dataset','Position',[0 340 90 20]);
datasetpopup = uicontrol(ph,'Style','popupmenu','String',datasets,'Value',1,'Position',[90 340 130 20],'Callback',@datasetpopup_Callback);
llooppopuptxt = uicontrol(ph,'Style','text','String','Loop Detection','Position',[220 340 130 20]);
llooppopup = uicontrol(ph,'Style','popupmenu','String',{'off','on'},'Value',1,'Position',[340 340 50 20],'Callback',@llooppopup_Callback);
earlytxt = uicontrol(ph,'Style','text','String','Early Stop','Position',[400 340 80 20]);
earlytermination = uicontrol(ph,'Style','popupmenu','String',{'off','No Improvement','Low Diversity', 'Both'},'Value',1,'Position',[480 340 100 20],'Callback',@earlytermination_Callback);
resettxt = uicontrol(ph,'Style','text','String','Reset (Early Stop)','Position',[390 310 130 20]);
reset = uicontrol(ph,'Style','popupmenu','String',{'off','on'},'Value',1,'Position',[530 310 50 20],'Callback',@reset_Callback);
maxtimetxt = uicontrol(ph,'Style','text','String','Max. time','Position',[0 310 90 20]);
maxtime = uicontrol(ph,'Style','popupmenu','String',{'Unlimited','1 min','2 min','5 min','10 min'},'Value',1,'Position',[90 310 130 20],'Callback',@maxtime_Callback);
ncitiesslidertxt = uicontrol(ph,'Style','text','String','# Cities','Position',[0 230 130 20]);
%ncitiesslider = uicontrol(ph,'Style','slider','Max',128,'Min',4,'Value',NVAR,'Sliderstep',[0.012 0.05],'Position',[130 230 150 20],'Callback',@ncitiesslider_Callback);
ncitiessliderv = uicontrol(ph,'Style','text','String',NVAR,'Position',[280 230 50 20]);
nindslidertxt = uicontrol(ph,'Style','text','String','# Individuals','Position',[0 200 130 20]);
nindslider = uicontrol(ph,'Style','slider','Max',1000,'Min',10,'Value',NIND,'Sliderstep',[0.001 0.05],'Position',[130 200 150 20],'Callback',@nindslider_Callback);
nindsliderv = uicontrol(ph,'Style','text','String',NIND,'Position',[280 200 50 20]);
genslidertxt = uicontrol(ph,'Style','text','String','# Generations','Position',[0 170 130 20]);
genslider = uicontrol(ph,'Style','slider','Max',1000,'Min',10,'Value',MAXGEN,'Sliderstep',[0.001 0.05],'Position',[130 170 150 20],'Callback',@genslider_Callback);
gensliderv = uicontrol(ph,'Style','text','String',MAXGEN,'Position',[280 170 50 20]);
mutslidertxt = uicontrol(ph,'Style','text','String','Pr. Mutation','Position',[0 140 130 20]);
mutslider = uicontrol(ph,'Style','slider','Max',100,'Min',0,'Value',round(PR_MUT*100),'Sliderstep',[0.01 0.05],'Position',[130 140 150 20],'Callback',@mutslider_Callback);
mutsliderv = uicontrol(ph,'Style','text','String',round(PR_MUT*100),'Position',[280 140 50 20]);
crossslidertxt = uicontrol(ph,'Style','text','String','Pr. Crossover','Position',[0 110 130 20]);
crossslider = uicontrol(ph,'Style','slider','Max',100,'Min',0,'Value',round(PR_CROSS*100),'Sliderstep',[0.01 0.05],'Position',[130 110 150 20],'Callback',@crossslider_Callback);
crosssliderv = uicontrol(ph,'Style','text','String',round(PR_CROSS*100),'Position',[280 110 50 20]);
elitslidertxt = uicontrol(ph,'Style','text','String','% elite','Position',[0 80 130 20]);
elitslider = uicontrol(ph,'Style','slider','Max',100,'Min',0,'Value',round(ELITIST*100),'Sliderstep',[0.01 0.05],'Position',[130 80 150 20],'Callback',@elitslider_Callback);
elitsliderv = uicontrol(ph,'Style','text','String',round(ELITIST*100),'Position',[280 80 50 20]);
selecttxt = uicontrol(ph,'Style','text','String','Selection','Position',[320 230 130 20]);
select = uicontrol(ph,'Style','popupmenu','String',{'stochastic univ. s.','fitness proportional','tournament 3','tournament 5'},'Value',1,'Position',[450 230 130 20],'Callback',@selection_Callback);
crossovertxt = uicontrol(ph,'Style','text','String','Crossover','Position',[320 200 130 20]);
crossover = uicontrol(ph,'Style','popupmenu', 'String',{'order','edge_recombination'}, 'Value',1,'Position',[450 200 130 20],'Callback',@crossover_Callback);
mutationtxt = uicontrol(ph,'Style','text','String','Mutation','Position',[320 140 130 20]);
mutation = uicontrol(ph,'Style','popupmenu', 'String',{'inversion','simple_inversion','scramble','swap_mutation','insert'}, 'Value',1,'Position',[450 140 130 20],'Callback',@mutation_Callback);
%inputbutton = uicontrol(ph,'Style','pushbutton','String','Input','Position',[55 10 70 30],'Callback',@inputbutton_Callback);
runbutton = uicontrol(ph,'Style','pushbutton','String','START','Position',[0 10 50 30],'Callback',@runbutton_Callback);

set(fh,'Visible','on');


    function datasetpopup_Callback(hObject,eventdata)
        dataset_value = get(hObject,'Value');
        dataset = datasets{dataset_value};
        % load the dataset
        data = load(['datasets/' dataset]);
        x=data(:,1);y=data(:,2);
        NVAR=size(data,1); 
        set(ncitiessliderv,'String',size(data,1));
        axes(ah1);
        plot(x,y,'ko') 
    end
    function llooppopup_Callback(hObject,eventdata)
        LOCALLOOP = get(hObject,'Value') -1;
    end
    function earlytermination_Callback(hObject,eventdata)
        EARLY = get(hObject,'Value') -1;
    end
    function reset_Callback(hObject,eventdata)
        RESET = get(hObject,'Value') -1;
    end
    function maxtime_Callback(hObject,eventdata)
        times = [60,120,300,600];
        MAXTIME = get(hObject,'Value') -1;
        if(MAXTIME > 0)
            MAXTIME = times(MAXTIME);
        end
    end
    function ncitiesslider_Callback(hObject,eventdata)
        fslider_value = get(hObject,'Value');
        slider_value = round(fslider_value);
        set(hObject,'Value',slider_value);
        set(ncitiessliderv,'String',slider_value);
        NVAR = round(slider_value);
    end
    function nindslider_Callback(hObject,eventdata)
        fslider_value = get(hObject,'Value');
        slider_value = round(fslider_value);
        set(hObject,'Value',slider_value);
        set(nindsliderv,'String',slider_value);
        NIND = round(slider_value);
    end
    function genslider_Callback(hObject,eventdata)
        fslider_value = get(hObject,'Value');
        slider_value = round(fslider_value);
        set(hObject,'Value',slider_value);
        set(gensliderv,'String',slider_value);
        MAXGEN = round(slider_value);
    end
    function mutslider_Callback(hObject,eventdata)
        fslider_value = get(hObject,'Value');
        slider_value = round(fslider_value);
        set(hObject,'Value',slider_value);
        set(mutsliderv,'String',slider_value);
        PR_MUT = round(slider_value)/100;
    end
    function crossslider_Callback(hObject,eventdata)
        fslider_value = get(hObject,'Value');
        slider_value = round(fslider_value);
        set(hObject,'Value',slider_value);
        set(crosssliderv,'String',slider_value);
        PR_CROSS = round(slider_value)/100;
    end
    function elitslider_Callback(hObject,eventdata)
        fslider_value = get(hObject,'Value');
        slider_value = round(fslider_value);
        set(hObject,'Value',slider_value);
        set(elitsliderv,'String',slider_value);
        ELITIST = round(slider_value)/100;
        GGAP = 1-ELITIST;
    end
    function crossover_Callback(hObject,eventdata)
        crossover_value = get(hObject,'Value');
        crossovers = get(hObject,'String');
        CROSSOVER = crossovers(crossover_value);
        CROSSOVER = CROSSOVER{1};
    end
    function mutation_Callback(hObject,eventdata)
        mutation_value = get(hObject,'Value');
        mutations = get(hObject,'String');
        MUTATION = mutations(mutation_value);
        MUTATION = MUTATION{1};
    end
    function selection_Callback(hObject,eventdata)
        methods = ["fitness_prop" "tournament3" "tournament5"];
        SELECTION = get(hObject,'Value') -1;
        if(SELECTION > 0)
            SELECTION = methods(SELECTION);
        else
            SELECTION = 'ranksus';
        end
    end
    function runbutton_Callback(hObject,eventdata)
        %set(ncitiesslider, 'Visible','off');
        set(nindslider,'Visible','off');
        set(genslider,'Visible','off');
        set(mutslider,'Visible','off');
        set(crossslider,'Visible','off');
        set(elitslider,'Visible','off');
        run_ga(x, y, NIND, MAXGEN, MAXTIME, NVAR, ELITIST, STOP_PERCENTAGE, PR_CROSS, PR_MUT, SELECTION, CROSSOVER, MUTATION, LOCALLOOP, EARLY, STALL_PERCENTAGE, RESET, ah1, ah2, ah3);
        end_run();
    end
    function inputbutton_Callback(hObject,eventdata)
        [x y] = input_cities(NVAR);
        axes(ah1);
        plot(x,y,'ko')
    end
    function end_run()
        %set(ncitiesslider,'Visible','on');
        set(nindslider,'Visible','on');
        set(genslider,'Visible','on');
        set(mutslider,'Visible','on');
        set(crossslider,'Visible','on');
        set(elitslider,'Visible','on');
    end
end

function [x, y] = input_cities(ncities)
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

