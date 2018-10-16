function [Fspread, Flim, showSpread] = CalSpread( t, F, varargin ) % sep,
%CalSpread Shows all fluorescence traces spaced out, for each individual experiment
Nexp = size( F, 1 ); Nroi = size( F{1},2 );
% Parse inputs
IP = inputParser;
addRequired( IP, 't', @iscell )
addRequired( IP, 'F', @iscell )
addOptional( IP, 'Event', {}, @iscell )
addOptional( IP, 'tshade', cell(Nexp,1), @iscell )
addParameter( IP, 'iosc', [], @isnumeric )
addParameter( IP, 'seti', 1:Nroi, @isnumeric )
addParameter( IP, 'setk', 1:Nexp, @isnumeric )
addParameter( IP, 'save', '', @ischar )
addParameter( IP, 'FixSep', 0, @isnumeric )
addParameter( IP, 'title', '', @ischar )
addParameter( IP, 'unit', 'Fluorescence', @ischar )
addParameter( IP, 'FontSize', 18, @isnumeric )
addParameter( IP, 'show', true, @islogical )
parse( IP, t, F, varargin{:} );
Event = IP.Results.Event; Nevent = cellfun(@numel, Event);
tshade = IP.Results.tshade;
iosc = IP.Results.iosc;
seti = IP.Results.seti;
setk = IP.Results.setk;
savefile = IP.Results.save;
FixSep = IP.Results.FixSep;
FigTitle = IP.Results.title;
FigUnit = IP.Results.unit;
FontSize = IP.Results.FontSize;
show = IP.Results.show;

% Estimate separation (optional)
if FixSep == 0
    range = [];
    for k = 1:Nexp
        range = [ max( F{k}(:,:) ) - min( F{k}(:,:) ), range ]; %#ok<*AGROW>
    end
    sep = prctile( range, 90 ); 
else
    sep = FixSep;
end
% Find time limits
t = cellfun( @single, t, 'UniformOutput', false ); 
tmin = min( cellfun(@min, t) );  tmax = max( cellfun(@max, t) ); 
tlim = [tmin, tmax];  
% Generate Fspread, Flim
Fspread = cell(Nexp,1); target = 0:sep:sep*(Nroi-1);  Flim = [0,0]; % Fmax = 0; 
for k = setk
    Fspread{k} = single( F{k}(:,:) );
    Fspread{k} = Fspread{k} + repmat(target-mean(Fspread{k}),size(F{k},1),1); % adjust each trace so that its mean = target
    Flim = [min(min(Fspread{k}(:),Flim(1))), max(max(Fspread{k}(:)), Flim(2))]; 
end
Flim = [floor(Flim(1))-2, ceil(Flim(2))+5];
% Make the plot
if show
    showSpread = figure('units','normalized','outerposition',[0 0 1 1]);
    shade = [0.7,0.7,0.7]; %[0.5,0.8,0.3];
    for k = setk
        if ~isempty( tshade{k} )
            for j = 1:size(tshade{k},1)
                area( tshade{k}(j,:), [Flim(2), Flim(2)], 'facecolor', shade, 'edgecolor', shade, 'basevalue', Flim(1) ); hold on;
            end
        end
        if isempty( iosc )
            plot( t{k}, Fspread{k}(:,seti) ); hold on; 
        else
            plot( t{k}, Fspread{k}(:,seti), 'b' ); hold on;
            for i = intersect( iosc, seti ) 
                plot( t{k}, Fspread{k}(:,i), 'r');
            end
        end
        if ~isempty( Event )
            for i = seti
                for j = 1:Nevent(k,i) 
                    plot( t{k}( [Event{k,i}(j).ind] ), Fspread{k}( [Event{k,i}(j).ind], i ), 'k' );
                    for z = 1:Event{k,i}(j).Nsub
                        plot( t{k}( [Event{k,i}(j).sub(z).rise] ), Fspread{k}( [Event{k,i}(j).sub(z).rise], i ), 'b.' );
                    end
                end
            end
        end
    end
    xlim(tlim); ylim(Flim); % axis tight;
    xlabel('Time (s)', 'FontSize',FontSize); ylabel( FigUnit, 'FontSize',FontSize ); title( FigTitle, 'FontSize',FontSize );
    set(gca,'TickLength',[0.005,0], 'TickDir','out', 'YtickLabel',[], 'box','off', 'FontSize',FontSize);
else 
    showSpread = []; 
end
% Save the figure (optional)
if show && ~isempty( savefile )
    fprintf('\nSaving %s.\n', [savefile,'.fig'] );
    saveas( showSpread, savefile, 'fig' );
end
end