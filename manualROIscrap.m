%% Select data to work on
clear; clc;
%mouse = 'DL67';  date = '170601';  run = 5; 
%mouse = 'DL89';  date = '171119';  run = 6; 
%mouse = 'DL89';  date = '171119';  run = 4;
%mouse = 'DL72'; date = '170508'; run = 1;
%mouse = 'DL115'; date = '180704'; run = 1;
%mouse = 'DL72'; date = '170614'; run = 4;
%mouse = 'DL122'; date = '181012'; run = 1; % Z-stacks
mouse = 'DL128'; date = '181003'; run = 1; % single plane


%% Convert sbx to tif
pmt = 1; sbxType = 'sbx'; binSize = 10;  maxFrame = 3000; 
singlePlaneSbx2tifAndy( mouse, date, run, pmt, binSize, sbxType, maxFrame ) % animalID, dateID, run, pmt, binSize, sbxType, maxFrame
%p = struct('type','sbxreg', 'server',[], 'startframe',0, 'optolevel',[], 'estimate',false );
%sbxFirstLastAndy(mouse, date, run, 300, 0 ) %sbxFirstLast(mouse, date, run, 300, 0 )

%% 
%clear; clc;
%mouse = 'DL89';  date = '171119';  run = 4;
%mouse = 'DL67';  date = '170601';  run = 6; 
%mouse = 'DL115'; date = '180704'; run = 1;
% Get simpcell data
sc = sbxLoad(mouse, date, run, 'simpcell');
speed = sbxSpeed(mouse, date, run ); speed = speed';

%% Extract PMT0 (green channel) signals from sbx
frameRate = 15.49;
% Get file path and metadata
dataPath = sbxPath(mouse, date, run, 'sbxreg'); 
dataFolder = fileparts( dataPath );  dataFolder = [dataFolder,'\'];
inf = sbxInfo(dataPath, true);
edges = sbxRemoveEdges(dataPath);
Nframe = inf.max_idx + 1; 

% Read sbxreg, remove edges, rearrange and rescale
tic
fprintf('Reading %s\n', dataPath)
movie = fread(inf.fid, inf.nsamples/2*Nframe, 'uint16=>uint16');
movie = reshape(movie, [inf.nchan inf.sz(2) inf.recordsPerBuffer Nframe]);
movie = squeeze(movie);
movie = movie(edges(1)+1:end-edges(2), edges(3)+1:end-edges(4),:);
movie = 65535 - permute(movie, [2,1,3] );
movieDim = size( movie );
toc

% Save/load summary projections
tic
meanProjPath = sprintf('%s%s_%s_run%d_MEAN.tif', dataFolder, mouse, date, run);
stdProjPath = sprintf('%s%s_%s_run%d_STD.tif', dataFolder, mouse, date, run);
maxProjPath = sprintf('%s%s_%s_run%d_MAX.tif', dataFolder, mouse, date, run);

SATopt = struct('overwrite',true, 'message',false );
if exist(meanProjPath, 'file')
    movieMean = loadtiff(meanProjPath);
else
    movieMean = mean( movie, 3 ); saveastiff( uint16(movieMean), meanProjPath, SATopt ); toc
end

if exist(stdProjPath, 'file')
    movieStd = loadtiff(stdProjPath);
else
    movieStd = std( double(movie), 0, 3 ); saveastiff( uint16(movieStd), stdProjPath, SATopt ); toc
end

if exist(maxProjPath, 'file')
    movieMax = loadtiff(maxProjPath);
else
    movieMax = max(movie,[],3); saveastiff( uint16(movieMax), maxProjPath, SATopt ); toc
end

% Show summary projections
figure('units','normalized','OuterPosition',[0,0,1,1]);
sp(3) = subplot(1,3,3); imshow( movieMax, [] ); title('Max');
sp(2) = subplot(1,3,2); imshow( movieStd, [] ); title('Std Dev');
sp(1) = subplot(1,3,1); imshow( movieMean , [] );  title('Mean');
linkaxes(sp, 'xy')
impixelinfo;


%% Import the ROIs from ImageJ
roiPath = sprintf('%s%s_%s_run%d_RoiSet.zip', dataFolder, mouse, date, run);
fprintf('Importing ROIs: %s \Nframe', roiPath);
roi = ReadImageJROI(roiPath);
Nroi = numel( roi );
padDist = 3;
% Build masks from the imported ROIs
blankBin = false( movieDim(1), movieDim(2) );  allROI = blankBin;
roiProps = struct('area',NaN, 'box',nan(1,4), 'cent',nan(1,2), 'clicks',nan(0,2), 'edge',nan(0,2), 'image',[], 'ind',[], 'length',NaN, 'mask',[], 'orientation',NaN, 'pix',[]);
for r = 1:Nroi
    tempMask = blankBin; 
    tempLength = 0;
    if strcmp( roi{r}.strType, 'PolyLine' )
        roiProps(r).clicks = roi{r}.mnCoordinates + repmat( [1,1], size(roi{r}.mnCoordinates,1), 1 ); 
        roiPix = zeros(0,2);
        for s = 1:(size(roi{r}.mnCoordinates,1)-1)
            [segPix, segLength] = lineSegPix( flip( roiProps(r).clicks(s,:) ), flip( roiProps(r).clicks(s+1,:) ), 0.4 );
            %[segPix, segLength] = lineSegPix( flip( roi{r}.mnCoordinates(s,:)+[1,1] ), flip( roi{r}.mnCoordinates(s+1,:)+[1,1] ), 0.4 ); % ImageJ uses [0,0], not [1,1] as first pixel coord
            tempLength = tempLength + segLength;
            roiPix = [roiPix; segPix];
        end
        roiProps(r).length = tempLength; 
        roiInd = sub2ind( movieDim(1:2), roiPix(:,1), roiPix(:,2) );
        roiInd = unique( roiInd );
        tempMask( roiInd ) = true;
        %imshow( tempMask ); pause; clf;
        % Pad the ROI
        tempDist = bwdist(tempMask);
        edgeInd = find( tempDist < padDist );
        roiProps(r).mask = tempMask;  roiProps(r).mask(edgeInd) = true;
    else
        % vnRectBounds = rectangular bounds of the ROI: ['nTop', 'nLeft', 'nBottom', 'nRight']. ImageJ uses [0,0], not [1,1] as first pixel coord
        [x,y] = meshgrid( 1:movieDim(2), 1:movieDim(1) );
        ellipseCent = [mean(roi{r}.vnRectBounds([2,4])), mean(roi{r}.vnRectBounds([1,3])) ];
        ellipseRad = [diff(roi{r}.vnRectBounds([2,4])), diff(roi{r}.vnRectBounds([1,3]))]/2;
        ellipseInd = find( ((x - ellipseCent(1))/(ellipseRad(1))).^2 + ((y - ellipseCent(2))/(ellipseRad(2))).^2 <= 1 );
        tempMask( ellipseInd ) = true; 
        %figure; imshow( tempMask ); pause; clf;
        %yRange = roi{r}.vnRectBounds([1,3])+[1,1];  xRange = roi{r}.vnRectBounds([2,4])+[1,1]; % 
        %tempMask( yRange(1):yRange(2), xRange(1):xRange(2) ) = true;
        roiProps(r).mask = tempMask;
    end
    tempProps = regionprops( roiProps(r).mask, 'Area','BoundingBox','Centroid','Image','Orientation','PixelIdxList','PixelList' );
    roiProps(r).area = tempProps.Area; 
    roiProps(r).box = round(tempProps.BoundingBox); %
    roiProps(r).cent = round(tempProps.Centroid); %
    roiProps(r).image = tempProps.Image; 
    roiProps(r).orientation = tempProps.Orientation; 
    roiProps(r).pix = tempProps.PixelList; 
    roiProps(r).ind = tempProps.PixelIdxList; 
    tempEdge = bwboundaries( roiProps(r).mask );
    roiProps(r).edge = [tempEdge{1}(:,2), tempEdge{1}(:,1)];
    allROI( roiProps(r).ind ) = true;
end

figure; imshow( allROI, [] ); impixelinfo;

%% Exclude overlapping pixels 

%% Overlay ROIs on the std dev projection
figure;
imshow( movieStd, [] ); hold all;
for r = 1:Nroi-1
    plot( roiProps(r).edge(:,1), roiProps(r).edge(:,2), '--' ); hold all;
end
plot( roiProps(Nroi).edge(:,1), roiProps(Nroi).edge(:,2), 'w--' );

%% Extract signals 
tic
T = (1/frameRate)*(0:Nframe-1);
F = zeros( Nframe, Nroi );

for z = 1:Nframe %0:(nf-1) %(nf-1) 
    tempFrame = movie(:,:,z);
    for r = 1:Nroi
        F(z,r) = mean( tempFrame( roiProps(r).ind ) );
    end
    if rem(z, 500) == 0, fprintf('%d/%d   ', z, Nframe); toc; end
    %figure; imshow( tempFrame, [] ); hold on; %plot( roiProps(r).edge(:,1), roiProps(r).edge(:,2), 'w--' )
end
f = F(:,1:Nroi-1) - repmat( F(:,Nroi), 1, Nroi-1 );
toc
clearvars movie tempFrame % edges

%% Plot the raw fluor traces
figure;
plot( T, F(:,Nroi), 'k' ); hold on; 
xlabel('Time (s)'); ylabel('F (raw)');
for r = 1:Nroi-1
    if rem(r,3) == 0, cla; end
    plot( T, F(:,Nroi), 'k' ); hold all; 
    plot( T , F(:,r) ); %hold all; % 1:Nframe
    pause;
end

%% Rescale dF/F, filter and show all traces
windowSize = round(0.2*frameRate);
b = (1/windowSize)*ones(1,windowSize);  a = 1;
Flow = filtfilt( b, a, F );
fLow = Flow(:,1:Nroi-1) - repmat( Flow(:,Nroi), 1, Nroi-1 );
% {
close all;
figure('Units','normalized','OuterPosition',[0,0,1,1]);
for r = 1:Nroi
    plot( T, F(:,r), 'b' ); hold on;
    plot( T, Flow(:,r), 'k', 'LineWidth',2 );
    pause;
    cla;
end
%}

%% convert to dF/Fo
dFF = zeros( Nframe, Nroi ); % background subtracted deltaF/Fo
for r = 1:Nroi-1
    [~, ~, Otemp] = deconvolveCa( fLow(:,r), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); % use deconvolveCa to estimate Fo
    dFF(:,r) = (fLow(:,r) - Otemp.b)/Otemp.b;
end
[~, ~, Oback] = deconvolveCa( Flow(:,Nroi), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); % use deconvolveCa to estimate Fo
dFF(:,Nroi) = (Flow(:,Nroi) - Oback.b)/Oback.b;
dFFlim = [floor( min(min(dFF)) ), ceil( max(max(dFF)) )]; %yMin = floor( min(min(dFF)) ); yMax = ceil( max(max(dFF)) );  

close all;
figure('Units','normalized','OuterPosition',[0,0,1,1]);
sp(1) = subplot(2,1,1); sp(2) = subplot(2,1,2); linkaxes(sp,'x');
for r = Nroi
    % Plot the results
    subplot(sp(1)); cla;
    if r < Nroi
        plot( T, fLow(:,r) );
    else
        plot( T, Flow(:,r));
    end
    %title( sprintf('ROI %d, SN =  %2.1f, Fo = %2.1f', r, Otemp.sn, Otemp.b ) );
    ylabel('f_L_o_w');
    subplot(sp(2)); cla;
    plot( T, dFF(:,r) );
    xlim([0, T(end)]); 
    ylabel('dF/Fo'); xlabel('Time (s)');
    pause;
end

CalSpread( {T}, {dFF} ); % , 'FontSize',24

%% Load the wheel speed data

figure;
sp(1) = subplot(2,1,2);
plot( T, speed(1:Nframe) );
%plot( T, sc.brainmotion(1:Nframe) );
sp(2) = subplot(2,1,1);
linkaxes(sp,'x')
for r = 1:Nroi-1
    plot( T, dFF(:,r) ); hold all;
    pause;
end


%% Show big picture, zoomed in and fluor trace for each ROI
close all;
figure('Units','normalized','OuterPosition',[0,0,1,1]);
sp(1) = subplot(2,2,1);
sp(2) = subplot(2,2,2);
sp(3) = subplot(2,2,3:4);
backEdge = bwboundaries( roiProps(Nroi).mask );
subplot(sp(1))
imshow( movieMean, [] ); hold on;
plot( backEdge{1}(:,2), backEdge{1}(:,1), '--w' );   
for r = 1:Nroi-1
    subplot(sp(1))   
    plot( roiProps(r).edge(:,1), roiProps(r).edge(:,2), '--' ); hold all;
    subplot(sp(2)); % cla;
    zoomIm = movieMean( round(roiProps(r).box(2)):round(roiProps(r).box(2)+roiProps(r).box(4)), ...
        round(roiProps(r).box(1)):round(roiProps(r).box(1) + roiProps(r).box(3) ));
    imshow(  zoomIm, []  ); hold on;
    if strcmp( roi{r}.strType, 'PolyLine' )
        plot( roiProps(r).clicks(:,1) - roiProps(r).box(1), roiProps(r).clicks(:,2) - roiProps(r).box(2), 'rx'  )
    end
    hold off;
    
    subplot(sp(3))
    plot( T, dFF(:,r) );
    %{
    if strcmp( roi{r}.strType, 'PolyLine' )
        plot( T, f(:,r) ); hold all; %plot( T, F(:,r) ); hold all;
    else
        plot( T, F(:,Nroi), 'k' );
    end
    %}
    xlabel('Time (s)'); ylabel('dF/Fo'); title( sprintf('ROI %d', r ) );
    hold on;
    xlim([T(1), T(end)]); ylim([dFFlim(1), dFFlim(2)]);
    impixelinfo;
    pause;
    cla;
end

%% Calculate correlations between each pair
%{
corrF = corrcoef( F ); corrF( corrF == 1 ) = NaN;
corrf = corrcoef( f ); corrf( corrf == 1 ) = NaN;
figure; 
subplot(1,2,2); imagesc( corrf ); axis square; title('Background Subtracted');
subplot(1,2,1); imagesc( corrF ); axis square; title('Background Included (last channel)');
%}

corrDFF = corrcoef( dFF ); corrDFF( corrDFF == 1 ) = NaN;
figure; 
subplot(1,2,1);
imagesc( corrDFF ); axis square; title('Background Subtracted'); colorbar; impixelinfo;
set(gca,'TickDir','out')

% Distribution of correlation values
subplot(1,2,2);
hist( corrDFF(:), -1:0.02:1 ); xlabel('Correlation'); ylabel('# of ROI pairs');
set(gca,'TickDir','out')

%% Examine pairs of ROIs with high correlations
[x, y] = ind2sub( [Nroi, Nroi], find( triu( corrDFF ) > 0.5 ) ); disp([x,y])
figure;
sp(1) = subplot(3,1,2); sp(2) = subplot(3,1,3);
linkaxes(sp,'x');
for q = 1:length(x)
    subplot(3,1,1); cla
    imshow( movieMean, [] ); hold on;
    plot( roiProps(x(q)).edge(:,1), roiProps(x(q)).edge(:,2), 'b--' ); hold all;
    plot( roiProps(y(q)).edge(:,1), roiProps(y(q)).edge(:,2), 'r--' ); hold all;
    subplot(sp(1)); cla
    plot( T, dFF(:,x(q)), 'b' );  ylabel('dF/Fo');
    title(sprintf('ROI %d and %d, correlation = %2.3f', x(q), y(q), corrDFF(x(q),y(q))) );
    hold on;
    subplot(sp(2)); cla
    plot( T, dFF(:,y(q)), 'r' ); hold on;
    xlabel('Time (s)'); ylabel('dF/Fo');
    xlim([T(1), T(end)]);
    pause; 
end


%% Filter and deconvolve the fluorescence signals
C = zeros(Nframe, Nroi); S = zeros(Nframe, Nroi);
for r = Nroi:-1:1 
    [Ctemp, Stemp, Otemp] = deconvolveCa( dFF(:,r), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); % , 'optimize_smin'
    C(:,r) = Ctemp; S(:,r) = Stemp; O(r) = Otemp;
end
fprintf('\nDeconvolution complete\n');

%% Deconvolve the background trace and show results
Sback = S(:,Nroi); Cback = C(:,Nroi); Oback = O(Nroi);
S99 = prctile(Sback(Sback>0), 99 );
figure;
hist( Sback(Sback>0), 0:10:1000 )
title( sprintf('99 percentile = %2.2f', S99) ); % sprintf('MAD = %2.2f', mad(Sback(Sback>0)) )

close all;
figure('Units','normalized','OuterPosition',[0,0,1,1]);
sp(1) = subplot(3,1,1); sp(2) = subplot(3,1,2); sp(3) = subplot(3,1,3);
linkaxes(sp,'x');
subplot(sp(1)) 
plot( T, dFF(:,Nroi) );
ylabel('dF/Fo'); title( sprintf('ROI %d  SN =  %2.1f, b = %2.1f', Nroi, Oback.sn, Oback.b ) );
set(gca,'TickLength',[0,0]);
subplot(sp(2));
plot( T, Cback );
ylabel('[Ca]');
xlim([T(1), T(end)]);
set(gca,'TickLength',[0,0]);
subplot(sp(3));
plot( T, Sback ); hold on; 
line([T(1), T(end)], S99*[1,1], 'Color', 'k');
xlabel('Time (s)'); ylabel('Pseudospike rate');
xlim([T(1), T(end)]);
set(gca,'TickLength',[0,0]);
title( sprintf('par = %2.2f, S_9_9 = %2.3f', Oback.pars, S99 ) );

%% Show the deconvolution results

Smax = ceil( max(max(S)) );

close all;
figure('Units','normalized','OuterPosition',[0,0,1,1]);
sp(1) = subplot(3,1,1); sp(2) = subplot(3,1,2); sp(3) = subplot(3,1,3);
linkaxes(sp,'x');
for r = 1:Nroi
    subplot(sp(1)); cla;
    plot( T, dFF(:,r) ); hold on;
    plot( T, dFF(:,Nroi), 'k' );
    
    ylabel('dF/Fo'); title( sprintf('ROI %d', r ) );
    set(gca,'TickLength',[0,0]); ylim([yMin, yMax]);
    
    subplot(sp(2)); cla;
    plot( T, C(:,r) ); hold on;
    plot( T, C(:,Nroi), 'k' );
    ylabel('[Ca]');
    xlim([T(1), T(end)]);
    set(gca,'TickLength',[0,0]);  ylim([0, yMax]);
    title( sprintf('SN =  %2.3f', O(r).sn ) );
    
    subplot(sp(3)); cla;
    plot( T, S(:,r) ); hold on;
    plot( T, S(:,Nroi), 'k' );
    line([T(1), T(end)], S99*[1,1], 'Color', 'k');
    xlabel('Time (s)'); ylabel('Pseudospike rate');
    xlim([T(1), T(end)]); ylim([0,Smax]);
    set(gca,'TickLength',[0,0]);
    title( sprintf('par = %2.2f', O(r).pars ) );
    pause;
end

%% Parameters by ROI
figure;
for r = 1:Nroi
    %plot( O(r).b, O(r).sn, '.' ); hold all;
    %xlabel('Background'); ylabel('Noise');
    plot( O(r).pars, O(r).sn, '.' ); hold all;
    xlabel('AR1 parameter'); ylabel('Noise');
end
axis square;

%% Correlation using S instead
corrS = corrcoef( S ); corrS( corrS == 1 ) = NaN;
figure; 
subplot(1,2,1);
imagesc( corrS ); axis square; title('Background Subtracted'); colorbar; impixelinfo;
set(gca,'TickDir','out')

% Distribution of correlation values
subplot(1,2,2);
hist( corrS(:), -1:0.02:1 ); xlabel('Correlation'); ylabel('# of ROI pairs');
set(gca,'TickDir','out')

%% Examine pairs of ROIs with high correlation in S
[x, y] = ind2sub( [Nroi, Nroi], find( triu( corrS ) > 0.5 ) ); disp([x,y])
figure;
sp(1) = subplot(3,1,2); sp(2) = subplot(3,1,3);
linkaxes(sp,'x');
for q = 1:length(x)
    subplot(3,1,1); cla
    imshow( movieMean, [] ); hold on;
    plot( roiProps(x(q)).edge(:,1), roiProps(x(q)).edge(:,2), 'b--' ); hold all;
    plot( roiProps(y(q)).edge(:,1), roiProps(y(q)).edge(:,2), 'r--' ); hold all;
    subplot(sp(1)); cla
    plot( T, S(:,x(q)), 'b' );  ylabel('S');
    title(sprintf('ROI %d and %d, correlation = %2.3f', x(q), y(q), corrS(x(q),y(q))) );
    hold on;
    subplot(sp(2)); cla
    plot( T, S(:,y(q)), 'r' ); hold on;
    xlabel('Time (s)'); ylabel('S');
    xlim([T(1), T(end)]);
    pause; 
end
%CalSpread( {T}, {S} );

%% Correlation between activity and running
corrSpeed = zeros(1,Nroi); pSpeed = zeros(1,Nroi);
for r = 1:Nroi
    [corrSpeed(r), pSpeed(r)] = corr( S(:,r), speed );
end

hist( corrSpeed(:), -1:0.02:1 ); xlabel('Correlation'); ylabel('# of ROIs');
set(gca,'TickDir','out')

plot( corrSpeed, pSpeed, '.' ); 
xlabel('Rho'); ylabel('pValue');
axis square;

figure;
sp(1) = subplot(2,1,1); sp(2) = subplot(2,1,2);
linkaxes(sp,'x');
subplot(sp(2)); cla
plot( T, speed, 'r' ); hold on;
xlabel('Time (s)'); ylabel('Speed');
for r = find( pSpeed < 0.05 & corrSpeed > 0 ) 
    subplot(sp(1)); cla
    plot( T, S(:,r), 'b' );  ylabel('S');
    title(sprintf('ROI %d, correlation = %2.3f, pVal = %2.3f', r, corrSpeed(r), pSpeed(r) ) );
    xlim([T(1), T(end)]);
    pause; 
end

%% Test performance of different deconvolutions
close all;
figure('Units','normalized','OuterPosition',[0,0,1,1]);
sp(1) = subplot(3,1,1); sp(2) = subplot(3,1,2); sp(3) = subplot(3,1,3);
linkaxes(sp,'x');
for r = 1:Nroi
    [Car1, Sar1, Oar1] = deconvolveCa( dFF(:,r), 'ar1', 'constrained', 'optimize_b', 'optimize_pars', 'smin',-3 ); % appears to yield best results
    %[Car2, Sar2, Oar2] = deconvolveCa( f(:,r), 'ar1', 'thresholded', 'optimize_b',  'optimize_pars', 'smin',-3 );
    %[Cexp, Sexp, Oexp] = deconvolveCa( f(:,r), 'exp2',  'foopsi', 'optimize_pars' );
    subplot(sp(1)) 
    plot( T, f(:,r) );
    ylabel('f'); title( sprintf('ROI %d', r ) );
    set(gca,'TickLength',[0,0]);
    
    subplot(sp(2));
    plot( T, Car1, 'b' ); hold on;
    plot( T, Car2, 'r' );
    %plot( T, Cexp, 'k--' );
    ylabel('[Ca]');
    xlim([T(1), T(end)]);
    set(gca,'TickLength',[0,0]);
    
    subplot(sp(3));
    plot( T, Sar1, 'b' ); hold on;
    plot( T, Sar2, 'r' );
    %plot( T, Sexp, 'k--' );
    xlabel('Time (s)'); ylabel('Pseudospike rate');
    xlim([T(1), T(end)]);
    set(gca,'TickLength',[0,0]);
    %title( sprintf('SN =  %2.1f, b = %2.1f', Otemp.sn, Otemp.b ) );
    pause;
end

%% 
plot( T, sc.brainmotion );
plot( T, sc.shearml );
plot( T, sc.shearap );
plot( T, sc.transml );
plot( T, sc.transap );
plot( T, sc.scaleml );
plot( T, sc.scaleap );