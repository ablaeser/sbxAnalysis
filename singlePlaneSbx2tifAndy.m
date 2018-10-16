function singlePlaneSbx2tifAndy(mouse, date, run, pmt, binSize, sbxType, maxFrame )   

    fprintf('\nUsing ablaeser version of singlePlaneSbx2tif\n')
    % default is just read green channel
    if nargin<4, pmt = 0; end
    if nargin<5, binSize = 1; end
    if nargin<6, sbxType = 'sbxreg'; end
    if nargin<7, maxFrame = Inf; end
    
    % Get file path and metadata
    inputPath = sbxPath(mouse, date, run, sbxType); 
    fprintf('Reading %s \n', inputPath)
    inf = sbxInfo(inputPath, true);
    edges = sbxRemoveEdges(inputPath);
    if ~isfield(inf, 'volscan') && length(inf.otwave)>1, error('This function is only used for single plate sbx frames.'); end

    % Set in to read the whole file if unset
    Nin = inf.max_idx + 1; % # of frames in the original movie
    Nout = floor(Nin/binSize); % # of frames to write out
    if Nout > maxFrame, Nout = maxFrame;  end
    tic;
    % Read the file and remove edges
    movie = fread(inf.fid, inf.nsamples/2*Nin, 'uint16=>uint16');
    movie = reshape(movie, [inf.nchan, inf.sz(2), inf.recordsPerBuffer, Nin]);
    movie = squeeze( movie(pmt+1,:,:,:) );
    movie = movie(edges(1)+1:end-edges(2), edges(3)+1:end-edges(4),:);
    movie = 65535 - permute(movie, [2,1,3] );  
    toc
    
    % Write the movie to tiff
    if pmt == 0, chn = 'greenChn'; else, chn = 'redChn'; end
    [inputDir,inputName] = fileparts( inputPath );
    outputPath = sprintf('%s\\%s_%s_%dFrameBin.tif', inputDir, inputName, chn, binSize  );
    if exist(outputPath, 'file') == 2, delete(outputPath); end
    fprintf('\nWriting %s ...\n', outputPath)
    SATopt = struct('overwrite',true, 'message',false, 'append',false, 'big',true );
    if binSize == 1
        outMovie = movie(:,:,1:Nout);
    else      
        outMovie = zeros( size(movie,1), size(movie,2), Nout );
        for z = 1:Nout %(Nout-1)  
            %saveastiff( uint16( mean(movie(:, :, binSize*(z-1)+1:binSize*z-1), 3) ), outputPath, SATopt);
            if rem(z, 500) == 0, fprintf('%d/%d   ', z, Nout-1); toc; end
            outMovie(:,:,z) = mean(movie(:, :, binSize*(z-1)+1:binSize*z-1), 3);
        end
    end
    saveastiff( uint16(outMovie), outputPath, SATopt);
    fprintf('Done\n')
    toc
end