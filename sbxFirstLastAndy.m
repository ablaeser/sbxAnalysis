function sbxFirstLastAndy(mouse, date, run, N, pmt, varargin)
%SBXFIRSTLAST Make TIFFs of the first and last 500 frames
    fprintf('\nUsing ablaeser version of sbxFirstLast\n')
    p = inputParser;
    addOptional(p, 'type', 'sbx');  %sbx file type to open: sbx, sbxreg, xyreg, sbxclean, demonsreg
    addOptional(p, 'server', []);  % Server name
    addOptional(p, 'startframe', 0, @isnumeric);
    addOptional(p, 'optolevel', []);
    addOptional(p, 'estimate', false);  % Whether to give an estimate of where the path would be if it does not exist
    parse(p, varargin{:});
    p = p.Results;
    
    if nargin < 4, N = 500; end
    if nargin < 5, pmt = 0; end
    
    path = sbxPath(mouse, date, run, p.type, 'pmt', pmt, 'server', p.server,'estimate',p.estimate);
    fprintf('\nRead %s\n', path );
    tic
    outputPath = sprintf('%s_first-last-%i.tif', path(1:strfind(path,'.')-1), N);
    if ~exist(outputPath, 'file')
        info = sbxInfo(path);
        firstN = sbxReadPMT(path, p.startframe, N, pmt, p.optolevel); % firstN = permute(firstN, [2,1,3]);
        lastN = sbxReadPMT(path, info.max_idx+1-N, N, pmt, p.optolevel);
        
        
        SATopt = struct('overwrite',true, 'message',false, 'append',false, 'big',true );
        %saveastiff( firstN,  outputPath, SATopt ); 
        %writetiff( firstN, outputPath, class(firstN));
        %{
        SATopt = struct('overwrite',false, 'message',false, 'append',true );
        for z = 1:N %(nf-1)  
            saveastiff( firstN(:,:,z), outputPath, SATopt);
            if rem(z, 100) == 0, fprintf('%d/%d   ', z, N); toc; end
        end
        %}
        
        saveastiff( cat(3, firstN, lastN), outputPath, SATopt );
        %writetiff(cat(3, firstN, lastN), outputPath, class(firstN));
        fprintf('\nWrote %s\n', outputPath );
        toc
    end
end
%{
figure;
for z = 1:size(firstN,3)
    imshow( lastN(:,:,z), [] ); %imshow( firstN(:,:,z), [] );
    %if z == 1,  impixelinfo;  end
    title( sprintf('z = %d', z ) );
    pause(0.1);
    cla;
end
%}
