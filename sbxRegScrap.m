Nframe = size(movie,3);
Nref = 50;
refIm = mean( movie(:,:, round( Nframe/2 )+[-Nref:Nref] ), 3 ); imshow( refIm, [] );
refFFT = fft2( refIm );

regMovie = zeros( size(movie,1), size(movie,2), size(movie,3) );
tic
for z = 1:Nframe
    %imshow( movie(:,:,z), [] );
    [regInfo, regFrameFFT] = dftregistration( refFFT, fft2( movie(:,:,z) ) );
    regMovie(:,:,z) = ifft2( regFrameFFT );
    if rem(z,500)==0, fprintf('z = %d  ',z); toc, end
    %imshowpair( movie(:,:,z), uint16(regMovie(:,:,z)), 'montage' );
end

regPath = 'C:\2pdata\DL128\181003_DL128\181003_DL128_run1\redReg.tif';
fprintf('Writing %s ...', regPath )
SATopt = struct('overwrite',true, 'message',false, 'append',false, 'big',true );
saveastiff( uint16(regMovie), regPath, SATopt ); fprintf('\nDone'); toc