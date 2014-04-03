function wht_noise_show( filename, fignum )

    % Default input
    if nargin == 1, fignum = 1; end

    % Load data
    D = load(filename);
    
    % Parameters
    seqlen = length(D.seq);
    whtlen = pow2(nextpow2(seqlen));
    
    % Information
    fprintf('Showing results in "%s"...\n', filename);
    fprintf('Sequence of length %d with %d repeats.\n', seqlen, D.rep );
        
    % Show results
    figure(fignum)
    set(gcf,'renderer','opengl');
    
        [x,y] = ndgrid(1:whtlen,1:seqlen);
        nucl  = ['A','T','G','C'];

        for i = 1:4 

            subplot(2,2,i)
            h=surf( x,y,squeeze(D.res(i,:,:)) );
            set(h,'edgecolor','none');
            
            xlabel('WH coefficient');
            ylabel('# SNV');
            zlabel('RMSE');
            
            title(nucl(i));
            colormap hsv;
            shading interp;
            axis tight;
            
        end
    
end