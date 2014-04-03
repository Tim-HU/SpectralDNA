function savename = wht_noise_study(seq,rep)

    % Transform sequence
    seq = seq2num(seq);
    
    % Parameters
    seqlen = length(seq);
    order  = nextpow2(seqlen);
    whtlen = pow2(order);
    
    % Number of repeats
    if nargin == 1
        rep = ceil(seqlen / 2);
    end
    
    % Precomputations
    seq_wht = dna_wht(seq);
    res     = zeros(4,whtlen,seqlen);

    % For an increasing number of mutations
    for n_mut = 1:seqlen
        
        dif = zeros(4,whtlen,rep);
        
        % Repeat the experiement rep times
        for s = 1:rep
            
            % Create mutated sequence
            loc = randi(seqlen,1,n_mut);
            snp = randi(3,1,n_mut);
            mut = seq;
            mut(loc) = mod( seq(loc) + snp, 4 );
            
            % Compute absolute difference of WH coefficients
            dif(:,:,s) = abs( seq_wht - dna_wht(mut) );
        end
        
        % Summarize results
        res(:,:,n_mut) = sqrt(mean( dif, 3 ));        
        fprintf( '%d/%d mutations done...\n', n_mut, seqlen );
        
    end
    disp('Done!');
    
    % Save results
    savename = sprintf('wht_noise_study_%s.mat',datestr(now,'mmmdd-HH:MM:SS'));
    save(savename,'seq','rep','seq_wht','res');
    
    % Show results
    wht_noise_show(savename);    
    
    %-------------------------------------------------------
    %-------------------------------------------------------
    
    function trf_ = dna_wht(seq_)        
        trf_ = zeros(4,whtlen);
        for n = 1:4
            trf_(n,:) = fwht( binarize(seq_,n-1), [], 'sequency' );
        end        
    end

    function bin_ = binarize(seq_,val_)
        bin_ = seq_ == val_;
        bin_ = 2*bin_ - 1;
    end
    
end