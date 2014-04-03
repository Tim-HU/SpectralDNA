function savename = wht_noise_study(seq,rep,noise_mode)

    % Transform sequence
    seq = seq2num(seq);
    
    % Parameters
    seqlen = length(seq);
    order  = nextpow2(seqlen);
    whtlen = pow2(order);
    
    % Default value for number of repeats.
    % Default value for noise_mode.
    if nargin == 1
        rep = ceil(seqlen / 2);
        noise_mode = 'mutate';
    elseif nargin == 2
        noise_mode = 'mutate';
    end
    
    % Precomputations
    seq_wht = dna_wht(seq);
    res     = zeros(4,whtlen,seqlen);

    if strcmp(noise_mode,'mutate')

        % For an increasing number of mutations
        for n_mut = 1:seqlen

            dif = zeros(4,whtlen,rep);

            % Repeat the experiement rep times
            for s = 1:rep

                % Create mutated sequence
                % loc = randi(seqlen,1,n_mut); 
                % randi changed to randsample because we want sample
                % without replacement.
                loc = randsample(seqlen,n_mut);
                snp = randi(3,1,n_mut);
                mut = seq;
                mut(loc) = mod( seq(loc) + snp, 4 );

                % Compute absolute difference of WH coefficients
                dif(:,:,s) = abs( seq_wht - dna_wht(mut) );
            end

            % Summarize results
            res(:,:,n_mut) = (mean( dif, 3 ));        
            fprintf( '%d/%d mutations done...\n', n_mut, seqlen );

        end

    elseif strcmp(noise_mode,'insert')
        % Insert increasing number of random individual nucleotides in
        % sequence, and due to sequence length is fixed, the last 1 bp is
        % removed at each 1bp insert event.
        
        % For an increasing number of inserts
        for n_mut = 1:seqlen

            dif = zeros(4,whtlen,rep);

            % Repeat the experiement rep times
            for s = 1:rep

                % Create inserted sequence
                loc = randsample(seqlen,n_mut);
                mut = seq;
                for each_loc = loc
                    mut((each_loc+1):seqlen) = mut(each_loc:(seqlen-1));
                    mut(each_loc) = randi(4);
                    %disp(each_loc);
                    %disp(mut);
                end

                % Compute absolute difference of WH coefficients
                dif(:,:,s) = abs( seq_wht - dna_wht(mut) );
            end

            % Summarize results
            res(:,:,n_mut) = (mean( dif, 3 ));        
            fprintf( '%d/%d insertions done...\n', n_mut, seqlen );

        end
        
    elseif strcmp(noise_mode,'shift')
        disp('shift not implemented yet.');
        
    else
        disp('Error: noise_mode argumant can only be eith mutate, insert or shift.');
    end
    disp('Done!');
    
    
    % Save results
    savename = sprintf('wht_noise_study_%s.mat',datestr(now,'mmmdd-HHMMSS'));
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