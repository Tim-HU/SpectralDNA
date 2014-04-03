function test_wht_noise()

    % Run `matlabpool open 4` beforehand for faster computation

    s = get_fasta_sequence('sequences.fasta',1);
    s_short = s(1:256)
    wht_noise_study( s_short, 25, 'insert' );

end
