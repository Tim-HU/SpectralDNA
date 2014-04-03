function test_wht_noise()

    % Run `matlabpool open 4` beforehand for faster computation

    s = get_fasta_sequence('sequences.fasta',1);
    wht_noise_study( s(1001:1128), 25 );

end