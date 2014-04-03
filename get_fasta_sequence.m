function seq = get_fasta_sequence( filename, k )

    % Load data
    % g = importdata(filename);
    % importdata works in matlab2013 only. fastaread allows backward
    % compatibility with matlab2012.
    g = fastaread(filename); 
    % Apply selector
    if nargin == 1
        glen = length(g);
        seq  = cell(1,glen);
        for i = 1:glen
            seq{i} = g(i).Sequence;
        end
    else
        seq = g(k).Sequence;
    end
end