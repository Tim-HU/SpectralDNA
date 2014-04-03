function num = seq2num(seq)
    codes = -ones(1,255);
    nucl  = ['a','t','g','c'];
    codes(nucl) = 0:3;
    codes(upper(nucl)) = 0:3;    
    num = codes(seq);
end