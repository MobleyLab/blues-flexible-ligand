function avgProb = blockAvg( prob , numBlock)
    len = length( prob);
    lenBlock = len/numBlock;
    for block = 1:numBlock
        
        for numSubBlock = 2:10
            lenSubBlock = lenBlock / numSubBlock ;
            blockData = prob( (block -1)*lenBlock + 1 : block*lenBlock );
            clear avgSubBlock
            for subBlock = 1:numSubBlock
                start= ceil((subBlock -1)*lenSubBlock )+1;
                last = ceil(subBlock*lenSubBlock );
                avgSubBlock( subBlock  ) = mean(  blockData(start : last)  );
            end
            avgBlock = mean( avgSubBlock );
            avgSubBlockStd( numSubBlock -1) = std(avgSubBlock )/sqrt(numSubBlock);
        end
        avgProb( block, : ) = [avgBlock ,  max( avgSubBlockStd )];
    end
    
end