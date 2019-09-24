function probLeft = getProb(data)
   size(data)
   idxL = 0; idxR = 0;
   for i = 1:length(data)
        if data(i,2) < -11
            left( idxL +1 ) = data(i,2); idxL = idxL+1;
        else
            right( idxR +1 ) = data(i,2); idxR = idxR+1;
        end
        probLeft(i) = idxL / ( idxL +idxR);
   end

end
