% My_stripzero
% Y. Yang, UCLA Physics & Astronomy
% First version date: 2015. 07. 13.

% output parameter: PadVol: padded N-d array
% input parameters: PadVol (padded N-d array)   
%                   orisize (size of N-d array after removing padded zeros)                  

% My_stripzero removes zero from the input array, so that 
% size(Vol) = orisize.

% My_stripzero preserves the center pixel position 
% at (n+1)/2 [in case of n=odd] or at n/2+1 [in case of n=even] 
% in each dimension.

% length(size(PadVol)) must be equal to length(orisize).
% orisize should be equal to or smaller than size of PadVol in all
% dimensions.



function Vol = My_stripzero(PadVol,orisize)
    if length(size(PadVol))~=length(orisize)
        fprintf(1,'input volume dimension and and paddedsize length does not match!\n');
    elseif sum(size(PadVol)<orisize) > 0
        fprintf(1,'paddedsize should be equal to or smaller than original volume in all dimensions!\n');
    else      
        currevalstr = 'PadVol(';
        for i=1:length(size(PadVol))
            if mod(size(PadVol,i),2) == 0
                if mod(orisize(i),2) == 0
                    startind = 1 + (size(PadVol,i)-orisize(i))/2;
                    endind = orisize(i) + (size(PadVol,i)-orisize(i))/2;
                else
                    startind = 1 + (size(PadVol,i)-orisize(i)+1)/2;
                    endind = orisize(i) + (size(PadVol,i)-orisize(i)+1)/2;
                end
            else
                if mod(orisize(i),2) == 0
                    startind = 1 + (size(PadVol,i)-orisize(i)-1)/2;
                    endind = orisize(i) + (size(PadVol,i)-orisize(i)-1)/2;
                else
                    startind = 1 + (size(PadVol,i)-orisize(i))/2;
                    endind = orisize(i) + (size(PadVol,i)-orisize(i))/2;
                end
            end
            currevalstr = sprintf('%s%d:%d',currevalstr,startind,endind);
            if i < length(size(PadVol))
                currevalstr = strcat(currevalstr,',');
            end
        end
        currevalstr = sprintf('%s)',currevalstr);
        Vol = eval(currevalstr);        
    end    
end