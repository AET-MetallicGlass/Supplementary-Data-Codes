% My_paddzero
% Y. Yang, UCLA Physics & Astronomy
% First version date: 2015. 07. 13.

% output parameter: PadVol: padded N-d array
% input parameters: Vol (original N-d array before padding)
%                   paddedsize (size of N-d array after padding)

% My_paddzero padds zero to the input array, so that
% size(PadVol) = paddedsize.

% My_paddzero preserves the center pixel position
% at (n+1)/2 [in case of n=odd] or at n/2+1 [in case of n=even]
% in each dimension.

% length(size(Vol)) must be equal to length(paddedsize).
% paddedsize should be equal to or larger than size of Vol in all
% dimensions.

% Minh Rose, UCLA Math 2018.10.12
% Make it single matrix (from double) to save memory

% Yao Yang, 2020.07.13
% Coherent Imaging Group, UCLA 
% Make the data type an input (default single)

function PadVol = My_paddzero(Vol,paddedsize,dtype)
if nargin < 3
    dtype = 'single';
end
if length(size(Vol))~=length(paddedsize)
    fprintf(1,'input volume dimension and and paddedsize length does not match!\n');
elseif sum(size(Vol)>paddedsize) > 0
    fprintf(1,'paddedsize should be equal to or larger than original volume in all dimensions!\n');
else
%     %         PadVol = zeros(paddedsize,'single');
    PadVol = zeros(paddedsize, dtype);
    currevalstr = 'PadVol(';
    for i=1:length(size(Vol))
        if mod(size(Vol,i),2) == 0
            if mod(paddedsize(i),2) == 0
                startind = 1 + (paddedsize(i)-size(Vol,i))/2;
                endind = size(Vol,i) + (paddedsize(i)-size(Vol,i))/2;
            else
                startind = 1 + (paddedsize(i)-size(Vol,i)-1)/2;
                endind = size(Vol,i) + (paddedsize(i)-size(Vol,i)-1)/2;
            end
        else
            if mod(paddedsize(i),2) == 0
                startind = 1 + (paddedsize(i)-size(Vol,i)+1)/2;
                endind = size(Vol,i) + (paddedsize(i)-size(Vol,i)+1)/2;
            else
                startind = 1 + (paddedsize(i)-size(Vol,i))/2;
                endind = size(Vol,i) + (paddedsize(i)-size(Vol,i))/2;
            end
        end
        currevalstr = sprintf('%s%d:%d',currevalstr,startind,endind);
        if i < length(size(Vol))
            currevalstr = strcat(currevalstr,',');
        end
    end
    currevalstr = sprintf('%s) = Vol;',currevalstr);
    eval(currevalstr);
    
end
end