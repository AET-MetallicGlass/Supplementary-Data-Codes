%%  my_ifft %%

%% Helper function to calculate inverse fft

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function realout = my_ifft(k)
realout =fftshift((ifftn(ifftshift(k))));
end