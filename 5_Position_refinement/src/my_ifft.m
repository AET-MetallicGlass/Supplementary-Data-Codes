%%  my_ifft %%

%% Helper function to calculate inverse fft

function realout = my_ifft(k)
realout =fftshift((ifftn(ifftshift(k))));
end