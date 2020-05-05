%%  my_fft %%
% Helper function to calculate forward fft

function kout = my_fft(img)
kout = fftshift(fftn((ifftshift(img))));
end
