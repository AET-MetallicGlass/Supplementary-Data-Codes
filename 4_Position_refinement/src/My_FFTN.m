function Y = My_FFTN(X)
Y = fftshift(fftn(ifftshift(X)));
end
