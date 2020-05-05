function v = fatom_vector(q,Z)

% fparams[79][0] =  9.61263398e-001 ; a1
% fparams[79][1] =  1.70932277e-001 ; b1
% fparams[79][2] =  3.69581030e+000 ; a2
% fparams[79][3] =  1.29335319e+001 ; b2
% fparams[79][4] =  2.77567491e+000 ; a3
% fparams[79][5] =  6.89997070e-001 ; b3
% fparams[79][6] =  2.95414176e-001 ; c1
% fparams[79][7] =  1.63525510e-001 ; d1
% fparams[79][8] =  3.11475743e-001 ; c2
% fparams[79][9] =  1.39200901e+000 ; d2
% fparams[79][10] =  1.43237267e-002 ; c3
% fparams[79][11] =  2.71265337e-002 ; d3

% a = [9.61263398e-001 3.69581030e+000 2.77567491e+000];
% b = [1.70932277e-001 1.29335319e+001 6.89997070e-001];
% c = [2.95414176e-001 3.11475743e-001 1.43237267e-002];
% d = [1.63525510e-001 1.39200901e+000 2.71265337e-002];

fpara = fparameters(Z);
a = [fpara(1) fpara(3) fpara(5)];
b = [fpara(2) fpara(4) fpara(6)];
c = [fpara(7) fpara(9) fpara(11)];
d = [fpara(8) fpara(10) fpara(12)];

num = numel(q);
v = zeros(1,num);

for hh = 1:num
% Lorenzians %
suml = sum( a./((q(hh)^2)+b) );

% Gaussians %
sumg = sum( c.*exp(-(q(hh)^2).*d) );

v(hh) = suml + sumg;
end

% Lorenzians %
% suml = sum( a./((q^2)+b) );
% Gaussians %
% sumg = sum( c.*exp(-(q^2).*d) );
% v = suml + sumg;

end