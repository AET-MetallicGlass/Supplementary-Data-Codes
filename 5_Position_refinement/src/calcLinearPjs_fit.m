function Projs = calcLinearPjs_fit(para,xdata)

angles  = xdata.angles;
model   = xdata.model;
atoms	= xdata.atoms;
Z_arr	= xdata.Z_arr;
Res     = xdata.Res;
volSize = xdata.volSize;

CropHalfWidth = xdata.CropHalfWidth;

htAr = para(1,:);
bfAr = para(2,:);

fixedfa_big = make_fixedfa_man(volSize, Res, Z_arr);

Projs = My_create_volProjs_from_model_exact_HB_fixedfa_fainput(...
    model,atoms,[htAr 0],[bfAr 1],volSize, Res, CropHalfWidth, angles',fixedfa_big);

end