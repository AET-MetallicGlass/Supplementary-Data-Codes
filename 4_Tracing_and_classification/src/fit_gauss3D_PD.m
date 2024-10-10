function [x, resnorm,residual] = fit_gauss3D_PD(init_guess, xdata, ydata, varargin)

%fit_gauss3D_PD(init_guess, xdata, ydata, fixed, lowerbound,upperboud)
% Y. Yang, UCLA, 2015
% This function is just after slight modification from
% original version by M. Bartels, UCLA, 2014
% To make the Gaussian always positively defined

%Fit 3D gaussian and allow fixed parameters as well as bounds

%set standard values for optional parameters
fixed=zeros(size(init_guess));
lowerbound= zeros(size(init_guess)) - Inf;
upperbound= zeros(size(init_guess)) + Inf;

%check optional parameters
if nargin>3
   fixed=cell2mat(varargin(1));
end

if nargin>4
   lowerbound=cell2mat(varargin(2));
end

if nargin>5
   upperbound=cell2mat(varargin(3));
end

%convert to double
init_guess    = double(init_guess);
ydata = double(ydata);

%store all initial parameters
init_guess_all = init_guess;

%only fit the variable parameters
init_guess = init_guess_all(fixed == 0);
lowerbound = lowerbound(fixed == 0);
upperbound = upperbound(fixed == 0);


opt=optimset('TolFun',1e-12);
opt=optimset(opt,'Display','off');

%do the fit
[x, resnorm,residual] = lsqcurvefit(@Fhelp, init_guess, xdata, ydata, lowerbound, upperbound, opt);

%add fixed parameters to the final result again
init_guess_all(fixed == 0) = x;
x = init_guess_all;

    
    function y = Fhelp(x, xdata)
        %helper function to deal with fixed parameters
        
        %merge with fixed parameters
        x_current = x;
        x = init_guess_all;
        x(fixed == 0) = x_current;
        
        %calculate the acutal function
        y=calc_gauss3D_PD(x,xdata);

    end

end



