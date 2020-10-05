function [Base, Corrected_Spectrum]=baseline_normal(Spectrum)
%Input
%-------
%Spectrum: vector of size (N*1)
%Output
%-------
%Base: Identified Baseline vector of size (N*1)
%Corrected_Spectrum: Corrected Spectrum vector of size (N*1)
l=length(Spectrum);  
lp=ceil(0.05*l);
initial_Spectrum=[ones(lp,1)*Spectrum(1) ; Spectrum ; ones(lp,1)*Spectrum(l)];
l2=length(initial_Spectrum);
S=initial_Spectrum;
n=1;
flag1=0;
while flag1==0
    
n=n+2;
    
i=(n-1)/2;
    
[Baseline, stripping]=peak_stripping(S,n);
A(i)=trapz(S-Baseline);
Stripped_Spectrum{i}=Baseline;
S=Baseline;
if i>3
    if A(i-1)<A(i-2) && A(i-1)<A(i)
        i_min=i-1;
        flag1=1;
    end 
end
end
Base=Stripped_Spectrum{i_min}; 
Corrected_Spectrum=initial_Spectrum-Base; Corrected_Spectrum=Corrected_Spectrum(lp+1:lp+l);
Base=Base(lp+1:lp+l);
end