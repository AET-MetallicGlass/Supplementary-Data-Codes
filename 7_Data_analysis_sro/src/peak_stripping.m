function [Baseline, stripping]=peak_stripping(Spectrum,Window)
stripping=0;
y=sgolayfilt(Spectrum,0,Window);
n=length(Spectrum);
Baseline=zeros(n,1);
for i=1:1:n
   if Spectrum(i)>y(i)
       stripping=1;
       Baseline(i)=y(i);
   else
       Baseline(i)=Spectrum(i);
   end
       
    
    
end
end