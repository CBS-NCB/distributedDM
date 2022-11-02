function [map_hor,map_vert,avrAmp] = func_MakeAngleMaps(CmplxMaps,visual_field,Threshold)
avrAmp = zeros(size(CmplxMaps{1}{1}));
maxAmp = zeros(size(CmplxMaps{1}{1}));
for j=1:2
    for i=1:2
        temp = mat2gray(log(abs(CmplxMaps{j}{i})+1));
        avrAmp = avrAmp + temp/4;
        maxAmp = max(maxAmp,temp);
    end
end
roibw = double(avrAmp>Threshold);
%%
for j=1:2
    for i=1:2
        temp = mat2gray(log(abs(CmplxMaps{j}{i})+1));
        tempPhi = angle((CmplxMaps{j}{i}));
        Cmplx_AmpCut{j}{i} = temp.*exp(1i*double(roibw).*tempPhi);
    end
end

%%
ang0 = angle(Cmplx_AmpCut{1}{1});
ang2 = angle(Cmplx_AmpCut{1}{2});
ang1 = angle(Cmplx_AmpCut{2}{1});
ang3 = angle(Cmplx_AmpCut{2}{2});

%% Segmentation (Ian Nauhaus)
delay_hor  = angle( exp(1i*ang0)+exp(1i*ang2) );
delay_vert = angle( exp(1i*ang1)+exp(1i*ang3) );
delay_hor  = delay_hor  + (pi/2)*(1-sign(delay_hor));
delay_vert = delay_vert + (pi/2)*(1-sign(delay_vert));
map_hor  = 0.5*(angle(exp(1i*(ang0-delay_hor)))  - angle(exp(1i*(ang2-delay_hor))));
map_vert = 0.5*(angle(exp(1i*(ang1-delay_vert))) - angle(exp(1i*(ang3-delay_vert))));
map_hor  = (map_hor  /(2*pi)) * visual_field(1);
map_vert = (map_vert /(2*pi)) * visual_field(2);
