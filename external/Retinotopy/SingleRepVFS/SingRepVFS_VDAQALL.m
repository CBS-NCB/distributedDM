global VDAQ
clear CmplxMaps

%% here, Load the data and do basic prep using Main_WideField.m for all the iexps
for iAzm = 1:p.nrepeats % for M150908-2, trial 3 and 16 is crap.
    for iAlt = 1:p.nrepeats
        for idx  = 1:length(Expts)
            iexp  = Expts(1,idx);
            VDAQ = VDAQ_ALL(idx);
            %set rep tensor as VDAQ.tensor
            for istim=1:p.nstim,
                VDAQ.tensor{istim} = fliplr(rot90(VDAQ_ALL(idx).alldata{istim,iAzm},3)); end
            Main_WideField_SingleRep; % Run the initial analyses
            [AbsMaps{idx}, AngleMaps{idx}, CmplxMaps(iAzm).C{idx}] = TensorFrequency([], myfreqs );
        end
        
        % Apply eye position (meridians) with respect to the screen
        EyeY = -17;
        EyeX = -10;
        
        % Direct and Reverse maps for altitude and azimuth
        ang0 = angle(CmplxMaps(iAzm).C{1}{1});
        ang2 = angle(CmplxMaps(iAzm).C{1}{2});
        ang1 = angle(CmplxMaps(iAzm).C{2}{1});
        ang3 = angle(CmplxMaps(iAzm).C{2}{2});
        
        % Segmentation (Ian Nauhaus)
        %Find delay as the angle between the vectors
        delay_hor  = angle( exp(1i*ang0)+exp(1i*ang2) );
        delay_vert = angle( exp(1i*ang1)+exp(1i*ang3) );
        
        %The delay always pushes the vectors counter clockwise. (This is simply mod(val, pi), MA20150727)
        delay_hor  = delay_hor  + (pi/2)*(1-sign(delay_hor));
        delay_vert = delay_vert + (pi/2)*(1-sign(delay_vert));
        
        % additional corrections
        angind = delay_hor < pi/2;
        delay_hor(angind) = abs(delay_hor(angind) - pi);
        angind = delay_vert < pi/2;
        delay_vert(angind) = abs(delay_vert(angind) - pi);
        
        %Use delay vector to calculate retinotopy.
        map_hor  = 0.5*(angle(exp(1i*(ang0-delay_hor)))  - angle(exp(1i*(ang2-delay_hor))));
        map_vert = 0.5*(angle(exp(1i*(ang1-delay_vert))) - angle(exp(1i*(ang3-delay_vert))));
        
        map_hor  = (map_hor  /(2*pi)) * visual_field(1);
        map_vert = (map_vert /(2*pi)) * visual_field(2);
        
        % get the visual field sign map
        [VFS(:,:,iAzm), VFS_thr(:,:,iAzm)] = getVisSign(map_hor, map_vert);
        
        % Save the file for segmentation analysis
        % save([DIRS.SaveDir 'AnalyzedRet'],'CmplxMaps','EyeX', 'EyeY','visual_field')
    end
end

%%
figure
subplot(221); imagesc(mean(VFS,3)), axis image; colorbar; title('mean VFS')
subplot(222); imagesc(std(VFS,[],3)), axis image; colorbar; title('sd VFS')
zScore = (VFS - mean(VFS(:)))/std(VFS(:));
subplot(223);imagesc(mean(zScore,3)./std(VFS,[],3)), axis image; colorbar; title('mu/sd')
subplot(224);imagesc(mean(zScore,3)), axis image; colorbar; title('z VFS')

% ZZ = zscore(VFS+1,0,3);
% figure; imagesc(mean(ZZ,3)); axis image; colorbar

%%
save([DIRS.SaveDir,'ComplexMaps_AllRep'],'CmplxMaps','-v7.3')
% save([DIRS.SaveDir,'VDAQALL'],'VDAQ_ALL','-v7.3')


















