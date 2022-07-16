% Close all figures and clear workspace
close all; clear; clc;
addpath(fullfile(pwd,'discrete.differencing'))
basedir='/Users/miladtavakolian/Library/CloudStorage/Box-Box/Purdue-IU-LiverTransplantEcho/LiverTransplant/Processing_6_21_22';

allBmodefiles=dir(fullfile(basedir,'Results_Bmode'));
alldGLSfiles=dir(fullfile(basedir,'Results_dGLS'));

%Cell with all folder names
Bmodenames={allBmodefiles.name};
dGLSnames={alldGLSfiles.name};

for ii= 4:length(allBmodefiles)

    %Find the ii'th Bmode foldername

    Bmodefoldername=allBmodefiles(ii).name;
%     Bmodefoldername = Bmodefoldername(~ismember({Bmodefoldername.name},{'.','..','.DS_Store'}));

    %See if that particular Bmodeefoldername is found in dGLSnames
    indx1=find(strcmp(dGLSnames,Bmodefoldername));

    %See if that particular CFIfoldername is found in dGLSnames
%     indx2=find(strcmp(dGLSnames,CFIfoldername));

    if ~isempty(indx1) %&& ~isempty(indx2)

        fprintf('\n Matching Files Found \n')

        [Bmodefoldername dGLSnames(indx1)]

        % Load the B-mode autosegmentation results
        filename_Bmode=fullfile([basedir,filesep,'Results_Bmode',filesep,Bmodefoldername,'/A4C Bmode/'],'*.mat');
        matfiles=dir(filename_Bmode);
        sgmnt = load(fullfile([basedir,filesep,'Results_Bmode',filesep,Bmodefoldername,'/A4C Bmode/'],matfiles.name));

        % Load the B-mode strain results
        filename_dGLS=fullfile([basedir,filesep,'Results_dGLS',filesep,dGLSnames{indx1},'/A4C Bmode/'],'*.mat');
        matfiles=dir(filename_dGLS);
        strn = load(fullfile([basedir,filesep,'Results_dGLS',filesep,dGLSnames{indx1},'/A4C Bmode/'],matfiles.name));

        %%
        % Determine sample indices using the function determine_number_of_beats
        % where
        %   eprime_v1 - first mitral annulus velocity from the structure sgmnt
        %   eprime_v2 - second mitral annulus velocity from the structure sgmnt
        %   dt  -   temporal resolution of the data
        %   HR  -   Heart rate, if it exists
        try
        if isfield(sgmnt.data.image.information,'HeartRate')
            %     [nBeats,sampInds] = determine_number_of_beats(sgmnt.eprime_v1,...
            %         sgmnt.eprime_v2,sgmnt.data.dt,sgmnt.data.image.information.HeartRate);
            [nBeats,sampInds] = determine_number_of_beats_v2(sgmnt.lvLength,...
                sgmnt.data.dt,sgmnt.data.image.information.HeartRate);
        else
            %     [nBeats,sampInds] = determine_number_of_beats(sgmnt.eprime_v1,...
            %         sgmnt.eprime_v2,sgmnt.data.dt,[]);
            [nBeats,sampInds] = determine_number_of_beats_v2(sgmnt.lvLength,...
                sgmnt.data.dt,[]);
        end
        
        C = strsplit(sgmnt.seg_inputs.patName,{'_'});

        sgmnt.laVolumeMOD2 = filloutliers(...
            filloutliers(filloutliers(sgmnt.laVolumeMOD,'pchip','quartiles'),...
            'pchip','movmedian',9),'pchip','movmedian',7);
% 
%          sgmnt.laVolumeMOD3 = filloutliers(...
%             filloutliers(filloutliers(sgmnt.laVolumeMOD,'spline','quartiles'),...
%             'spline','movmedian',9),'spline','movmedian',7);

                 sgmnt.laVolumeMOD3 = filloutliers(...
            filloutliers(filloutliers(sgmnt.laVolumeMOD,'linear','quartiles'),...
            'pchip','movmedian',9),'pchip','movmedian',7);

         sgmnt.lvVolumeMOD3 = filloutliers(...
            filloutliers(...
            filloutliers(sgmnt.lvVolumeMOD,'pchip','quartiles'),...
            'pchip','movmedian',9),...
            'pchip','movmedian',7);
        
        % Using the sample indices, find the relevant cardiac biomarkers from
        % B-mode data
        % Pre-allocate memory for autosegmentation LV & LA volumes
        bmrkr.lv.edv = zeros(nBeats,1);
        bmrkr.lv.esv = zeros(nBeats,1);
        bmrkr.la.edv = zeros(nBeats,1);
        bmrkr.la.esv = zeros(nBeats,1);
        for n = 1:nBeats
            bmrkr.lv.edv(n) = sgmnt.lvVolumeMOD(sampInds(n));
            bmrkr.lv.esv(n) = min(sgmnt.lvVolumeMOD(sampInds(n):sampInds(n+1)-1));
            bmrkr.la.edv(n) = sgmnt.laVolumeMOD3(sampInds(n));
            bmrkr.la.esv(n) = max(sgmnt.laVolumeMOD3(sampInds(n):sampInds(n+1)-1));
        end
        
figure;hold on;
plot(linspace(0,sgmnt.seg_inputs.num_time-1,sgmnt.seg_inputs.num_time)...
    *sgmnt.seg_inputs.dt,sgmnt.laVolumeMOD(1:end,1),'--','LineWidth',3);

plot(linspace(0,sgmnt.seg_inputs.num_time-1,sgmnt.seg_inputs.num_time)...
    *sgmnt.seg_inputs.dt,sgmnt.laVolumeMOD2(1:end,1),'--','LineWidth',3);
% plot(linspace(0,sgmnt.seg_inputs.num_time-1,sgmnt.seg_inputs.num_time)...
%     *sgmnt.seg_inputs.dt,sgmnt.laVolumeMOD3(1:end,1),'-.','LineWidth',2);
plot(linspace(0,sgmnt.seg_inputs.num_time-1,sgmnt.seg_inputs.num_time)...
    *sgmnt.seg_inputs.dt,sgmnt.laVolumeMOD3(1:end,1),'-','LineWidth',2);

grid on; set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold');
% leg = legend('LA Volume','LA Volume_{filloutliers}','LA Volume_{new filloutliers}'); 
% leg = legend('LA Volume','LA Volume2','LA Volume3','LA Volume4');
leg = legend('LA Volume','LA Volume_{pchip filloutliers}','LA Volume_{linear filloutliers}');
set(leg,'Location','best'); xlabel('Time, second','FontSize',18,'FontWeight','bold');
ylabel('LA Volume (Method of Disk), mL','FontSize',18,'FontWeight','bold');
title(C(1:2)); hold off
savedir = fullfile(basedir,filesep,'LVA_Size_Plots',filesep,'LA',filesep);
print(gcf,'-dpng',[savedir,sgmnt.seg_inputs.patName(1:end-0),'_LA.png'],'-r200'); 
% savefig(gcf,[savedir,sgmnt.seg_inputs.patName(1:end-0),'_LA.fig']);
close(gcf);
figure;hold on;
plot(linspace(0,sgmnt.seg_inputs.num_time-1,sgmnt.seg_inputs.num_time)...
    *sgmnt.seg_inputs.dt,sgmnt.lvVolumeMOD(1:end,1),'--','LineWidth',3);
plot(linspace(0,sgmnt.seg_inputs.num_time-1,sgmnt.seg_inputs.num_time)...
    *sgmnt.seg_inputs.dt,sgmnt.lvVolumeMOD3(1:end,1),'-.','LineWidth',2);
grid on; set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold');
leg = legend('LV Volume','LV Volume_{filloutliers}'); set(leg,'Location','best');
xlabel('Time, second','FontSize',18,'FontWeight','bold');
ylabel('LV Volume (Method of Disk), mL','FontSize',18,'FontWeight','bold');
title(C(1:2)); hold off
savedir = fullfile(basedir,filesep,'LVA_Size_Plots',filesep,'LV',filesep);
print(gcf,'-dpng',[savedir,sgmnt.seg_inputs.patName(1:end-0),'_LV.png'],'-r200'); 
% savefig(gcf,[savedir,sgmnt.seg_inputs.patName(1:end-0),'_LV.fig']);
close(gcf);
        % Compute stroke volume and ejection fraction
        bmrkr.lv.sv = abs(bmrkr.lv.edv-bmrkr.lv.esv);
        bmrkr.lv.ef = bmrkr.lv.sv./bmrkr.lv.edv*100;
        bmrkr.la.sv = abs(bmrkr.la.esv-bmrkr.la.edv);
        bmrkr.la.ef = bmrkr.la.sv./bmrkr.la.esv*100;
        % Create average mitral annulus velocity signal
        temp.MAV    = 0.5*sgmnt.eprime_v1 + 0.5*sgmnt.eprime_v2;
        % Pre-allocate memory for mitral annulus measurements
        bmrkr.lv.sprime = zeros(nBeats,1);
        bmrkr.lv.eprime = zeros(nBeats,1);
        bmrkr.lv.aprime = zeros(nBeats,1);
        for n = 1:nBeats
            temp.mav = -temp.MAV(sampInds(n):sampInds(n+1)-1)*100;
            % Compute sprime directly from the signal
            bmrkr.lv.sprime(n) = max(temp.mav);
            % Find in peaks from mitral annular velocity waveform
            [pks,locs]  =   findpeaks(temp.mav);
            % Reject peaks less than the standard deviation of the current MAV
            % signal
            locs(pks < mean(temp.mav(temp.mav >0))) = [];
            pks(pks < mean(temp.mav(temp.mav >0))) = [];
            % Find the distance between peaks
            dlocs = diff(locs);
            % Set rejection threshold
            thrsh = round(0.1*numel(sampInds(n):sampInds(n+1)-1));
            % Reject peaks
            pks = pks(logical([1 (dlocs(:)' > thrsh)]));
            [pks,inds]     =   sort(pks,'descend');
            if numel(inds) == 1
                inds = cat(1,inds,2);
                pks  = cat(1,pks,temp.mav(end));
            end
            inds = sort(inds(1:2),'ascend');
            % First indexed peak is eprime, second is aprime
            bmrkr.lv.eprime(n) = pks(inds(1));
            bmrkr.lv.aprime(n) = pks(inds(2));
        end

        % Pre-allocate memory for strain and strain rate measurements
        bmrkr.lv.GLSpeak    = zeros(nBeats,1);
        bmrkr.lv.GLSRs      = zeros(nBeats,1);
        bmrkr.lv.GLSRe      = zeros(nBeats,1);
        bmrkr.la.GLSpeak    = zeros(nBeats,1);
        bmrkr.la.GLSRs      = zeros(nBeats,1);
        bmrkr.la.GLSRe      = zeros(nBeats,1);
        for n = 1:nBeats
            % Compute the max strains
            bmrkr.lv.GLSpeak(n) = -(min(strn.LV.dGLS_corrected(sampInds(n):sampInds(n+1)-1)-strn.LV.dGLS_corrected(sampInds(n))));
            bmrkr.la.GLSpeak(n) =  -(min(strn.LA.dGLS_corrected(sampInds(n):sampInds(n+1)-1)-strn.LA.dGLS_corrected(sampInds(n))));
            % Compute the peak diastolic strain rate
            [bmrkr.lv.GLSRe(n),indLV] =  max(strn.LV.dGLSr_corrected(sampInds(n):sampInds(n+1)-1));
            [bmrkr.la.GLSRe(n),indLA] =  max(strn.LA.dGLSr_corrected(sampInds(n):sampInds(n+1)-1));
            % Compute the peak systolic strain rate
            bmrkr.lv.GLSRs(n) = -min(strn.LV.dGLSr_corrected(1:indLV));
            bmrkr.la.GLSRs(n) =  -min(strn.LA.dGLSr_corrected(1:indLA));
        end
        % Determine sample indices using the function determine_number_of_beatse
        % where
        %   eprime_v1 - first mitral annulus velocity from the structure sgmnt
        %   eprime_v2 - second mitral annulus velocity from the structure sgmnt
        %   dt  -   temporal resolution of the data
        %   HR  -   Heart rate, if it exists
%         if isfield(cfi.data.image.information,'HeartRate')
%             [nBeats,sampInds] = determine_number_of_beats(cfi.data.eprime_v1,...
%                 cfi.data.eprime_v2,cfi.data.dt,cfi.data.image.information.HeartRate);
%         else
%             [nBeats,sampInds] = determine_number_of_beats(cfi.data.eprime_v1,...
%                 cfi.data.eprime_v2,cfi.data.dt,[]);
%         end
%         end
        save(fullfile(basedir,'Postprocessing_corrected',sprintf('%s.mat',Bmodefoldername)));
        catch
        end
    else
        fprintf('\n No Match Foud \n')
        continue;
    end
end

