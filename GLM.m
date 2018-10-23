% Performs a GLM on task phases as predictors and brain area activation as
% output variable
% calculates max activation (voxel values) within each brain region, then takes the avge activation b/w participants in each region 
% per experimental stage (reg value)

% TODO: 
% check if error term needs to be divided by sqrt(n) or if std() already takes care of that through normalization
% (look up second parameter value)
function general_linear_model()
close all;
startpath='/zifnas/polina.litvak/'
%addpath(startpath)
addpath(genpath(startpath))
data_path = [startpath 'output/01/data_4zs_18/']
output_path=[data_path 'GLM_analysis/'];

%if ~exist(output_path, 'dir') 
    mkdir(output_path)

    files=dir([data_path '*.mat']); nfiles=size(files,1);

    betas=[];
    group='C';
    idx_group_split=1; % remember where in the 3rd dim we switch from HC to SCZ
    %betas = zeros(ts_file.PLRNN.regs+1*size(Y,2), 4);
    for k= 1:nfiles
       file=files(k).name;
       vpn=file(end-12:end-8);

       if group ~= vpn(1)
          group = vpn(1);
          idx_group_split=k ;%remember from where group switches
       end

       ts_file=load(strcat(data_path,file));
       Y=ts_file.PLRNN.data';
       Yhat=Y;

       %Z transform all the voxel measurements
       for i=1:size(Y,2)
            Yhat(:,i)=(Y(:,i)-mean(Y(:,i)))/std(Y(:,i));
       end
       regs=ts_file.PLRNN.regs;
       regs_codes=[1 2 3 7 8 9];

       if size(Y,1) ~= size(regs,2)
           fprintf('**************** file = %s has mismatched data and regs size ************', file);
           regs=regs(1:size(Y,1));
       end

       T=length(Y);
       M=length(regs_codes);

       X=zeros(T,M+1);
       X(:,1)=ones(T,1);

       for m=1:M
          log_regs=regs == regs_codes(m); %X(2)=1=control,X(3)=2=training,X(4)=3=testing,X(5)=7=delay_train_test,X(6)=8=delay_cont_train,X(7)=9=score presnt.
          %fprintf('**************** k = %f\n', k);
          X(:,m+1)=log_regs'; 
       end
       %B=(X'*X)^-1*X'*Y; %solve for parameters
       B=pinv(X'*X)*X'*Yhat;
       betas=cat(3,betas,B);   %concat the betas of all hc and scz respectively
    end

    betas=betas(2:end,:,:); %disregard first row - this is the offset

    % get unique labels and index of next value of ROI
    % avge voxel values within region
    [label, index] = unique(ts_file.PLRNN.ROI.labs');
    label=strrep(label, '_', ' ');
    label=strrep(label, 'rmR', 'R');
    label=strrep(label, 'rmL', 'L');
    label{7}='Zs';

    for z=1:size(betas,3) % loop over all participants and collect max activation per brain region for each experiment stage
        for k=1:size(betas,1) % various experiment stages
           for m=1:length(index)-1 % brain region labels
                l_bound=index(m); % the breakdown of voxels per brain region
                r_bound=index(m+1)-1;
                max_beta_in_region(k,m,z)=max(abs(betas(k,l_bound:r_bound,z)));
           end
           l_bound=index(m+1);
           max_beta_in_region(k,m+1,z)=max(abs(betas(k,l_bound:end,z)));
        end
    end


    mean_max_activation(:,:,1)=mean(max_beta_in_region(:,:,1:idx_group_split-1),3); %HC
    mean_max_activation(:,:,2)=mean(max_beta_in_region(:,:,idx_group_split:end),3); %SCZ

    err(:,:,1)=std(max_beta_in_region(:,:,1:idx_group_split-1),0,3); %HC
    err(:,:,2)=std(max_beta_in_region(:,:,idx_group_split:end),0,3); %SCZ
    %keyboard


    %get mean max activation from latent states
    Zs_data=Zs_GLM();

    mean_max_activation(:,7,1)=Zs_data.mean_max_activation(:,1);
    mean_max_activation(:,7,2)=Zs_data.mean_max_activation(:,2);

    regions=ts_file.PLRNN.ROI.labs;
    save([output_path 'GLM_calc_output.mat'],'mean_max_activation','max_beta_in_region', 'Zs_data','idx_group_split', 'regions')
% else
%     dat=load([output_path 'GLM_calc_output.mat'])
%     max_beta_in_region=dat.max_beta_in_region
%     mean_max_activation=dat.mean_max_activation
%     Zs_data=dat.Zs_data
%     idx_group_split=dat.idx_group_split
%     regions=dat.regions
% end

thresh=.05;
plot_titles={'ctrl','train','test','WM delay','ctrl delay','score'}; 
% STATS
for stage=1:6
    for ROI=1:3
        ROI_name=string(label(ROI)); % 'L Hippocampus','L brodmann area 46','L brodmann area 9','R Hippocampus','R brodmann area 46','R brodmann area 9'
        stage_name=string(plot_titles(stage));
           
        %HC -left ROI vs Zs 
        [h,p,~,stats]=ttest2(squeeze(max_beta_in_region(stage,ROI,1:idx_group_split-1)),squeeze(Zs_data.max_beta_in_region(stage,1:idx_group_split-1))'); 
        txt1=sprintf('HC vs Zs,Left %s, %s ***** T(%d)=%2.2f, p=%2.3f',string(label(ROI)),stage_name,stats.df,stats.tstat,p);
        %HC - right ROI vs Zs
        [h,p,~,stats]=ttest2(squeeze(max_beta_in_region(stage,ROI+3,1:idx_group_split-1)),squeeze(Zs_data.max_beta_in_region(stage,1:idx_group_split-1))');
        txt2=sprintf('HC vs Zs,Right %s, %s ***** T(%d)=%2.2f, p=%2.3f',string(label(ROI+3)),stage_name,stats.df,stats.tstat,p);
        %SCZ - left ROI vs Zs
        [h,p,~,stats]=ttest2(squeeze(max_beta_in_region(stage,ROI,idx_group_split:end)),squeeze(Zs_data.max_beta_in_region(stage,idx_group_split:end))'); 
        txt3=sprintf('SCZ vs Zs,Left %s, %s ***** T(%d)=%2.2f, p=%2.3f',string(label(ROI)),stage_name,stats.df,stats.tstat,p);
        %SCZ - right ROI vs Zs
        [h,p,~,stats]=ttest2(squeeze(max_beta_in_region(stage,ROI+3,idx_group_split:end)),squeeze(Zs_data.max_beta_in_region(stage,idx_group_split:end))');
        txt4=sprintf('SCZ vs Zs,Right %s, %s ***** T(%d)=%2.2f, p=%2.3f',string(label(ROI+3)),stage_name,stats.df,stats.tstat,p);
        %HC - left ROI vs SCZ 
        
        
        [h,p,~,stats]=ttest2(squeeze(max_beta_in_region(stage,ROI,1:idx_group_split-1)),squeeze(max_beta_in_region(stage,ROI,idx_group_split:end)));
        txt5=sprintf('HC vs SCZ, Left %s, %s ***** T(%d)=%2.2f, p=%2.3f',string(label(ROI)),stage_name,stats.df,stats.tstat,p);
        if p < thresh
            disp(txt5);
        end
        %HC - right ROI vs SCZ 
        [h,p,~,stats]=ttest2(squeeze(max_beta_in_region(stage,ROI+3,1:idx_group_split-1)),squeeze(max_beta_in_region(stage,ROI+3,idx_group_split:end))); 
        txt6=sprintf('HC vs SCZ, Right %s, %s ***** T(%d)=%2.2f, p=%2.3f',string(label(ROI+3)),stage_name,stats.df,stats.tstat,p);
        if p < thresh
            disp(txt6);
        end
        
%         disp(txt1);
%         disp(txt2);
%         disp(txt3);
%         disp(txt4);
         
       
         disp('===================================================================');
    end
    disp('===================================================================');
    
   %HC vs SCZ Zs 
   [h,p,~,stats]=ttest2(squeeze(Zs_data.max_beta_in_region(stage,1:idx_group_split-1))',squeeze(Zs_data.max_beta_in_region(stage,idx_group_split:end))'); 
   txt=sprintf('HC vs SCZ in latent states, %s ***** T(%d)=%2.2f, p=%2.3f',stage_name,stats.df,stats.tstat,p);
   if p < thresh
            disp(txt);
   end
   %disp(txt); 
     
end

keyboard
   
thresh=.05;
for stage=1:6
   stage_name=string(plot_titles(stage));
   for other_stage=stage+1:6 
        o_stage_name=string(plot_titles(other_stage));
        % HC - compare betwee reg stages - e.g. copare ctrl to train,test,wm -delay,ctr-delay,score, then train to test, wm-delay,ctr-delay,score and so on

        [h,p,~,stats]=ttest(squeeze(Zs_data.max_beta_in_region(stage,1:idx_group_split-1))',squeeze(Zs_data.max_beta_in_region(other_stage,1:idx_group_split-1))');
        if p < thresh
            txt=sprintf('within HC group, Zs, %s vs %s ***** T(%d)=%2.2f, p=%2.3f',stage_name,o_stage_name,stats.df,stats.tstat,p);
            disp(txt);

        end

   end
   disp('--------------');
end
keyboard
for stage=1:6
   stage_name=string(plot_titles(stage));
   for other_stage=stage+1:6 
        o_stage_name=string(plot_titles(other_stage));
        % HC - compare betwee reg stages - e.g. copare ctrl to train,test,wm -delay,ctr-delay,score, then train to test, wm-delay,ctr-delay,score and so on

        [h,p,~,stats]=ttest(squeeze(Zs_data.max_beta_in_region(stage,idx_group_split:end))',squeeze(Zs_data.max_beta_in_region(other_stage,idx_group_split:end))');
        if p < thresh
            txt=sprintf('within SCZ group, Zs, %s vs %s ***** T(%d)=%2.2f, p=%2.3f',stage_name,o_stage_name,stats.df,stats.tstat,p);
            disp(txt);

        end

   end
   disp('--------------');
end


%PLOTS

sreg=getLabels(regions); %get short region labels
nplots=3; rplots=2; n=1; ifig=1;
f1(ifig)=figure('color','white'); hold on; hold on
set(f1(ifig),'PaperOrientation','Portrait','PaperUnits','normalized','PaperPosition', [0 0 1 1]);

for ix=1:size(mean_max_activation,1) %loop over the experimental stages (regs), ouput subplot per stage      
    subplot(rplots,nplots,n); hold on; box on;
    x=1:size(mean_max_activation,1)+1;
    title(plot_titles(ix));
    
    bpcombined = [mean_max_activation(ix,:,1)', mean_max_activation(ix,:,2)'];
    b = bar(x, bpcombined, 'grouped');
    set(b(1), 'FaceColor','b', 'EdgeColor','b') ;%HC
    set(b(2), 'FaceColor','r', 'EdgeColor','r') ;%SCZ
    
    width = b.BarWidth;
    offset = ((width + 0.5) / 10); % 0.5 is approximate net width of white spacings per group
    
    % plot max activation per region for each participant
    for iy=1:size(index) %loop over brain regions    
        plot(iy-offset,squeeze(max_beta_in_region(ix,iy,1:idx_group_split-1)),'ob'); %HC
        plot(iy+offset,squeeze(max_beta_in_region(ix,iy,idx_group_split:end)),'or'); %SCZ
        hc=errorbar(iy-offset,mean_max_activation(ix,iy,1),err(ix,iy,1),'b');  
        scz=errorbar(iy+offset,mean_max_activation(ix,iy,2),err(ix,iy,2),'r');
    end
    
    %plot latent states activation
    plot(iy+1-offset,Zs_data.max_beta_in_region(ix,1:idx_group_split-1),'ob'); %HC
    plot(iy+1+offset,Zs_data.max_beta_in_region(ix,idx_group_split:end),'or'); %SCZ
   
    hc=errorbar(iy+1-offset,Zs_data.mean_max_activation(ix,1),std(Zs_data.max_beta_in_region(ix,1:idx_group_split-1)),'b');  
    scz=errorbar(iy+1+offset,Zs_data.mean_max_activation(ix,2),std(Zs_data.max_beta_in_region(ix,idx_group_split:end)),'r');
     
    legend([hc scz], {'HC' 'SCZ'});
    xticks([1 2 3 4 5 6 7]); xticklabels(label);
    set(gca,'XTickLabelRotation',45);
    n=n+1;
end
hold off;


%plot latent states HC againt SCZ for the different regs

nplots=2; rplots=1; n=1; ifig=1;
f1(ifig)=figure('color','white'); hold on; hold on
set(f1(ifig),'PaperOrientation','Portrait','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
for y=1:2
    x=1:size(mean_max_activation,1);
    title('Zs HC vs SCZ for different cognintive stages (Regs)');
    subplot(2,1,n); hold on; box on;

    b = bar(x, Zs_data.mean_max_activation(:,y));

    if y==1 %HC
        set(b, 'FaceColor','b', 'EdgeColor','b');
        title('Zs HC different cognintive stages (Regs)');
        for iy=1:size(mean_max_activation,1) %loop over regs  
            vals=Zs_data.max_beta_in_region(iy,1:idx_group_split-1);
            
            nstd=4;
            
            ox=vals>mean(vals)+nstd*std(vals) | vals<mean(vals)-nstd*std(vals);
            
            vals(ox)=[]; 
            
            plot(iy,Zs_data.max_beta_in_region(iy,1:idx_group_split-1),'or'); %HC
            hc=errorbar(iy,mean(vals),std(vals),'r'); 
        end
    else  %SCZ
        set(b, 'FaceColor','r', 'EdgeColor','r');
        title('Zs SCZ different cognintive stages (Regs)');
        for iy=1:size(mean_max_activation,1) %loop over regs  
            vals=Zs_data.max_beta_in_region(iy,idx_group_split:end);
            
            nstd=4;
            
            ox=vals>mean(vals)+nstd*std(vals) | vals<mean(vals)-nstd*std(vals);
            
            vals(ox)=[]; 
            
            plot(iy,vals,'ob'); %SCZ
            scz=errorbar(iy,mean(vals),std(vals),'b');
        end
    end

    width = b.BarWidth;
    offset = ((width + 0.5) / 10); % 0.5 is approximate net width of white spacings per group

    xticks([1 2 3 4 5 6]); xticklabels(plot_titles);
    set(gca,'XTickLabelRotation',45);
    n=n+1;
end;

%plot latent states HC againt SCZ for the different regs
figure();hold on; box on;
x=1:size(mean_max_activation,1);
title('Zs HC vs SCZ for different cognintive stages (Regs)');

bpcombined = [Zs_data.mean_max_activation(:,1), Zs_data.mean_max_activation(:,2)];
b = bar(x, bpcombined, 'grouped');
set(b(1), 'FaceColor','b', 'EdgeColor','b') ;%HC
set(b(2), 'FaceColor','r', 'EdgeColor','r') ;%SCZ

width = b.BarWidth;
offset = ((width + 0.5) / 10); % 0.5 is approximate net width of white spacings per group

% plot max activation per region for each participant
for iy=1:size(mean_max_activation,1) %loop over regs  
    plot(iy-offset,Zs_data.max_beta_in_region(iy,1:idx_group_split-1),'ob'); %HC
    plot(iy+offset,Zs_data.max_beta_in_region(iy,idx_group_split:end),'or'); %SCZ
    hc=errorbar(iy-offset,Zs_data.mean_max_activation(iy,1),std(Zs_data.max_beta_in_region(iy,1:idx_group_split-1)),'b');  
    scz=errorbar(iy+offset,Zs_data.mean_max_activation(iy,2),std(Zs_data.max_beta_in_region(iy,idx_group_split:end)),'r');
end
legend([hc scz], {'HC' 'SCZ'});
xticks([1 2 3 4 5 6]); xticklabels(plot_titles);
set(gca,'XTickLabelRotation',45);

keyboard


keyboard


clearvars -except betas nfiles
end
