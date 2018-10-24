 %goes through all fixed points in zs and picks out only the ev that
 %are b/w 1 and -1 : max|eig(A+W.*D)| < 1 -----> these are the stable fixed points!
 %then plots all eigen values for the identified fixed points

function ev_plot()
    clear all; close all; clc
   
    path = '/zifnas/polina.litvak/output/09/triton_02/PLRNN/'; 
    addpath('/zifnas/polina.litvak/programs/toolbox/PLRNN2/PLRNN_fMRI_sp_simultan_AWCh')
    
    plot_titles={'1z','2z','3z','4z','5z','6z','7z', '8z','9z','10z'}; 
    nplots=1; rplots=1; %n=1; ifig=1;
   
    for i=3:10
        path_model=[path num2str(i) 'z/']; %iterate directories and plot for the diff number of latent states 
        data=calcEV(path_model); 
        ev_hc=data.ev_hc; %all fixed points, without restrictions
        ev_scz=data.ev_scz;
        maxabs_ev_hc=data.maxabsev_hc; %i use this for distrib. plot generation
        maxabs_ev_scz=data.maxabsev_scz; %one outlyer excl
        subjects=data.subjects;
   
        %subplot(rplots,nplots,iplot); hold on; box on
        figure(); hold on;
        ax = gca;                        % gets the current axes
        ax.XAxisLocation = 'origin';     % sets them to zero
        ax.YAxisLocation = 'origin';     % sets them to zero
       
        %unit circle
        n = 1000; %Define number of points on circle
        theta = linspace(0, 2*pi, n);
        j = sqrt(-1);
        x = @(theta) cos (theta) + j*sin (theta);
        plot(real(x(theta)), imag(x(theta)))

        hc_imag=imag(ev_hc);
        scz_imag=imag(ev_scz);

        hc_real=real(ev_hc);
        ind_max_hc=(hc_real==max(hc_real))
        mx_hc=find(ind_max_hc==1)

        scz_real=real(ev_scz);
        ind_max_scz=(scz_real==max(scz_real))
        mx_scz=find(ind_max_scz==1)

        eigs_hc=complex(hc_real(mx_hc),hc_imag(mx_hc))
        eigs_scz=complex(scz_real(mx_scz),scz_imag(mx_scz))
      
        hold off;
        
        filename=sprintf('/home/polina.litvak/plots/eig.png');
        print(gcf,filename,'-dpng','-r600'); 
        
        % fit distributions %    
        nbins=15; series1 = maxabs_ev_hc; series2 = maxabs_ev_scz;
        [series1] = hist(series1);
        [series2] = hist(series2);
        
        figure (); hold on;
        h1=histfit(maxabs_ev_hc,nbins,'kernel');
        set(h1(2),'color','b')
        delete(h1(1));

        h2=histfit(maxabs_ev_scz,nbins, 'kernel');
        yt2=get(gca,'YTick');
        ytickformat('%.2f');
        set(h2(2),'color','r')
        delete(h2(1));
        
        xlim([0,3]);
        xlabel('|eigenvalue|')
        ylabel('rel. frequency')
        set(gca, 'YTick', yt2, 'YTickLabel', sprintf('%.2f\n',yt2/numel(maxabs_ev_scz)))
        
        % Put up legend.
        legend1 = sprintf('  HC mean=%.2f, SE=%.3f ', mean(maxabs_ev_hc),std(maxabs_ev_hc));
        legend2 = sprintf('SCZ mean=%.2f, SE=%.3f ', mean(maxabs_ev_scz),std(maxabs_ev_scz));
        legend({legend1, legend2});   
        
%       [h,p,~,stats]=ttest2(maxabs_ev_hc,maxabs_ev_scz);
%       txt=sprintf('T(%d)=%2.2f, p=%2.3f',stats.df,stats.tstat,p);
%       [h,np] = kstest(maxabs_ev_hc) %normality check HC
%       [h,np] = kstest(maxabs_ev_scz) %normality check SCZ

        %[h,p,stats]=ranksum(maxabs_ev_hc,maxabs_ev_scz);
        [h,p,stats]=kstest2(maxabs_ev_hc,maxabs_ev_scz)
        txt=sprintf('Z=%2.2f, p=%2.3f',stats,p);
        
        filename=sprintf('/home/polina.litvak/Pictures/jun12/eig-distrib.png');
        print(gcf,filename,'-dpng','-r600'); 
    end
end

function calc=calcEV(pat)
    %load in fMRI data-set
    str=['*Z*.mat'];
    files=dir([pat str]);
    nfiles=size(files,1);
    subjects=cell(1,nfiles);
    ev_hc=[];ev_scz=[];out_hc=[];out_scz=[];maxabs_ev_hc=[];maxabs_ev_scz=[];stbl_ev_hc=[];stbl_ev_scz=[];subjects=[];
    
    %loop over different subjects
    for k=1:nfiles
        %load results of PLRNN estimation
        file=[ files(k).name];
        filename=file(1:end-4);
        vpn=get_vpn(file);
        dat=load([pat file]);
             
        %keyboard
        W=dat.nets{1}.W;
        A=dat.nets{1}.A;
        h=dat.nets{1}.h;
        
        %ev=eig(A+W);
        
        %go through all fixed points in zs and pick out only the ev that
        %are b/w 1 and -1 : max|ev| < 1 -----> these are the stable fixed points!
        
        [zs,ev,R]=AllFP_PLRNN(A,W,h);
        ev_mat=cell2mat(ev)
        ind_stbl= find(sum(abs(real(ev_mat)) < 1.0) == length(h)); %length(h) is the num of states, find eigvals for each fp are < 1
        
        % store stable fps
        if length(ind_stbl) > 0
            ev_st=ev_mat(:,ind_stbl);
            if strfind(vpn,'C')== 1
                stbl_ev_hc=[stbl_ev_hc max(abs(real(ev_st)))];
            else
                stbl_ev_scz=[stbl_ev_scz max(abs(real(ev_st)))];
            end
        end
        
        %all fixed points now, not just the stable ones 
        ind= find(sum(abs(real(ev_mat)) < 100.0) == length(h)); %find eigvals for each fp 
        if length(ind) > 0
            ev_st=ev_mat(:,ind);
            if strfind(vpn,'C')== 1
                ev_hc=[ev_hc ev_mat];
                maxabs_ev_hc=[maxabs_ev_hc max(abs(real(ev_st)))];
            else
                ev_scz=[ev_scz ev_mat];
                maxabs_ev_scz=[maxabs_ev_scz max(abs(real(ev_st)))];
            end
            subjects=[subjects string(vpn)];
        end
    end
    num_stbl_hc=length(stbl_ev_hc);num_stbl_scz=length(stbl_ev_scz);
    num_total_hc=length(ev_hc);num_total_scz=length(ev_scz);
    
    percent_stable_fp_hc=num_stbl_hc/num_total_hc*100
    percent_stable_fp_scz=num_stbl_scz/num_total_scz*100
    percent_stbl_fp_all= (num_stbl_hc+num_stbl_scz)/(num_total_hc+num_total_scz)*100
    
    calc.subjects=subjects;
    calc.ev_hc=ev_hc;
    calc.ev_scz=ev_scz;
    calc.maxabsev_hc=maxabs_ev_hc;
    calc.maxabsev_scz=maxabs_ev_scz;
end


function [zs,ev,R]=AllFP_PLRNN(A,W,h)
  % z_t = A z_t-1 + W phi(z_t-1,0) + h

  M=length(h);

  % move through all hyper-cubes
  R=zeros(1,2^M);
  n=0; zs=cell(1); ev=cell(1);
  for i=0:2^M-1
      d=dec2bin(i)-48; d=[zeros(1,M-length(d)) d];

      D=repmat(d,M,1);

      W0=W.*D;
      z1=-(A+W0-eye(M))^-1*h;
      d1=(z1>0);
      R(i+1)=sum(abs(d1-d'));
      if R(i+1)==0 
          n=n+1; zs{n}=z1;
          ev{n}=eig(A+W0);
      end;
  end;
end

function vpn=get_vpn(file)
    if isempty(strfind(file, 'Z10'))
        vpn=file(end-23:end-19);
    else
        vpn=file(end-32:end-28);
    end
end
