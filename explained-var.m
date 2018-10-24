function explained_var()
    %comparing explained variance for a different number of states in terms of F-Test
    %excluding subjects whose MSE > 3 stddev (mean sq. error is calculated
    %b/w all subjects for a given state space model - e.g. all subj for 1z)
    % this version doesn't exclude the same people from states that are
    % being compared (when comparing 1z to 2z, there should be a union of
    % outlyers - explained_var_pairwise_outlyrs.m does this correctly)
    
    clear all; close all; clc
   
    path = '/zifnas/polina.litvak/output/09/triton_02/PLRNN/'; 
    addpath('/zifnas/polina.litvak/programs/toolbox/PLRNN2/PLRNN_fMRI_sp_simultan_AWCh')
   
    plot_titles={'1z','2z','3z','4z','5z','6z','7z', '8z','9z','10z'}; 
    
    subj_err_zs=[];
    outliers_index={};
    %if ~exist([path 'exp_var.mat'], 'file') 
    for i=1:10
        path_model=[path num2str(i) 'z/'];

        data=calcResiduals(path_model); 
        outliers_index{end+1} = data.outliers_index;
        err = data.err(:);
        vals{i}=err;

        subj_errors = data.err_per_subj';
        subj_err_zs=vertcat(subj_err_zs,subj_errors(:)');
     end

    % get sum of total error per number of states
    for i=1:length(vals)
        sum_error(i)=sum(vals{i},1)
        mse(i)=mean(cell2mat(vals(i)));
        stdv(i)=std(cell2mat(vals(i)));
    end 
    
    mean_err=[];
    mean_std=[];
    without_outlyrs={};
    for i=1:10
               
        nstd=3;
        ox=subj_err_zs(i,:)>mean(subj_err_zs(i,:))+nstd*std(subj_err_zs(i,:)) | subj_err_zs(i,:)<mean(subj_err_zs(i,:))-nstd*std(subj_err_zs(i,:));
       
        %x(ox)=[]; g(ox)=[];      
        temp=subj_err_zs(i,:);
        temp(ox)=[];
        
        num_subj(i)=length(temp); % number of remaining subject after unstable ones are exc

        
        txt=sprintf('nusm of excl subjects for  %d are %d', i,sum(ox));
        disp(txt);
        
        mean_err=vertcat(mean_err,mean(temp));
        mean_std=vertcat(mean_std,std(temp));
        without_outlyrs=vertcat(without_outlyrs,temp);
    end
    keyboard;
    
    % exp variance plot - total error as func of number of states
    figure; hold on;
    bar(1:10,mean_err)
    errorbar(1:10,mean_err,mean_std,'.')
    for i=1:10
        sbj_err=without_outlyrs(i);
        plot(i,cell2mat(sbj_err),'or');
    end
    title('1-step forward prediction');
    xticks(1:10);xticklabels({'1z','2z','3z','4z','5z','6z','7z','8z','9z','10z'});
    xlabel('# latent states'); ylabel('mean squared error');
      
    figure(); x=1:10;
    %plot(x,sum_error,'b','LineWidth',2);
    plot(x,mse,'b','LineWidth',2);
    xticks(x);xticklabels({'1z','2z','3z','4z','5z','6z','7z','8z','9z','10z'});
    set(gca,'XTickLabelRotation',30);
    xlabel('number of states'); ylabel('squared error');
    title('1-step forward prediction');
    
    % p-value comparison b/w sqrd error for full and restricted models
    for i=1:9 %iterate states
        err_full = vals{i+1};%full.err(:);
        err_restr = vals{i};%restr.err(:);

        ss_full = sum(err_full); %sum of squares
        ss_restr = sum(err_restr);

        zfull=i+1;%full.num_states; 
        zrestr=i;%restr.num_states;
        n=data.num_fmri_series; 

        %A W h B Gamma Sigma ---> count params in these matrices to get degrees of freedom
        q= zfull+zfull^2-zfull+ zfull+ zfull*n+ n + zfull; %num of params full model
        q= num_subj(i+1)*q;
        r= zrestr+zrestr^2-zrestr+ zrestr+ zrestr*n+ n + zrestr;
        r= num_subj(i)*r;
        dof1=q-r;
        dof2=length(err_restr)-q-1;

        ratio=((ss_restr-ss_full)/dof1)/(ss_full/dof2);
        pv=1-fcdf(ratio,dof1,dof2);

        txt=sprintf('F(%d,%d)=%2.3f, p=%2.3f',dof1,dof2,ratio,pv);
        disp(txt);
    end
    keyboard
    %save([path 'exp_var.mat'],'mean_max_activation','max_beta_in_region', 'Zs_data','idx_group_split', 'regions')
end

%--------------------------------------------------------------------------
function calc=calcResiduals(pat)
    %load in fMRI data-set
    str=['*Z*.mat'];
    files=dir([pat str]);
    nfiles=size(files,1);

    err_all={};
    ev_all=[];
    subjects={};
    
    %loop over different subjects
    for k=1:nfiles
        %load results of PLRNN estimation
        file=[ files(k).name];
        filename=file(1:end-4);
        vpn=file(end-23:end-19);
        dat=load([pat file]);
        X=dat.X';
        LL=dat.LLall{1}(end); %log likelihood per individual - unused atm
        
        W=dat.nets{1}.W;
        A=dat.nets{1}.A;
        h=dat.nets{1}.h;
        Ezi=dat.EzAll{1};
        
        [p,T]=size(Ezi);

        B_est=dat.nets{1}.B;
        M_est=dat.nets{1}.M;
        rp=dat.rp; %noise covariates

        %ind=PLRNN_getBestEstimate(dat);
        Ezi_plrnn=dat.EzAll{1}; %state estimate
        [p,T]=size(Ezi_plrnn);

        % do one step prediction
        Ezi_pred=zeros(p,T);
        for t=2:T
            zt=Ezi_plrnn(:,t-1); % use Zt-1 to predict Zt_hat ?
            zt_pred=A*zt+W*max(0,zt)+h;
            Ezi_pred(:,t)=zt_pred;
        end
%       2-delta t forward pred        
%       for t=3:T
%           zt-2=Ezi_plrnn(:,t-2); % use Zt-1 to predict Zt_hat ?
%           zt-1=A*zt-2+W*max(0,zt-2)+h;
%           zt_pred=A*zt-1+W*max(0,zt-1)+h;
%           Ezi_pred(:,t)=zt_pred;
%       end
        
        %Ezi_pred=Ezi_plrnn; 
    
        TR=dat.PLRNN.preprocess.RT; %scan repitition time
        M=size(Ezi_plrnn,1); %number of states
        H=getConvolutionMtx(M,T,TR);
        Ezi=reshape(Ezi_pred,1,numel(Ezi_pred));
        hZ=H*Ezi';
        hEzi=reshape(hZ,M,T);

        % get x_pred from Ezi_pred
        X_pred=(B_est*hEzi+M_est*rp')';
        
        %if abs(ev) < 1.5 %exclude subjects with unstable eigen values
            err_subj=(X-X_pred).^2;
            err_all{k}=err_subj;
            err_per_subj(k)= sum(sum(err_subj));
%           ev_all=vertcat(ev_all, ev);
            subjects{k}=string(vpn)+' ';
%       else
%           txt=sprintf('path %s, subject %s has ev=%2.3f', pat,string(vpn),ev);
%           disp(txt); 
%       end
    end
   
    err_all_vect=[];
    nstd=3;
    indx=find(err_per_subj>mean(err_per_subj)+nstd*std(err_per_subj) | err_per_subj<mean(err_per_subj)-nstd*std(err_per_subj))
    
    err_all(indx)=[];
    for i=1:size(err_all,2)
        temp=err_all(i);
        temp_mat=cell2mat(temp);
        err_all_vect=vertcat(err_all_vect,temp_mat(:));
    end
    calc.err_per_subj=err_per_subj;
    calc.err=err_all_vect;
    calc.num_states=p;
    calc.num_subjects=nfiles;
    calc.num_fmri_series=size(X,2);
    calc.subjects=subjects;
    calc.ev=ev_all;
    calc.outliers_index=indx;
end

%--------------------------------------------------------------------------
function H=getConvolutionMtx(no_states,no_trials,TR)
    M=no_states;
    T=no_trials;
    %get hrf for convolution
    hrf=spm_hrf(TR);

    %%%new version
    n=1;
    for i=1:length(hrf)
        zhrf(n)=hrf(i); n=n+1;
        for j=2:M
            zhrf(n)=0; n=n+1;
        end
    end
    Hconv=convmtx(zhrf',M*T);

    %Hconv=Hconv';
    H=Hconv(1:T*M,:);
end
