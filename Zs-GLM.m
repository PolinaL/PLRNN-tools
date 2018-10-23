%todo -move all constant properties into a separate class
function [Zs]=Zs_GLM()
    zifnas='/zifnas/polina.litvak';
    %pati=[zifnas '/output/01/_01/'];
    pati=[zifnas '/output/09/triton_02/PLRNN/4z/']
    %pati=[zifnas '/output/01/_04/PLRNN/'];

    %load in fMRI data-set
    noZ=size(Zs, 1); %noZ=4;
    files=dir([pati '*.mat']);
    nfiles=size(files,1);

    betas=[];
    group='C';
    idx_group_split=1; % remember where in the 3rd dim we switch from HC to SCZ

    %loop over different fMRI subjects
    for k=1:nfiles
        %load results of PLRNN estimation
        file=[ files(k).name];
        filename=file(1:end-4);
        vpn=file(end-23:end-19);
        
        if group ~= vpn(1)
           group = vpn(1);
           idx_group_split=k; %remember from where group switches
        end
        
        dat=load([pati file]);
        %ind=PLRNN_getBestEstimate(dat);
        
        Y=dat.EzAll{1};
        
        
        %Y=dat.EzAll{ind}'; %state estimate
        
        regs=dat.PLRNN.regs;
        regs_codes=[1 2 3 7 8 9];
        
        if size(Y,2) ~= size(regs,2);
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
           X(:,m+1)=log_regs' ;
        end
        B=pinv(X'*X)*X'*Y'; %solve for parameters, pinv=Moore-Penrose pseudo inverse 
        betas=cat(3,betas,B) ;  %concat the betas of all hc and scz respectively
    end
    betas=betas(2:end,:,:); %disregard first row - this is the offset

    
    for z=1:size(betas,3) % loop over all participants and collect max activation per state for each experiment stage
        for k=1:size(betas,1) % various experiment stages
            max_beta_in_region(k,z)=max(abs(betas(k,:,z))); %because we don't know which of the states (z1,z3,z3..) encode this brain region the most, we take the max
        end
    end

    mean_max_activation(:,1)=mean(max_beta_in_region(:,1:idx_group_split-1),2); %HC
    mean_max_activation(:,2)=mean(max_beta_in_region(:,idx_group_split:end),2); %SCZ
    
    Zs.mean_max_activation=mean_max_activation;
    Zs.max_beta_in_region=max_beta_in_region;
end
