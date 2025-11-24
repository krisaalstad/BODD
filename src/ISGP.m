function [thetaout,Neff]=ISGP(Xin,Yin,sigy,hpri)
        %% Importance Sampling based GP hyperparameter MAP optimization
        Nd=size(Xin,2);
        Ns=size(Xin,1);
        
        Ne=1e4; % Ensemble members 

        dologit=0; % Assuming 0,1 bounds for now
        if dologit
            % Transformations for parameters
            glogit=@(x,a,b) log(((x-a)./(b-a))./(1-(x-a)./(b-a))); % Transforms from physical to unbounded (Gaussian) space
            %gexpit=@(xt,a,b) a+(b-a)./(1+exp(-xt)); % Transforms back from unbounded to physical space
            lsf=glogit(hpri.sfm,0,1)+hpri.sfs.*randn(Ne,1);
            lls=glogit(hpri.lsm,0,1)+hpri.lss.*randn(Ne,1);
        else
            lsf=log(hpri.sfm)+hpri.sfs.*randn(Ne,1); % Log of signal std prior
            lls=log(hpri.lsm)+hpri.lss.*randn(Ne,Nd);% Log of lengthscales prior
        end

        %lmlh=-0.5.*y'*alpha-sum(LLi);

        %tic;
        K=zeros(Ne,Ns,Ns); % Roughly 100 by 100 by 100
        for d=1:Nd
            Xd=Xin(:,d);
            %Kd=zeros(Ne,Ns,Ns); % Roughly 100 by 100 by 100
            Dsq=((Xd-Xd').^2);
            Kd=repmat(Dsq,[1 1 Ne]);
            Kd=permute(Kd,[3 1 2]);
            Kd=Kd./(exp(lls(:,d)).^2);
            K=K+Kd;
        end
        %disp('set up covs')
        %toc;
        %tic;
        K=(exp(lsf).^2).*exp(-0.5.*K);
        K=permute(K,[2 3 1]);
        lml=zeros(Ne,1); % Log marginal likelihood
        R=(sigy.^2).*eye(Ns);
        for j=1:Ne
            Cj=K(:,:,j)+R;
            Lj=chol(Cj,'lower');
            aj=(Lj')\(Lj\Yin);
            llj=-0.5.*(Yin'*aj)-trace(log(Lj)); % Ignore consant term
            lml(j)=llj;
        end
        %toc;
        %disp('ens chol')
        
        lmax=max(lml);
        lml=lml-lmax;
        w=exp(lml);
        w=w./sum(w);

        %{
        close all;
        histogram(w)
        k=waitforbuttonpress();
        %}
        Neff=1/sum(w(:).^2);
        
        here=find(w==max(w(:)),1,'first');
        if dologit
            gexpit=@(xt,a,b) a+(b-a)./(1+exp(-xt)); % Transforms back from unbounded to physical space
            ls=gexpit(lls,0,1)';
            sf=gexpit(lsf,0,1)';
        else
            ls=exp(lls)';
            sf=exp(lsf)';
        end
        thetaout=[ls(:,here); sf(here)];


    end
