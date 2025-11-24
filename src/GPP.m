function [Ypm,Yps]=GPP(Xt,Yt,Xp,kpars)
        %% Gaussian Process Prediction.
        % Perform Gaussian Process prediction using fixed hyperparameters.
        % Assumes an ARD-SE(RB) Kernel with associated parameters (kpars):
        % kpars.l= Lengthscales
        % kpars.sf = Signal standard deviation
        % kpars.sy = Observation noise standard deviation

        Nd=size(Xt,2); % Number of input dimensions
        Ntr=size(Xt,1);
        Npr=size(Xp,1);
        Ktt=zeros(Ntr,Ntr);
        Kpt=zeros(Npr,Ntr);
        Kpp=zeros(Npr,Npr); % Huge matrix... save as vector (squareform?)
        for d=1:Nd
            Xtd=Xt(:,d);
            Xpd=Xp(:,d);
            ld=kpars.l(d);
            Ktt=Ktt+(((Xtd-Xtd').^2)./ld^2);
            Kpt=Kpt+(((Xpd-Xtd').^2)./ld^2);
            Kpp=Kpp+(((Xpd-Xpd').^2)./ld^2);
        end
        sf=kpars.sf;
        sy=kpars.sy;
        
        Ktt=(sf^2).*exp(-0.5.*Ktt);
        Kpt=(sf^2).*exp(-0.5.*Kpt); % ktt is contained in kpt, can speed up by
        %slicing with tinds
        Kpp=(sf^2).*exp(-0.5.*Kpp); % Both the above contained in this...
        C=Ktt+(sy.^2).*eye(Ntr);

        % Based on Alg 4.4 in Hennig PNum book
        [R,pdflag]=chol(C,'upper'); % Consider using nearest SPD catch on chol
        if any(pdflag)
            %size(C)
            %disp(C)
            %disp(sf)
            %disp(sy)
            %disp(kpars.l)
            [Lm, DMCm, Pm]=modchol_ldlt(C);
            Cspd=Pm'*Lm*DMCm*Lm'*Pm;
            R=chol(Cspd);
        end
        al=R\((R')\Yt);
        Ypm=Kpt*al;
        L=Kpt/R;
        %toc;
        Yps=sqrt(diag(Kpp)- sum(L.*L,2)); 
        
        %max(abs(fsout-fsout2))
        %cholt=toc;
        %fprintf('\n chol time %4.2f \n',cholt);
        

    end
