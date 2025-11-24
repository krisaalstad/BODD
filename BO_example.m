clearvars;

% Example Bayesian optimization of a degree day model 
% This could be CryoGrid just it would take longer to run.
% Seems like this generally converges in less than 20 iterations,
% so it could be a nice way to calibrate parameters without the need
% for gradients (but also without uncertainty quantification) in an
% automatic way to avoid manual tuning.
% Similar ideas could also be used to generate spatial maps by exploring
% uncertain regions with a GP either as an alternative to or a
% post-processing step of clustering.

% Since it is easier to visualize this is just for 1 parameter, but
% this could also be a multidimensional parameter space (since the model 
% parameters to be calibrated are the inputs to a Gaussian process or neural network)

addpath('src');

% Some example data from ERA5 near Ny-Ålesund
load('dailyTP.mat');
t=f.t;
P=f.P;
T=f.T;
clear f;
tstart=min(t);
tend=max(t);

%% Synthetic truth run
% Also to generate (synthetic) observations
ctrue=1; % Snowfall faCtor
a=4;
b=0;
Tr=2;
Ts=-2;

Dtrue=DDM(t,T,P,a,b,ctrue,Tr,Ts);

tf=0.01;
Nt=numel(t);
Nk=round(tf*Nt);
rinds=sort(randsample(Nt,Nk));
thesetrain=zeros(Nk,1,'logical');
thesetrain(rinds)=1;
% Synthetic observations
sigy=10;
Dobs=Dtrue(thesetrain)+sigy*randn(Nk,1);
Dobs=max(Dobs,0);
tobs=t(thesetrain);


figure(1);
plot(t,Dtrue,'-k'); hold on;
scatter(tobs,Dobs,150,'.','MarkerEdgeColor',[0.8 0 0]);
ylabel('SWE [m]','Interpreter','Latex','FontSize',20);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.001, 0.005]);
set(groot, 'defaultAxesTickLabelInterpreter','Latex'); % Latex or tex
set(gca,'TickLabelInterpreter','Latex');
set(gca,'GridLineStyle','-');
set(gca,'FontSize',16);
drawnow;
xlim([datenum(tstart) datenum(tend)]);
datetick('x','keepticks','keeplimits');

%% Bayesian optimization 

% Initial guess
Nstart=2;
c=5.*rand(1,Nstart);
c=c';

% Search grid
cg=0.1:0.01:5;
cg=cg';

sigy=1e-1; % Noise term, in the same units as the output variable y
% Note that this is not the same as sigy for the syn obs, instead
% it corresponds to the so-called "nugget effect"

hpri.sfm=0.5; hpri.sfs=1.4; 
hpri.lsm=0.5; hpri.lss=1.4;

nits=20;
for k=1:nits
    D=DDM(t,T,P,a,b,c,Tr,Ts);
    Dp=D(thesetrain,:);
    err=Dp-Dobs;
    rmse=sqrt(mean(err.^2));
    rmse=rmse';
    
    X=c;
    Nf=size(X,2);

    Xp=cg;
    
    % Just do min/max scaling on inputs
    Xmin=min(Xp,[],1);
    Xmax=max(Xp,[],1);
    Xp=(Xp-Xmin)./(Xmax-Xmin);
    X=(X-Xmin)./(Xmax-Xmin);

    y=rmse;
    ym=mean(rmse);
    ys=std(rmse);
    y=(y-ym)./ys;

    % Training
    sigys=sigy./ys;
    [thetaout,Neff]=ISGP(X,y,sigys,hpri);

    kp.l=thetaout(1:Nf);
    kp.sf=thetaout(end);
    kp.sy=sigys;
    [ypm,yps]=GPP(X,y,Xp,kp);

    % Rescale
    ypred=ypm.*ys+ym;
    yps=yps.*ys;


    nsig=3;

    % Acquisition method
    domin=1;
    % Pick new point as that with the potentially lowest RMSE (including GP uncertainty)
    if domin % Find minimum
        optimist=ypred-nsig*yps; % The most optimistic surface for RMSE at each input
        acquire=optimist==min(optimist);
    else % Maximum variance search
        varis=yps.^2;
        acquire=varis==max(varis);
    end
    cnew=cg(acquire); 
    cnew=min(cnew); % In case multiple minima, just pick one for now
    


    fh=figure(2); clf;
    fh.WindowState = 'maximized';
    ystr='SWE RMSE';
    percplot(cg,[ypred-nsig.*yps ypred+nsig.*yps],[0 0 0.8],0.1,'-','none');
    hold on;
    plot(cg,ypred,'LineWidth',2,'Color',[0 0 0.8 0.4]);
    hold on;
    %plot(c,rmse,'-','LineWidth',3,'Color',[0.5 0.5 0.5]);
    scatter(c,rmse,500,'.',...
        'MarkerEdgeColor',[0.8 0 0],...
        'MarkerFaceAlpha',0.4);
    scatter(cnew,ypred(cg==cnew),2000,'.',...
        'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceAlpha',0.8);
    scatter(cnew,ypred(cg==cnew),500,'.',...
        'MarkerEdgeColor',[0.8 0.8 0],...
        'MarkerFaceAlpha',0.8);
    xlabel('Snowfall factor $c$','Interpreter','Latex','FontSize',20)
    ylabel(ystr,'Interpreter','Latex','FontSize',20);
    title(sprintf('Iteration %d (press any key for next iteration)',k),'Interpreter','Latex',...
        'FontSize',20);
    set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.001, 0.005]);
    set(groot, 'defaultAxesTickLabelInterpreter','Latex'); % Latex or tex
    set(gca,'TickLabelInterpreter','Latex');
    set(gca,'GridLineStyle','-');
    set(gca,'FontSize',16);
    pause(0.1);
    drawnow;
    %}

    c=[c; cnew]; % Update c with the new point for next pass
 
    dowait=waitforbuttonpress();




end


% Pick optimum as location with lowest mean RMSE prediction
optis=ypred==min(ypred);
copt=cg(optis);
fprintf('\n Found copt=%4.2f, true c=%4.2f \n',...
    copt,ctrue)





    





