function [D,fSCA]=DDM(t,T,P,a,b,c,Tr,Ts)
%% DDM: A simple degree day snowmelt model
% The first dimension of the forcing is assumed to be space, while the second dimension
% is time. 
% Note that this model can also be run in ensemble mode by combining the
% ensemble and spatial dimensions in the first dimension using ensemble chunks.
% That is, with Ns spatial grid cells and Ne ensemble members the first
% dimension is Ng=Ne x Ns and the first Ne elements are the ensemble members
% for grid cell 1, while elements (Ne+1):2*Ne are the corresponding 
% ensemble members for grid cell 2 and so on. 

% Parameters:
% a  = Degree day factor lit. values: 2.5 to 11.6 mm/d/K. 
% b  = Air temperature bias 
% c  = Snowfall scaling

% Hyperparameters:
% Tr = Pure rainfall temperature threshold (Default=2 deg C)
% Ts = Pure snowfall temperature threshold (Default=-2 deg C, Ts<=Tr)
Tmelt=0; % Temperature threshold for snowmelt (degrees C). Uncertainty in this is baked into bT

% Specify ensemble and spatial chunks separately, then merge in the DDM sim
% e.g. input NNc as Ngc x 4 [Ngc : Number of spatial cells per chunk]
% This could be done by broadcasting so that e.g. Tj grows from
% Ng x 1 to Ng x Ne after the b subtraction. Transposes can then help
% track the order of ensemble and spatial dimensions. Only need
% to store the full spatio-temporal state ensemble for D (SWE) and fSCA, 
% not the forcing
%P=sum(Pin(NNc,:))

% Pinj=double(fc.P(:,j)).fc.P_sf*;

Nt=numel(t); % Number of time steps.
Ne=size(c,1);


D=zeros(Nt,Ne);
Dj=zeros(Ne,1);

for j=1:Nt
    Pj=P(j);
    Tj=T(j);
    Tj=Tj-b; % Apply bias correction to temperature.
    frj=(Tj-Ts)./(Tr-Ts); % Rain fraction
    frj=max(frj,0); 
    frj=min(frj,1);
    fsj=1-frj; % Snow fraction
    Sj=c.*fsj.*Pj; % Snowfall 
    fddj=max(Tj-Tmelt,0); % Freezing degree day
    ablationj=a.*fddj;
    netaccj=Sj-ablationj;
    Dj=max(Dj+netaccj,0);
    D(j,:)=Dj;
end



end
