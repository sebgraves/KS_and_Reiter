function [param,glob] = setup(param,glob,options)

%% Create idiosyncratic productivity matrix and grid
egrid   = [0,1];
egrid0  = egrid;
P_e     = [1 - param.pi01, param.pi01; param.pi10, 1- param.pi10];

%% Create state space for k
kgrid           = nodeunif(glob.n(1),glob.kmin.^glob.curv,glob.kmax.^glob.curv).^(1/glob.curv); 
kgrid0          = kgrid;

%% Create state space for K
Kgrid           = nodeunif(glob.n(3),glob.Kmin,glob.Kmax);
Kgrid0          = Kgrid;

%% Create state space for A
Agrid           = nodeunif(glob.n(4),glob.Amin,glob.Amax);
Agrid0          = Agrid;

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',kgrid,0,glob.spliorder(1)},...
                         {'spli',egrid,0,glob.spliorder(2)},...
                         {'spli',Kgrid,0,glob.spliorder(3)},...
                         {'spli',Agrid,0,glob.spliorder(4)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);
Phi             = funbas(fspace,s);

%% Find actual grid given that fspace has changed points for splines
kgrid           = unique(s(:,1));
egrid           = unique(s(:,2));
Kgrid           = unique(s(:,3));
Agrid           = unique(s(:,4));

Nk              = size(kgrid,1);
Ne              = size(egrid,1);
NK              = size(Kgrid,1);
NA              = size(Agrid,1);

%% Compute quadrature nodes and weights
[eta,f]         = qnwnorm(glob.Na,0,param.sigA^2);
%% Compute expectations matrix
Emat            = kron(kron(kron(speye(NA),speye(NK)),P_e),speye(Nk))*kron(speye(Ns),f');

%% Compute fine grid for solving accurately
kgridf          = nodeunif(glob.nf(1),glob.kmin.^glob.curv,glob.kmax.^glob.curv).^(1/glob.curv); 
Nkf             = size(kgridf,1);
sf              = gridmake(kgridf,egrid,Kgrid,Agrid);
Nsf             = size(sf,1);
%Phif            = funbas(fspace,sf);
%Ematf           = kron(kron(kron(speye(NA),speye(NK)),P_e),speye(Nkf))*kron(speye(Nsf),f');

%% Compute other one-time matrices
sa              = [kron(s(:,1),ones(glob.Na,1)),kron(s(:,2),ones(glob.Na,1)),kron(s(:,3),ones(glob.Na,1)),g_fun(kron(s(:,4),ones(glob.Na,1)),kron(ones(Ns,1),eta),param,glob)];
Phia            = funbas(fspace,sa);
%saf             = [kron(sf(:,1),ones(glob.Na,1)),kron(sf(:,2),ones(glob.Na,1)),kron(sf(:,3),ones(glob.Na,1)),g_fun(kron(sf(:,4),ones(glob.Na,1)),kron(ones(Nsf,1),eta),param,glob)];
%Phiaf           = funbas(fspace,saf);

Phie            = splibas(egrid,0,glob.spliorder(2),s(:,2));
%Phief           = splibas(egrid,0,glob.spliorder(2),sf(:,2));
PhiA            = splibas(Agrid0,0,glob.spliorder(4),s(:,4));
%PhiAf           = splibas(Agrid0,0,glob.spliorder(4),sf(:,4));

fspacek         = fundef({'spli',kgridf,0,1});

%% Compute fine grid on individual states for tracking distribution
s_ind           = gridmake(kgridf,egrid);
Ns_ind          = size(s_ind,1);

%% Compute Qe matrix for updating distribution
Qe              = kron(P_e,ones(Nkf,1));

%% Declare additional global variables
glob.P_e    = P_e;
glob.Ns     = Ns;
glob.Phi    = Phi;
glob.egrid  = egrid;
glob.kgrid  = kgrid;
glob.Kgrid  = Kgrid;
glob.Agrid  = Agrid;
glob.Nk     = Nk;
glob.Ne     = Ne;
glob.NK     = NK;
glob.NA     = NA;
glob.eta    = eta;
glob.f      = f;
glob.Emat   = Emat;
glob.kgridf = kgridf;
glob.Nkf    = Nkf;
glob.Nsf    = Nsf;
%glob.Phif   = Phif;
%glob.Ematf  = Ematf;
glob.s_ind  = s_ind;
glob.Ns_ind = Ns_ind;
glob.Qe     = Qe;
glob.Phia   = Phia;
%glob.Phiaf  = Phiaf;
glob.fspace = fspace;
glob.s      = s;
glob.sf     = sf;
glob.Kgrid0 = Kgrid0;
glob.Agrid0 = Agrid0;
glob.kgrid0 = kgrid0;
glob.Phie   = Phie;
glob.PhiA   = PhiA;
%glob.Phief  = Phief;
%glob.PhiAf  = PhiAf;
glob.fspacek = fspacek;
