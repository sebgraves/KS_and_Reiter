function [param,glob] = setup(param,glob,options)

%% Create idiosyncratic productivity matrix and grid
egrid   = [0,1];
egrid0  = egrid;
P       = [1 - param.pi01, param.pi01; param.pi10, 1- param.pi10];

%% Create state space for k
kgrid           = nodeunif(glob.n(1),glob.kmin.^glob.curv,glob.kmax.^glob.curv).^(1/glob.curv);
kgrid0          = kgrid;

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',kgrid,0,glob.spliorder(1)},...
                         {'spli',egrid,0,glob.spliorder(2)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Find actual grid given that fspace has added points for cubic spline
kgrid           = s(s(:,2)==s(1,2),1);
egrid           = s(s(:,1)==s(1,1),2); % This is picking one of the other state and then 
Nk              = size(kgrid,1);
Ne              = size(egrid,1);

%% Compute expectations matrix
Phi             = funbas(fspace,s);
PnK             = kron(P,speye(Nk));    % This is for expected value fn

%% Compute fine frid for histogram
kgridf          = nodeunif(glob.nf(1),glob.kmin.^glob.curv,glob.kmax.^glob.curv).^(1/glob.curv);
Nkf             = size(kgridf,1);
egridf          = egrid;                % As doing discrete, same grid for z in histogram
Nef             = size(egridf,1);
sf              = gridmake(kgridf,egridf);
Nsf             = size(sf,1);
Phif            = funbas(fspace,sf);
PnKf            = kron(P,speye(Nkf));

s_ind           = gridmake(kgridf,egrid);


%% Compute QZ matrix for approximating stationary distribution
Qe              = kron(P,ones(Nkf,1));
%% Create other one time only basis matrices
Phi2            = splibas(egrid0,0,glob.spliorder(2),s(:,2));
Phi2f           = splibas(egrid0,0,glob.spliorder(2),sf(:,2));
%% Declare additional global variables
glob.kgrid0     = kgrid0;
glob.kgrid      = kgrid;
glob.kgridf     = kgridf;
glob.egrid0     = egrid0;
glob.egrid      = egrid;
glob.egridf     = egridf;
glob.P          = P;
glob.Qe         = Qe;
glob.Ne         = Ne;
glob.Nk         = Nk;
glob.Nkf        = Nkf;
glob.sf         = sf;
glob.Nsf        = Nsf;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;
glob.Phi2       = Phi2;
glob.Phi2f      = Phi2f;
glob.Phi        = Phi;
glob.Phif       = Phif;
glob.PnK        = PnK;
glob.PnKf       = PnKf;
glob.s_ind      = s_ind;
