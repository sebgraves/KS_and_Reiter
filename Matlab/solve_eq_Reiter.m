function eq_Reiter = solve_eq_Reiter(eq_SS,param,glob,options);

% Create X_SS vector
X_SS = [eq_SS.c;eq_SS.L;0;eq_SS.c];

%% Calculate F1,F2,F3,F4
[F1,F2,F3,F4] = take_derivatives(X_SS,param,glob,options);

%% Convert for gensys and solve linear rational expectations problem
G0  = -F1;
G1  = F2;
C   = zeros(size(X_SS,1),1);
Psi = F3;
Pi  = F4;

[A_gen, ~, B_gen, ~, ~, ~, ~, eu] = gensys(G0, G1, C, Psi, Pi);

eq_Reiter.A     = A_gen;
eq_Reiter.B     = B_gen;
eq_Reiter.eu    = eu;
eq_Reiter.X_SS  = X_SS;