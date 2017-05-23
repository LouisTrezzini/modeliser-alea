disp("===== Partie 2 =====");

disp("== Question 2 ==");

n=10; // Nombre de pages

function show_adj(Adj,diameters)
  [lhs,rhs]=argn(0);
  if rhs < 2 then diameters = 30*ones(1,n);end
  graph = mat_2_graph(sparse(Adj),1,'node-node');
  graph('node_x')=300*cos(2*%pi*(1:n)/(n+1));
  graph('node_y')=300*sin(2*%pi*(1:n)/(n+1));
  graph('node_name')=string([1:n]);
  graph('node_diam')= diameters;
  //graph('node_color')= 1:n;
  //show_graph(graph);
  rep=[1 1 1 1 2 2 2 2 2 2 2 2 2];
  plot_graph(graph,rep);
endfunction

Adj=grand(n,n,'bin',1,0.2);
// show_adj(Adj);

// Construction de la matrice de transition P
// associée à une matrice d'adjacence.
// Pss: transition d'origine,
// P: matrice de google
// z: vecteur de teleportation
// d: vecteur vaut 1 si le degré vaut zero et 0 sinon

function [P,Pss,Pprim,d,z,alpha]=google(Adj) // TODO
    z = grand(1, n, 'unf', 0, 1);
    z = z / sum(z);

    deg = sum(Adj,'c');
    d = double(deg == 0);

    N = deg;
    K = find(N == 0);
    N(K) = %inf;

    Pss = Adj ./ (N * ones(1, n));

    Pprim = Pss;
    Pprim(K, :) = ones(size(K, '*'), 1) * z;

    alpha = 0.8;
    P = alpha * Pprim + (1 - alpha) * ones(n, 1) * z;
endfunction

[P,Pss,Pprim,d,z,alpha]=google(Adj);
// verification que P est stochastique

disp(sum(P,'c'));

disp("== Question 3 == TODO");

dz = (d + (1 - alpha) * (1 - d))';
x = rand(n, 1)
y1 = P' * x;
y2 = alpha * Pss' * x + z' * (dz * x);

disp(y1 - y2);

disp("== Question 4 == FIXME");

// pi^T est un vecteur propre de P^T associé à la valeur propre 1
[R, diagevals]=spec(P');
I = find(abs(diag(diagevals) - 1) < 1.e-7);
pi = real(R(:, I))';
pi = pi / sum(pi);

disp(clean(pi * P - pi));

// xbasc(); show_adj(Adj, int(300 * pi));

disp("== Question 5 ==");

function [pi]=pi_iterative(P)
  p = ones(n, 1);
  while %t
    pn = P' * p;
    if norm(pn - p, %inf) < 10 * %eps then break;end
    p = pn;
  end
  pi = p';
  pi = pi / sum(pi);
endfunction

pi = pi_iterative(P);
disp(clean(pi * P - pi));

disp("== Question 6 == FIXME");

function [pi]=pi_iterative_sparse(Pss, d, z, alpha)
  p = ones(n, 1);
  while %t
    dz = (d + (1 - alpha) * (1 - d))';
    pn = alpha * Pss' * p + z' * (dz * p);
    if norm(pn - p, %inf) < 10 * %eps then break; end
    p = pn;
  end
  pi = p';
  pi = pi / sum(pi);
endfunction

pi = pi_iterative_sparse(Pss, d, z, alpha);
disp(clean(pi * P - pi));

disp("===== Partie 3 =====");

disp("== Question 7 == TODO");

m = n / 2;

// TODO

function controls=get_control(n)
  if n == 1 then
      controls = [0;1];
  else
    A = get_control(n - 1);
    controls = [zeros(size(A,'r'),1),A;ones(size(A,'r'),1),A];
  end
endfunction

controls = get_control(m);
// 2^5 controles
ncont = size(controls,'r');

// On suppose que l'on peut changer les liens sortants des
// deux premieres pages

costopt = -%inf;
kopt = [];

for k = 1:ncont
  for j = 1:ncont
    control1 = controls(k,:);
    control2 = controls(j,:);
    Adjc=Adj;
    Adjc(1:2, m + 1:$) = [control1 ; control2];
    [P,Pss,d,z,alpha] = google(Adjc);
    // pi = pi_iterative_sparse(Pss, d, z, alpha);
    pi = pi_iterative(P);
    cost = sum(pi(1:m));
    if cost > costopt then
        Adjopt=Adjc;
        piopt=pi;
        costopt = cost;
        kopt = [k, j];
    end
  end
end

// La matrice d'adjacence du graphe apr�s optimisation
/*xset('window',1);show_adj(Adjc,300*int(pi));*/

// ---- FIXME

disp("===== Partie 4 =====");

function y=r(x)
  y=x^2
endfunction

n = 4;
P = rand(n, n)
pr = sum(P,'c');
P = P ./ (pr * ones(1, n));

disp(P);

// on suppose ici que les etats de la chaine sont 1:n

disp("== Question 8 ==");

function cerg=ergodique_markov_T(T,P)
  Pk = eye(n,n);
  couts = r((1:n)');
  cerg = couts;
  for i = 1:T
    Pk = P * Pk;
    cerg = cerg + Pk * couts;
  end
  cerg = cerg / T;
endfunction

function [cerg,pi]=ergodique_markov(P)
  pi = pi_iterative(P);
  couts = r((1:n)');
  cerg = pi * couts;
endfunction

// test
T=100000; CT = ergodique_markov_T(T,P);
[c,pi]=ergodique_markov(P);

disp(c - CT);

// Le noyau de P-I est engendré par ones(n,1)
[x0,K]=linsolve(P - eye(n,n), zeros(n,1));

disp("== Question 9 ==");

// le projecteur spectral sur Espace propre associé a 1
Pr = ones(n,1) * pi; // [pi;pi;pi;....]
A = P - eye(n,n);    // A -Id
S = Pr - inv(Pr - A);
// vérifier que S*Pr et Pr*S sont nuls
disp("S * Pr");
disp(clean(S*Pr));
disp("Pr * S");
disp(clean(Pr*S));
// A*w  + R - c= 0
// A*c         = 0
R = r((1:n)');
// vérifions que w=-S*R et c=Pr*R sont solution du systeme linaire
w = -S * R;
c = Pr * R;
disp("A * w + R - c");
disp(A * w + R - c);
disp("A * c");
disp(A * c);
// Noter que w n’est pas unique, on peut rajouter à w les elts du noyau de A

// Montrons inversement que c doit être egal à Pr*R
// Pr*A est nul
disp("Pr * A");
disp(Pr * A);
// on doit donc avoir
// Pr*R - Pr*c = 0 et A*c =0
// en sommant
// Pr*R = (Pr-A)*c
// c = (Pr-A)^-1 *Pr*R
// c = (Pr-S)*Pr*R = Pr*Pr*R -S*Pr*R = Pr*R
// car Pr est un projecteur Pr^2 = Pr et S*Pr = 0
disp(clean(Pr^2-Pr));
disp(clean(S*Pr));
// conclusion c doit valoir Pr*R
// on le vérifie avec linsolve
[x0,K]=linsolve([A,-eye(n,n);zeros(n,n),A],[R;zeros(n,1)]);
// on vérifie bien que e = Pr*R

P1 = rand(n,n);
pr = sum(P1,'c');
P1 = P1 ./ (pr * ones(1, n));

z = grand(1, n, 'unf', 0, 1);
z = z / sum(z);

alpha = 0.8;

P = alpha * P1 + (1 - alpha) * ones(n, 1) * z;

// les couts Rm(i,j)
Rm = grand(n,n,'unf',0,1);

disp("== Question 10 ==");

// On le vérifie numeriquement
// trouver la solution de
// w = alpha*P1*w + sum(P.*Rm,'c')

w = inv(eye(n, n) - alpha * P1) * (sum(P .* Rm, 'c'))

// calcul de c
c = (1-alpha)*z*w

// (w,c) solution du pb ergodique ?

disp("w + c - (P*w + sum(P.*R,'c'))");
disp(w + c - (P*w + sum(P.*R,'c')));

// Maintenant on peut utiliser une méthode itérative

disp("== Question 11 ==");

function w=iterative_c(tol)
  b = sum(P .* R, 'c');
  w = ones(n, 1);
  while %t
    wn = alpha * P1 * w + b
    if norm(wn - w) < tol then break;end
    w = wn;
  end
endfunction

w = iterative_c(10 * %eps);
// calcul de c
c = (1 - alpha) * z * w

// (w,c) solution du pb ergodique ?

disp("w + c - (P*w + sum(P.*R,'c'))");
disp(w + c - (P*w + sum(P.*R,'c')));
