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
scf(1);
show_adj(Adj);

// Construction de la matrice de transition P
// associée à une matrice d'adjacence.
// Pss: transition d'origine,
// P: matrice de google
// z: vecteur de teleportation
// d: vecteur vaut 1 si le degré vaut zero et 0 sinon

function [P,Pss,Pprim,d,z,alpha]=google(Adj)
    z = grand(1, n, 'unf', 0, 1);
    z = z ./ sum(z);

    d0 = sum(Adj, 'c');
    d1 = double(d0 == 0); // d à remplacer à la fin.
    d2 = find(d1 == 1);
    d3 = zeros(n, 1);
    d3(d2) = %inf;
    d4 = sum(d1);

    Pss = Adj ./ ((d0 + d3) * ones(1, n));

    Pprim = Pss;

    Pprim(d2, :) = ones(d4, 1) * z;
    alpha = 0.8;
    P = alpha * Pprim + (1 - alpha) * (ones(n, 1) * z);
    d = d1;
endfunction

[P,Pss,Pprim,d,z,alpha]=google(Adj);
// verification que P est stochastique

disp(sum(P,'c'));

disp("== Question 3 ==");

w = (alpha * (d - ones(n, 1)) + ones(n, 1))';
x = rand(n, 1)
y1 = P' * x;
y2 = alpha * Pss' * x + z' * (w * x);

disp(clean(y1 - y2));

disp("== Question 4 ==");

// pi^T est un vecteur propre de P^T associé à la valeur propre 1
[R, diagevals]=spec(P');
index = find(abs(diag(diagevals - 1)) < 1.e-10); // indice dans la liste R des vecteurs propres associés à la valeur propre 1.
pi = real(R(:, index))';
pi = pi / sum(pi);

scf(2);
show_adj(Adj_opt,300*int(pi_opt));

disp(pi);
disp(clean(pi * P - pi));

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
disp(pi);
disp(clean(pi * P - pi));

disp("== Question 6 ==");

function [pi]=pi_iterative_sparse(Pss, d, z, alpha)
  p = ones(n, 1);
  while %t
    w = (alpha * (d - ones(n, 1)) + ones(n, 1))';
    pn = alpha * Pss' * p + z' * (w * p);
    if norm(pn - p, %inf) < 10 * %eps then break; end
    p = pn;
  end
  pi = p';
  pi = pi / sum(pi);
endfunction

pi = pi_iterative_sparse(Pss, d, z, alpha);
disp(clean(pi * P - pi));

disp("===== Partie 3 =====");

disp("== Question 7 ==");

p=2;
m = n / 2;

direction_matrix = repmat([0;1], [2**(m-1), 1]);
for k = 1 : m - 1
  val = repmat(cat(1, zeros(2**k,1), ones(2**k,1)), [2**(m-k-1), 1]);
  direction_matrix = cat(2, val, direction_matrix);
end

p_rank = -%inf;

// Valable que pour p = 2
for p1 = 1 : 2**m
  for p2 = 1 : 2**m
    Adj_temp = Adj;
    Adj_temp(1:2,m+1:$) = [direction_matrix(p1, :); direction_matrix(p2, :)];
    [P_temp,Pss_temp,Pprim_temp,d_temp,z_temp,alpha_temp]=google(Adj_temp);
    pi_temp = pi_iterative_sparse(Pss_temp, d_temp, z_temp, alpha_temp);
    p_rank_temp = sum(pi_temp(1:m));
    if p_rank_temp > p_rank
      Adj_opt = Adj_temp;
      pi_opt = pi_temp;
      p_rank = p_rank_temp;
    end
  end
end

scf(3);
show_adj(Adj_opt,300*int(pi_opt));

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
  cerg = Pk * couts;
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

disp(CT);
disp(c);
disp((c - CT) / c);

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
disp(clean(A * w + R - c));
disp("A * c");
disp(clean(A * c));
// Noter que w n’est pas unique, on peut rajouter à w les elts du noyau de A

// Montrons inversement que c doit être egal à Pr*R
// Pr*A est nul
disp("Pr * A");
disp(clean(Pr * A));
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

w = inv(eye(n, n) - alpha * P1) * (sum(P .* Rm, 'c'));

// calcul de c
c = (1-alpha)*z*w;

// (w,c) solution du pb ergodique ?

disp("w + c - (P * w + sum(P .* Rm))");
disp(clean(w + c - (P * w + sum(P .* Rm,'c'))));

// Maintenant on peut utiliser une méthode itérative

disp("== Question 11 ==");

function w=iterative_c(tol)
  b = sum(P .* Rm, 'c');
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

disp("w + c - (P * w + sum(P .* Rm))");
disp(clean(w + c - (P * w + sum(P .* Rm,'c'))));

disp("== Question 12 ==");
