// Dans cette partie on vérifie des théorèmes ergodiques 
// On se donne une fonction cout 

function y=r(x) 
  y=x^2
endfunction 

// une matrice stochastique 

n=4;
P=rand(n,n)
pr=sum(P,'c');
P = P ./ (pr*ones(1,n));

// On vérifie ici le théorème ergodique 
// 1/T sum_0^{T-1} E(r(X_t)) -> pi*r
// et la limite ne dépends pas du point initial 
// on suppose ici que les etats de la chaine sont 1:n

function cerg=ergodique_markov_T(T,P)
  cerg=0;
  R=r((1:n)'); 
  Pk=eye(n,n);
  cerg= R;
  for i=1:T
    Pk=P*Pk;
    cerg = cerg + Pk*R;
  end
  cerg=cerg/T;
endfunction 

// calcule de pi la mesure invariante 
function [cerg,pi]=ergodique_markov(P)
  [X,diagevals]=spec(P')
  I=find( abs(diag(diagevals)- 1) < 1.e-7);
  pi=X(:,I)';
  pi=pi/sum(pi);
  R=r((1:n)'); 
  cerg= pi*R;
endfunction 

// test 
T=100000; CT=ergodique_markov_T(T,P);
[c,pi]=ergodique_markov(P);

c-CT

// On veut vérifier maintenant que c est aussi solution du 
// système linéaire 
// (P-I)*w + R - c = 0
// (P-I)*c = 0 
// 
// Quelques vérifications 
// On verifie que (P-I)*c = 0 veut dire que c = c1*ones(n,1);
// Le noyau de P-I est engendré par ones(n,1) 
[x0,K]=linsolve(P- eye(n,n),zeros(n,1));

// le projecteur spectral sur Espace propre associé a 1
Pr = ones(n,1)*pi; // [pi;pi;pi;....] 
A = P-eye(n,n);    // A -Id 
S = Pr - inv(Pr-A)
// vérifier que S*Pr et Pr*S sont nuls 
clean(S*Pr)
clean(Pr*S)
// A*w  + R - c= 0
// A*c         = 0
R=r((1:n)'); 
w= -S*R;
c= Pr*R;
// vérifions qu'ils sont solution du systeme linaire 
A*w + R -c 
A*c
// On peut rajouter à w les elts du noyau de A 

// Montrons inversement que c doit être Pr*R 
// Pr*A est nul 
Pr*A 
// on doit donc avoir 
// Pr*R - Pr*c = 0 et A*c =0 
// en sommant 
// Pr*R = (Pr-A)*c
// c = (Pr-A)^-1 *Pr*R 
// c = (Pr-S)*Pr*R = Pr*Pr*R -S*Pr*R = Pr*R 
// car Pr est un projecteur Pr^2 = Pr et S*Pr = 0
clean(Pr^2-Pr)
clean(S*Pr)
// conclusion c doit valoir Pr*R 
// on le vérifie avec linsolve 
[x0,K]=linsolve([A,-eye(n,n);zeros(n,n),A],[R;zeros(n,1)]);
// on vérifie bien que e = Pr*R 

// On cherche donc ici à résoudre 
// w + c*ones(n,1) = Pw + R;
// On se place dans le cas ou R(i) = sum_j P(i,j)*R(i,j)
//
// w_i + c = sum_j P(i,j) w(j) + P(i,j)*R(i,j) 
// 
// on suppose que P = alpha*P1 + (1-alpha)*e*z

P1=rand(n,n)
pr=sum(P1,'c');
P1 = P1 ./ (pr*ones(1,n));

z=grand(1,n,'unf',0,1);
z=z/sum(z);

alpha = 0.8;
 
P = alpha*P1 + (1-alpha)*ones(n,1)*z;

// les couts R(i,j)
R = grand(n,n,'unf',0,1);

// On cherche un point fixe de Tbar 
// w = alpha*P1*w + P.*R

w = inv(eye(n,n) - alpha*P1)*(sum(P.*R,'c'))

// calcul de c 
c = (1-alpha)*z*w 

// (w,c) solution du pb ergodique ? 

w + c - (P*w + sum(P.*R,'c'))

// Maintenant on peut utiliser une méthode itérative 

function w=iterative_c(tol)
  R1=sum(P.*R,'c');
  w=ones(n,1);
  while %t 
    wn = alpha*P1*w + R1
    if norm(wn-w) < tol then break;end
    w=wn;
  end
endfunction

w=iterative_c(10*%eps);
// calcul de c 
c = (1-alpha)*z*w 

// (w,c) solution du pb ergodique ? 

w + c - (P*w + sum(P.*R,'c'))
  
    


























