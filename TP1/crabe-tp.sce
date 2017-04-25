clear

// Loi normale
function [x]=normale(y,m,s2)
  x=%e^(-(y-m).^2/2/s2)/sqrt(2*%pi*s2)
endfunction;

function [proba]=test_chi2(N,p0)
  n=sum(N);// taille de l'echantillon observe
  // calcul de zeta_ n
  zeta_n=n*sum(((N/n-p0).^2)./p0);
  // nombre de degres de liberte  (= nombre de classes dans N-1)
  d= length(N)-1;
  // on calcule la proba pour un chi 2 à d-1 degres d'etre superieur a zeta
  [p,q]=cdfchi("PQ",zeta_n,d);
  proba=q;
endfunction;

// Ouvrir le fichier de données (nombre de crabes par intervalle)
x=fscanfMat('crabe.txt');
x=x;

// intervalles
y=.580+.002+.004*[0:28];
yM=y+.002;
ym=y-.002;
Max=25;

crabe = zeros(sum(x));

xSum = cumsum(x);
for index= 1:29 - 1
    crabe(xSum(index) : xSum(index + 1)) = y(index);
end

// Dessiner la loi normale correspondante
X = linspace(0.58, 0.7);
histplot(y, crabe);

crabeNormal = normale(X, mean(crabe), stdev(crabe)^2);
plot(X, crabeNormal, 'r');

alpha = test_chi2(crabe, normale(crabe, mean(crabe), stdev(crabe)^2))

// Données
Npop = 3;

pi0=[1; 3; 2]/2/2;
pi=pi0;
mu=[.57; .67; 0.6];
s2=[1; 1; 1]/10000;

rho=ones(Npop, 1000);

// Algorithme EM pour les crabes
//------------------------------

N=1000;
iterMax = 1000;
R=zeros(iterMax + 1, 3 * Npop);
R(1, :) = [mu; pi; s2]';

for iter=1:iterMax
    rho_p = zeros(Npop, N);

    for k = 1:Npop
        rho_p(k, :) = (pi(k) * normale(crabe, mu(k), s2(k)))';
    end

    norme = sum(rho_p, 'r');
    for k = 1:Npop
        rho_p(k, :) = rho_p(k, :) ./ norme;
    end

    // Etape E
    pi_star = mean(rho_p, 'c');

    // A0 = pi_star' * log(pi);
    // A1 = rho_p(1, :) *  log(normale(crabe, mu(1), s2(1));
    // A1 = rho_p(2, :) *  log(normale(crabe, mu(2), s2(2));
    // Q = N * A0 + A1 + A2;

    // Etape M
    mu_star = rho_p * crabe ./ sum(rho_p, 'c');
    crabe_mean = zeros(N, Npop);
    for k = 1:Npop
        crabe_mean(:, k) = crabe - mu_star(k);
    end
    s2_star = sum(rho_p .* (crabe_mean .^ 2)', 'c') ./ sum(rho_p, 'c');

    R(iter + 1, :) = [mu_star; pi_star; s2_star]';
    mu = mu_star;
    pi = pi_star;
    s2 = s2_star;
end

// Affichages
f = scf();
histplot(y, crabe);
crabe_normal = zeros(Npop, size(X, 2));
for k = 1:Npop
    crabe_normal(k, :) = pi(k) * normale(X, mu(k), s2(k));
end
plot(X, [sum(crabe_normal, 'r'); crabe_normal]');
