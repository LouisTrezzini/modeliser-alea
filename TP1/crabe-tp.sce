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

histplot(y, crabe);

crabeNormal = normale(y, mean(crabe), stdev(crabe)^2);
plot(y, crabeNormal, 'r');

alpha = test_chi2(crabe, normale(crabe, mean(crabe), stdev(crabe)^2))

// Données
pi0=[1; 3]/2/2;
pi=pi0;
mu=[.57; .67];
s2=[1 ;1]/10000;

rho=ones(2,1000);

// Algorithme EM pour les crabes
//------------------------------

N=1000;
iterMax = 1000;
R=zeros(iterMax+1, 5);
R(1, :) = [mu(1); mu(2); pi(1); s2(1); s2(2)]';

for iter=1:iterMax
    rho_p = zeros(2, N);
    rho_p(1, :) = (pi(1) * normale(crabe, mu(1), s2(1)) ./ (pi(1) * normale(crabe, mu(1), s2(1)) + pi(2) * normale(crabe, mu(2), s2(2))))';
    rho_p(2, :) = (pi(2) * normale(crabe, mu(2), s2(2)) ./ (pi(1) * normale(crabe, mu(1), s2(1)) + pi(2) * normale(crabe, mu(2), s2(2))))';

    // Etape E
    pi_star = mean(rho_p, 'c');

    // A0 = pi_star' * log(pi);
    // A1 = rho_p(1, :) *  log(normale(crabe, mu(1), s2(1));
    // A1 = rho_p(2, :) *  log(normale(crabe, mu(2), s2(2));
    // Q = N * A0 + A1 + A2;

    // Etape M
    mu_star = rho_p * crabe ./ sum(rho_p, 'c');
    crabe_mean = zeros(N, 2);
    crabe_mean(:, 1) = crabe - mu_star(1);
    crabe_mean(:, 2) = crabe - mu_star(2);
    s2_star = sum(rho_p .* (crabe_mean .^ 2)', 'c') ./ sum(rho_p, 'c');

    R(iter + 1, :) = [mu_star(1); mu_star(2); pi_star(1); s2_star(1); s2_star(2)]';
    mu = mu_star;
    pi = pi_star;
    s2 = s2_star;
end

// Affichages
f = scf();
histplot(y, crabe);
crabeNormal1 = pi(1) * normale(y, mu(1), s2(1));
crabeNormal2 = pi(2) * normale(y, mu(2), s2(2));
plot(y, crabeNormal1, 'r');
plot(y, crabeNormal2, 'g');
plot(y, crabeNormal1 + crabeNormal2, 'c');
