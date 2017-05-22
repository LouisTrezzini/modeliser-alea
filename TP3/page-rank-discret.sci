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

Adj=grand(n,n,'bin',1,0.2);show_adj(Adj);

// Construction de la matrice de transition P
// associée à une matrice d'adjacence.
// Pss: transition d'origine,
// P: matrice de google
// z: vecteur de teleportation
// d: vecteur vaut 1 si le degré vaut zero et 0 sinon

function [P,Pss,Pprim,d,z,alpha]=google(Adj)
  deg = sum(Adj,'c');
  d = double( deg ==0)
  z=grand(1,n,'unf',0,1);
  z=z/sum(z);
  N=deg;
  K=find(N==0);
  N(K) = %inf;
  Pss= Adj ./ (N*ones(1,n));
  Pprim = Pss;
  Pprim(K,:) = ones(size(K,'*'),1)*z;
  alpha = 0.8;
  P= alpha*Pprim + (1-alpha)*ones(n,1)*z;
endfunction

[P,Pss,Pprim,d,z,alpha]=google(Adj);
// verification que P est stochastique

sum(P,'c')

dz = ( d + (1-alpha)*(1-d))';
x=rand(n,1)
y1=P'*x;
y2=alpha*Pss'*x + z'*(dz*x);

y1-y2

[X,diagevals]=spec(P');
I=find( abs(diag(diagevals)- 1) < 1.e-7);
pi= real(X(:,I))';
pi= pi/sum(pi);

clean(pi*P - pi)

xbasc();show_adj(Adj,int(300*pi));

function [pi]=pi_iterative()
  p=ones(n,1);
  while %t
    pn = P'*p;
    if norm(pn-p,%inf) < 10*%eps then break;end
    p = pn;
  end
  pi = p';
  pi= pi/sum(pi);
endfunction

pi=pi_iterative();
clean(pi*P - pi)

function [pi]=pi_iterative_sparse()
  p=ones(n,1);
  while %t
    pn = alpha*Pss'*p + z'*(dz*p);
    if norm(pn-p,%inf) < 10*%eps then break;end
    p = pn;
  end
  pi = p';
  pi= pi/sum(pi);
endfunction

pi=pi_iterative_sparse();
clean(pi*P - pi)

// maximization du << page rank >>
m=n/2;

// les liens des m premiers noeuds peuvent etre chang�s
// et les m-premiers liens sont interdits i.e dans la
// matrice d'adjacence on change Adj(1:m,m+1:$)
// get_controls renvoit une matrix ncontxn
// chaque ligne repr�sentant un choix possible pour n-liens

function controls=get_control(n)
  if n==1 then controls=[0;1];
  else
    A=get_control(n-1);
    controls = [zeros(size(A,'r'),1),A;ones(size(A,'r'),1),A];
  end
endfunction

controls = get_control(m);
// 2^5 controles
ncont = size(controls,'r');

// On suppose que l'on peut changer les liens sortants des
// deux premieres pages

costopt = -%inf;
kopt =[];

for k=1:ncont
  for j=1:ncont
    control1=controls(k,:);
    control2=controls(j,:);
    Adjc=Adj; Adjc(1:2,m+1:$)=[control1;control2];
    [P,Pss,d,z]=google(Adjc);
    // calcul de pi
    p=ones(n,1);
    while %t
      pn = P'*p;
      if norm(pn-p,%inf) < 10*%eps then break;end
      p = pn;
    end
    pi = p';
    pi = pi/sum(pi);
    cost = sum(pi(1:m));
    if cost > costopt then Adjopt=Adjc;piopt=pi; costopt = cost; kopt = [k,j];end
  end
end

// La matrice d'adjacence du graphe apr�s optimisation
xset('window',1);show_adj(Adjc,300*int(pi));
