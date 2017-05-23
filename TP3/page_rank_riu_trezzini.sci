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

show_adj(Adj);  
 
// Construction de la matrice de transition P  
// associée à une matrice d'adjacence.  
// Pss: transition d'origine,  
// P: matrice de google  
// z: vecteur de teleportation  
// d: vecteur vaut 1 si le degré vaut zero et 0 sinon  
 
function [P,Pss,Pprim,d,z,alpha, d0, d3]=google(Adj)
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
 
[P,Pss,Pprim,d,z,alpha, d, f]=google(Adj);  
// verification que P est stochastique  
sum(P,'c')

x = rand(n, 1)  
y1 = P' * x;  
y2 = alpha * Pss' * x + alpha * z' * ((d - ones(n, 1))' * x) + z' * sum(x);
y1-y2
