n=10; // Nombre de pages 

// 1- Matrice d'adjacence du graphe des pages 
//    B(i,j) = 1 si il existe un lien i -> j 

Adj=grand(n,n,'bin',1,0.2);show_adj(Adj);

// construction de la matrice de google associée 

[P,Pss,Pprim,d,z,alpha]=google(Adj);
// verification que P est stochastique 

sum(P,'c')

// iteration sur les valeurs 
// P = alpha*Prim + (1-alpha)e*z

// problème découplés x par x 
// La matrice R pour maximiser le page rank des 
// m premieres pages 
m = n/2;
Rm= zeros(n,n);
Rm(1:m,:)=1;

x=1;

// itérations sur les valeurs 

w=ones(n,1);
// w etant fixe on cherche le controle qui maximise Tbar 

controls=get_control(m);
ncont = size(controls,'r');

P1=Pprim;
P1opt= P1;

w=ones(n,1);
while %t then 
  // calcul d'un nouveau controle a w fixé 
  // i.e un nouveau controle par etat x=1:m
  for x=1:m
    // l'etat x 
    cost = - %inf ;
    for k=1:ncont 
      pnux= P1(x,:);
      control=controls(k,:);
      pnux(m+1:$)=control;
      if sum(pnux) < %eps then continue;end
      pnux= pnux/sum(pnux);
      costn = alpha*pnux*w + (alpha*pnux+(1-alpha)*z)*Rm(x,:)';
      if costn > cost then popt = pnux;cost = costn; optcontrol= control;end
    end
    P1opt(x,:)= pnux;
  end
  // calcul d'un nouveau w 
  // w = alpha*P1*w + sum((alpha*P1+ (1-alpha)*ones(n,1)*z)*Rm,'c');
  b = sum((alpha*P1opt+ (1-alpha)*ones(n,1)*z)*Rm,'c');
  wn = inv(eye(n,n)- alpha*P1opt)*b;
  if norm(w-wn) < 1.e-5 then printf("stop \n"); break;end 
  printf("norm(w-wn)=%f\n",norm(w-wn));
  w=wn;
end

// trouver la matrice d'adjacence optimisée 

  













