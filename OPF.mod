#limpiar memoria 
reset;
#iniciar en modo modelo
model;
#declarra conjuntos
set OB;
				
set OL within 1..1000 cross OB cross OB;

#declarar parametros 
#parametros del sistema
param Sbase;

#parametros de las barras
param Nombre{OB} symbolic;
param Tb{OB};
param V0{OB};
param th0{OB};
param Pd{OB};
param Qd{OB};
param Pg0{OB};
param Qg0{OB};
param Vg{OB};
param Qgmax{OB};
param Qgmin{OB};
param gsh{OB};
param bsh{OB};
param Vmax{OB};
param Vmin{OB};
param Pgmax{OB};
param Pgmin{OB};

#parametros de los ramales
param Tr{OL};
param r{OL};
param x{OL};
param bshl{OL};
param a{OL};
param amax{OL};
param amin{OL};
param fi{OL};
param Smax{OL};
param g{OL};
param b{OL};

#Declarar Variables
var V{OB};                  #Voltaje en la barra i
var th{OB};                 #Anguki en la barra i con respecto a
var Pg{OB};                 #Potencia activa generada en la barra
var Qg{OB};                 #potencia reactiva generada en la barra
var Pde{OL};                #flujo de potencia activa desde la barra
var Qde{OL};                #flujo de potencia reactiva desde la
var Ppa{OL};                #flujo de potencia activa desde la barra
var Qpa{OL};                #flujo de potencia reactiva desde la

#funcion objetivo 
minimize generacion_slack:                           #el objetivo es minimizar                      
        sum{i in OB:Tb[i]==3}(Pg[i]);


subject to balance_activa{i in OB}:	
Pg[i]+gsh[i]*V[i]^2-Pd[i]-sum{(k,i,j) in OL}(Pde[k,i,j])
-sum{(k,j,i) in OL}(Ppa[k,j,i])=0;

subject to balance_reactiva{i in OB}:
Qg[i]+bsh[i]*V[i]^2-Qd[i]-sum{(k,i,j) in OL}(Qde[k,i,j])
-sum{(k,j,i) in OL}(Qpa[k,j,i])=0;

#flujo de potencia activa desde la barra i hacia la barra j
subject to flujo_activa_de {(k,i,j) in OL}:
Pde[k,i,j]=g[k,i,j]*a[k,i,j]^2*V[i]^2-
a[k,i,j]*V[i]*V[j]*g[k,i,j]*cos(th[i]-th[j]+fi[k,i,j])
-a[k,i,j]*V[i]*V[j]*b[k,i,j]*sin(th[i]-th[j]+fi[k,i,j]);

#flujo de potencia reactiva desde la barra i hacia la j
subject to Flujo_reactivo_de{(k,i,j) in OL}:
Qde[k,i,j]=-(b[k,i,j]+bshl[k,i,j])*a[k,i,j]^2*V[i]^2
-a[k,i,j]*V[i]*V[j]*g[k,i,j]*sin(th[i]-th[j]+fi[k,i,j])
+a[k,i,j]*V[i]*V[j]*b[k,i,j]*cos(th[i]-th[j]+fi[k,i,j]);

#flujo de potencia activa desde la barra j hacia la barra i
subject to flujo_activa_pa {(k,i,j) in OL}:
Ppa[k,i,j]=g[k,i,j]*V[j]^2-
a[k,i,j]*V[i]*V[j]*g[k,i,j]*cos(th[i]-th[j]+fi[k,i,j])+
a[k,i,j]*V[i]*V[j]*b[k,i,j]*sin(th[i]-th[j]+fi[k,i,j]);

#flujo de potencia reactiva desde la barra j hacia la i
subject to flujo_reactivo_de{(k,i,j) in OL}:
Qpa[k,i,j]=-(b[k,i,j]+bshl[k,i,j])*V[j]^2
+a[k,i,j]*V[i]*V[j]*g[k,i,j]*sin(th[i]-th[j]+fi[k,i,j])
+a[k,i,j]*V[i]*V[j]*b[k,i,j]*cos(th[i]-th[j]+fi[k,i,j]);

#Limites potencias generadas
subject to limite_Pgenerada {i in OB}:
     Pgmin[i] <= Pg[i] <= Pgmax[i];

subject to limite_Qgenerada {i in OB}:
     Qgmin[i] <= Qg[i] <= Qgmax[i]; 

#Limites de voltaje 
subject to limite_voltaje {i in OB}:
     Vmin[i] <= V[i] <= Vmax[i];	
	 
#limites de potencia en los ramales
subject to limite_ramales_de{(k,i,j) in OL}:
             Pde[k,i,j]^2 + Qde[k,i,j]^2<=Smax[k,i,j]^2;

subject to limite_ramales_pa{(k,i,j) in OL}:
             Ppa[k,i,j]^2 + Qpa[k,i,j]^2 <= Smax[k,i,j]^2;		 
/*
#Restricciones adicionales:
#perdidas de activa y reactiva
subject to perdida_activa_Ppa {(k, i, j) in OL}:
    Ppa[k, i, j] = g[k, i, j] * (V[j]^2) -
    a[k, i, j] * V[i] * V[j] * g[k, i, j] * cos((th[i] - th[j]) + fi[k, i, j]) +
    a[k, i, j] * V[i] * V[j] * b[k, i, j] * sin((th[i] - th[j]) + fi[k, i, j]);

subject to perdida_reactiva_Qpa {(k, i, j) in OL}:
    Qpa[k, i, j] = -(b[k, i, j] + bshl[k, i, j]) * (V[j]^2) +
    a[k, i, j] * V[i] * V[j] * g[k, i, j] * sin((th[i] - th[j]) + fi[k, i, j]) +
    a[k, i, j] * V[i] * V[j] * b[k, i, j] * cos((th[i] - th[j]) + fi[k, i, j]);

#limites de voltaje 
subject to limite_voltaje{i in OB}:
            Vmin[i]<= V[i]<= Vmax[i];

#limites potencia generadas
subject to limite_potencia_generada{i in OB};
            Qmin[i] <= Qg[i] <= Qmax[i];
*/			
		 

/*
#desfase angular
             amin[k,i,j] <= th[i]-th[j] <= amax[k,i,j];
*/