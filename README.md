# Memòria

### Treball previ

Ens diuen que hem d'utilitzar el [mètode de diferències finites](https://en.wikipedia.org/wiki/Finite_difference_method) Per resoldre:

 _Equació de Poisson_


$$
\ u_{xx}(x,y)+u_{yy}(x,y)=f(x,y) ; \forall (x,y)\in D\equiv [a,b]\times[c,d]\in \mathbb{R}
$$
El que hem de fer primer de tot és tractar el domini $D$, per tal de convertir-lo en un domini discretitzat, així doncs considerem la graella de punts $(x_i,y_j)$ resultants de dividir $[a,b]$ en $n $ subintervals de mida $hx$ i $[c,d]$ en $m$ de mida $hy$. D'aquesta manera trobarem només els $u_{i,j}=u(x_i,y_j)$, dels punts de la graella.

Fixem-nos amb quines són les variables i quines les incògnites del problema: coneixem la $f(x,y)$, l'interval $D$ i $u(x,y)$ a $\partial D$ i volem trobar tots els $u_{i,j}$ de la graella.

Tenim per una banda l'equació de Poisson que ens relaciona $u_{xx}, u_{yy}$ amb $f$ i per altra banda les diferències centrals finites que ens relacionen numèricament $f''$ amb $f$.

 _Diferències centrals finites_
$$
p''(x_0)\approx\frac{p(x_0+h)-2p(x_0)+p(x_0-h)}{h^2}
$$


Si apliquem això a les columnes de la nostra graella tindre una aproximació a $u_{xx}$ en el punt $(x_i,y_j)$, en funció de $u_{i-1,j}$,$ u_{i,j}$ i $ u_{i+1,j}$, i si ho fem a les columnes tindrem  $u_{yy}$ en funció de $u_{i,j-1}$,$ u_{i,j}$ i $ u_{i,j+1}$.
$$
u_{xx}(x_i,y_j)\approx\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h^2}
\\
u_{yy}(x_i,y_j)\approx\frac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{h^2}
$$
i ara,
$$
u_{xx}(x,y)+u_{yy}(x,y)=f(x,y) \Rightarrow \frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h^2}+\frac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{h^2}\approx f(x_i, y_j)
\\
$$
Per tant tenim la següent equació de 5 incògnites, amb $h$ i $f$ coneguts. 
$$
u_{i+1,j}+u_{i-1,j}-4u_{i,j}+u_{i,j+1}+u_{i,j-1}\approx h^2\times f(x_i, y_j)
$$
Ara si suposem que $m$ i $n$ són més grans de 5, tindrem un sistema d'equacions d'on podrem resoldre numèricament $u_{i,j} $  $\forall (i,j)$; $0<i<n$, $0<j<m$, amb $f(x_i,y_i)$ coneguda. Volem expressar aquest sistema d'equacions com un sistema de la forma, $A\cdot u=b$. Per fer-ho hem de posar els valors de la graella $(m-1)\times (n-1)$ que representa $D$ en format de vector. De manera que $u$ és:

$$
\qquad u = 
\begin{pmatrix}
u_{1,1} &
u_{1,2} &
\dots &
u_{1,m-1} &
u_{2,1} &
u_{2,2} &
\dots &
u_{2,m-1} &
\dots &
u_{n-1,1} &
u_{n-1,2} &
\dots &
u_{n-1,m-1} 
\end{pmatrix}^\top
$$
El vector $b$ també de mida $(m-1)\cdot(n-1)$ que l'expressem en forma de matriu $(m-1)\times(n-1)$
$$
b = 
\begin{pmatrix}
f_{1,1}-g_{0,1}-g_{1,0} & f_{1,2}-g_{0,2} & \dots & f_{1,m-2}-g_{0,m-2} & f_{1,m-1}-g_{0,n-1}-g_{1,m} \\
f_{2,1}-g_{2,0} & f_{1,2} & \dots & f_{2,m-2} & f_{2,m-1}-g_{2,m} \\
\vdots  &\vdots &\ddots & \vdots & \vdots \\
f_{n-2,1}-g_{n-2,0} & f_{n-2,2} & \dots & f_{n-2,m-2} & f_{n-2,m-1}-g_{n-2,m} \\
f_{n-1,1}-g_{n,1}-g_{n-1,0} & f_{n-1,2}-g_{n,2} & \dots & f_{n-1,m-2}-g_{n,m-2} & f_{n-1,m-1}-g_{n,m-1}-g_{n-1,m} \\
\end{pmatrix}
$$

Veiem que els extrems tenim una incògnita resolta,ja que sabem el valor de la frontera $g$ en $\partial D$. També el convertim en vector:
$$
b = 
\begin{pmatrix}
f_{1,1}-g_{0,1}-g_{1,0} &
f_{1,2}-g_{0,2} &
\dots &
f_{2,1}-g_{2,0} &
f_{1,2} &
\dots &
f_{n-1,m-1}-g_{n,m-1}-g_{n-1,m} 
\end{pmatrix}^\top
$$
I finalment $A$ una matriu de $$(m-1)\cdot(n-1)\times(m-1)\cdot(n-1) $$. On tant $B$ i $I$ tenen $(m-1) $ columnes i $ (n-1)$ files.
$$
\begin{align*}
A = 
\begin{pmatrix}
B & I &0& \dots &0\\
I & B & I  &\dots  & 0 \\
0 & I & B   &\dots  & 0 \\
\vdots  &\ddots & \ddots & \ddots & \vdots\\
0  &\dots & 0 & I & B\\
\end{pmatrix} ; \quad 
B = 
\begin{pmatrix}
-4 & 1 & 0 & 0 & \dots &0\\
1 & -4 & 1 & 0 &  \dots & 0 \\
0 & 1 & -4 & 1  & \ddots & \vdots \\
0 & 0 & 1 & -4   & \ddots & 0 \\
\vdots  &\vdots &\ddots & \ddots & \ddots & 1\\
0 & 0  &\dots &0 & 1 & -4\\
\end{pmatrix}
I = 
\begin{pmatrix}
1 & 0 & \dots & 0 \\
0 & 1 & \dots & 0 \\
\vdots  &\vdots &\ddots &\vdots\\
0 & 0 & \dots & 1 \\
\end{pmatrix} 
\end{align*}
$$
Demostrat a l'annex que podem representar així el sistema d'equacions plantejat.

### Anàlisi de resultats

