# Memòria

`Bernat Esquirol Juanola`

### Treball previ

Ens diuen que hem d'utilitzar el [mètode de diferències finites](https://en.wikipedia.org/wiki/Finite_difference_method) Per resoldre:

 _Equació de Poisson_


$$
\ u_{xx}(x,y)+u_{yy}(x,y)=f(x,y) ; \forall (x,y)\in D\equiv [a,b]\times[c,d]\in \mathbb{R}
$$
El que hem de fer primer de tot és tractar el domini $D​$, per tal de convertir-lo en un domini discretitzat, així doncs considerem la graella de punts $(x_i,y_j)​$ resultants de dividir $[a,b]​$ en $n ​$ subintervals de mida $hx​$ i $[c,d]​$ en $m​$ de mida $hy​$. D'aquesta manera trobarem només els $u_{i,j}=u(x_i,y_j)​$, dels punts de la graella, amb $0<i<n ​$ , $0<j<m ​$ (no incloem els extrems ja que ja tenim els seus valors).

Fixem-nos amb quines són les variables i quines les incògnites del problema: coneixem la $f(x,y)$, l'interval $D$ i $u(x,y)$ a $\partial D$ i volem trobar tots els $u_{i,j}$ de la graella.

Tenim per una banda l'equació de Poisson que ens relaciona $u_{xx}, u_{yy}$ amb $f$ i per altra banda les diferències centrals finites que ens relacionen numèricament $f''$ amb $f$.

_Diferències centrals finites_
$$
p''(x_0)\approx\frac{p(x_0+h)-2p(x_0)+p(x_0-h)}{h^2}
$$


Si apliquem això a les columnes de la nostra graella tindrem una aproximació a $u_{xx}$ en el punt $(x_i,y_j)$, en funció de $u_{i-1,j}$,$ u_{i,j}$ i $ u_{i+1,j}$, i si ho fem a les columnes tindrem  $u_{yy}$ en funció de $u_{i,j-1}$,$ u_{i,j}$ i $ u_{i,j+1}$.
$$
u_{xx}(x_i,y_j)\approx\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h_x^2}
\\
u_{yy}(x_i,y_j)\approx\frac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{h_y^2}
$$
i ara,
$$
u_{xx}(x,y)+u_{yy}(x,y)=f(x,y) \Rightarrow \frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h_x^2}+\frac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{h_y^2}\approx f(x_i, y_j)
\\
$$
Per tant tenim la següent equació de 5 incògnites, amb $h$ i $f$ coneguts. 
$$
\begin{equation} \label{eq:1}
\frac{1}{h_x^2}u_{i+1,j}+\frac{1}{h_x^2}u_{i-1,j}-\Big(\frac{2}{h_x^2}+\frac{2}{h_y^2}\Big)u_{i,j}+\frac{1}{h_y^2}u_{i,j+1}+\frac{1}{h_y^2}u_{i,j-1}\approx  f(x_i, y_j)
\tag{1}
\end{equation}
$$
Ara si suposem que $m$ i $n$ són més grans de 5, tindrem un sistema d'equacions d'on podrem resoldre numèricament $u_{i,j} $  $\forall (i,j)$; $0<i<n$, $0<j<m$, amb $f(x_i,y_i)$ coneguda. Volem expressar aquest sistema d'equacions com un sistema de la forma, $A\cdot u=b$. Per fer-ho hem de posar els valors de la graella $(n-1)\times (m-1)$ que representa $D$ en format de vector. De manera que $u$:
$$
u = 
\begin{pmatrix}
u_{1,1} & u_{2,1}  & \dots &u_{n-1,1}\\
u_{1,2} & u_{2,2}  & \dots &u_{n-1,2}\\
\vdots &\vdots &\ddots &\vdots \\
u_{1,m-1} & u_{1,m-1} & \dots &u_{n-1,m-1}
\end{pmatrix}
$$
Ens quedarà com a vector $(m-1)\cdot(n-1)$:
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
Definim el vector $b$ com el els valors de l'equació en cada una de les incògnites, que l'expressem en forma de matriu $(m-1)\times(n-1)$:
$$
b = 
\begin{pmatrix} 
f_{1,1}-\frac{g_{1,0}}{h_y^2}-\frac{g_{0,1}}{h_x^2} & f_{2,1}-\frac{g_{2,0}}{h_y^2} & \dots & f_{n-2,1}-\frac{g_{n-2,0}}{h_y^2} & f_{n-1,1}-\frac{g_{n-1,0}}{h_y^2}-\frac{g_{n,1}}{h_x^2} \\
f_{1,2}-\frac{g_{0,2}}{h_x^2} & f_{2,1} & \dots & f_{n-2,2} & f_{n-1,2}-\frac{g_{n,2}}{h_x^2} \\
\vdots  &\vdots &\ddots & \vdots & \vdots \\
f_{1,m-2}-\frac{g_{0,m-2}}{h_x^2} & f_{2,m-2} & \dots & f_{n-2,m-2} & f_{n-2,m-1}-\frac{g_{n,m-2}}{h_x^2} \\
f_{1,m-1}-\frac{g_{1,m}}{h_y^2}-\frac{g_{0,m-1}}{h_x^2} & f_{2,m-1}-\frac{g_{2,m}}{h_y^2} & \dots & f_{n-2,m-1}-\frac{g_{n-2,m}}{h_y^2} & f_{n-1,m-1}-\frac{g_{n-1,m}}{h_y^2}-\frac{g_{n,m-1}}{h_x^2} \\
\end{pmatrix}
$$

Veiem que a les vores tenim una incògnita resolta (o dues als extrems), ja que sabem el valor de la frontera $g$ en $\partial D$. També el convertim en vector $(m-1)\cdot(n-1)$:
$$
b = 
\begin{pmatrix}
f_{1,1}-\frac{g_{0,1}}{h_x^2}-\frac{g_{1,0}}{h_y^2} &
f_{1,2}-\frac{g_{0,2}}{h_x^2} &
\dots &
f_{2,1}-\frac{g_{2,0}}{h_y^2} &
f_{2,2} &
\dots &
f_{n-1,m-1}-\frac{g_{n-1,m}}{h_y^2}-\frac{g_{n,m-1}}{h_x^2}
\end{pmatrix}^\top
$$
D'aquesta manera podem definir a partir de $(\ref{eq:1})$ la matriu $A$ que relaciona incògnites $u$ amb valors $f$ una matriu de $$(m-1)\cdot(n-1)\times(m-1)\cdot(n-1) $$. On tant $B$ i $I$ tenen $(m-1) $ columnes i $ (n-1)$ files.
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
-\Big(\frac{2}{h_y^2}+\frac{2}{h_x^2}\Big) & \frac{1}{h_y^2} & 0 & 0 & \dots &0\\
\frac{1}{h_y^2} & -\Big(\frac{2}{h_y^2}+\frac{2}{h_x^2}\Big) & \frac{1}{h_y^2} & 0 &  \dots & 0 \\
0 & \frac{1}{h_y^2} & -\Big(\frac{2}{h_y^2}+\frac{2}{h_x^2}\Big) & \frac{1}{h_y^2}  & \ddots & \vdots \\
0 & 0 & \frac{1}{h_y^2} & -\Big(\frac{2}{h_y^2}+\frac{2}{h_x^2}\Big)   & \ddots & 0 \\
\vdots  &\vdots &\ddots & \ddots & \ddots & \
\frac{1}{h_y^2}\\
0 & 0  &\dots &0 & \frac{1}{h_y^2} & -\Big(\frac{2}{h_y^2}+\frac{2}{h_x^2}\Big)\\
\end{pmatrix}
I = 
\begin{pmatrix}
\frac{1}{h_x^2} & 0 & \dots & 0 \\
0 & \frac{1}{h_x^2} & \dots & 0 \\
\vdots  &\vdots &\ddots &\vdots\\
0 & 0 & \dots & \frac{1}{h_x^2} \\
\end{pmatrix} 
\end{align*}
$$
Es pot comprovar que $A\times u=b$, ho fem amb la primera fila de $A$:

$$
\begin{pmatrix}
-\Big(\frac{2}{h_y^2}+\frac{2}{h_x^2}\Big) & \frac{1}{h_y^2} & 0 & \dots &0&\frac{1}{h_x^2} &0&\dots &0
\end{pmatrix}
\begin{pmatrix}
u_{1,1} \\ u_{1,2}  \\u_{1,3} \\\vdots \\u_{1,m-1} \\u_{2,1} \\ u_{2,2} \\\vdots\\ u_{n-1,m-1}
\end{pmatrix}
=
f_{1,1}-\frac{g_{0,1}}{h_x^2}-\frac{g_{1,0}}{h_y^2}
$$

$$
-\Big(\frac{2}{h_y^2}+\frac{2}{h_x^2}\Big)u_{1,1}+\frac{u_{1,2}}{h_y^2}+\frac{u_{2,1}}{h_x^2}=f_{1,1}-\frac{g_{0,1}}{h_x^2}-\frac{g_{1,0}}{h_y^2} 
$$

Això és exactament el sistema que teniem.

### Anàlisi de resultats

A la pràctica s'ha implementat aquests vectors com a matrius $m-1\times n-1$ amb l'ennumeració natural $(1,...n-1)$. Els algoritmes són ja adaptats al problema i no funcionals per qualsevol matriu.La matriu d'incògnites és $w$. De manera que $u_{1,1}=w[1][1]$. El punt $u_{1,1}$ és l'extrem el punt més pròxim a $(a, c)$,  $u_{1,j}$ recorra l'eix

### Proves

He comprovat amb la graella $n=6; m=5$ que el valor de $\omega=1.4$ és el que aconsegueix minimitzar l'error (tolerance=$1\times10^{-10}$) i el nombre d'iteracions (max=$1000$) fetes. Per tant farem les comprovacions modificant $n$ i $m$ amb aquest valor.

|             | 5x6     |           | 10x12     |           | 25x30     |                | 50x60     |                |
| ----------- | ------- | --------- | --------- | --------- | --------- | -------------- | --------- | -------------- |
| $J$         | 145     | tolerance | 562       | tolerance | max       | 0.0000463177   | max       | 0.0010650960   |
| $GS$        | 76      | tolerance | 291       | tolerance | max       | _0.0000000008_ | max       | _0.0001008817_ |
| $SOR_{1,4}$ | _41_    | tolerance | _208_     | tolerance | max       | 0.0000000033   | max       | 0.0001877131   |
|             | **6x5** |           | **12x10** |           | **30x25** |                | **60x50** |                |
| $J$         | 116     | tolerance | 455       | tolerance | max       | 0.0000138072   | max       | 0.0008877161   |
| $GS$        | 61      | tolerance | 237       | tolerance | max       | 0.0000000195   | max       | 0.0002667446   |
| $SOR_{1,4}$ | _41_    | tolerance | _169_     | tolerance | _987_     | tolerance      | max       | _0.0000864609_ |

|           | $SOR_{1,4}$ |
| --------- | ----------- |
| **8x25**  | 754         |
| **10x20** | 508         |
| **14x14** | 299         |
| **20x10** | *247*       |
| **25x8**  | 273         |

### Conclusions

Si mirem els resultats obtinguts en podem treure dos conclusions:

1. $SOR_{1,4}$ és clarament l'algoritme més eficient, només en el cas concret de tenir una $m$ gran i una $n$ més petita està per sota de $GS$, com veiem a la primera taula.
2. Tots els algoritmes tendeixen a tenir menys error quan la $n$ i $m$ provoquen una graella de quadrats enlloc de rectangles. És a dir, com que tenim un rectangle de $2\times1$ a l'exemple, el més eficient possible és triar una $n=2m$. A la segona taula podem veure com de totes les combinacions que tenen el mateix nombre de caselles dins la graella ($\approx200$) la que troba la solució amb menys iteracions és la que compleix $n=2m$.

Per tant per tenir màxima rapidesa de convergència hauriem d'utilitzar el mètode SOR amb $\omega=1.4$ en una graella de $2m\times m$, en l'exemple donat.