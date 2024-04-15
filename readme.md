多元二次函数表示为

$f(x) = x^TAx+Bx+C$​

求极小值

$f'(x)=2Ax+B=0$

即求解

$Ax=-\frac{1}{2}B$

当 $\theta$ 固定时

$E(V) = E_S(V)+\lambda_LE_L(V) +\lambda_BE_B(V)$

$E(V)=[E_S,E_L,E_B]*\begin{bmatrix}1\\\lambda_L\\\lambda_B\end{bmatrix}$



$E_S(V) = \frac{1}{N_S}\sum||(A_q(A_q^TA_q)^{-1}A^{T}_{q}-I)V_q||^2$

令$A_q=A_q(A_q^TA_q)^{-1}A^{T}_{q}-I$

则$E_S(V)=\frac{1}{N_S}\sum(A_q*V_q)^T(A_qV_q)$

$E_S(V)=\sum\frac{1}{N_S}*(V_q^T(A^T_qA_q)V_q)$

$C = R\hat{e}(\hat{e}^T\hat{e})^{-1}\hat{e}^TR^T-I$

$E_L(V) = \frac{1}{N_L}\sum||C_j*e_{j}||^2$

$C_j$ 是线段j旋转的矩阵，$e_{j}$为$j$在所属的桶内的输出的双线性插值线段,$\hat{e}$为初始线段。

同样的

$E_L(V) = \sum\frac{1}{N_L}e_{j}^T(C_j^TC_j)*e_{j}$

e可以由V进行双线性插值得到$e_j=e_j'F_L$

$$\begin{bmatrix}
{\hat{x}_0}&{\hat{y}_0}&\hat{x}_0\hat{y}_0&1&0&0&0&0\\
{0}&{0}&0&0&{\hat{x}_0}&{\hat{y}_0}&\hat{x}_0\hat{y}_0&1\\
{\vdots}&{\vdots}&{\vdots}&{\vdots}&{\vdots}&{\vdots}&{\vdots}&{\vdots}&\\{0}&{0}&0&0&{\hat{x}_3}&{\hat{y}_3}&\hat{x}_{3}\hat{y}_3&1\\
\end{bmatrix}*F_L=\begin{bmatrix}x_0\\y_0\\\vdots\\{y_3}\end{bmatrix}$$

令左边这个为$k$则$e_j=e_j'k^{-1}*V_q$

更进一步令$K_j=e_j'k_j^{-1}$

则$E_L(V)=\sum\frac{1}{N_L}*(K_jV_q)^T(C_j^TC_j)K_jV_q$​

则$E_L(V)=\sum\frac{1}{N_L}*V_q^T(K_j^TC_j^TC_jK_j)V_q$

对于$E_B(V)$ 由于存在平方，需要求解的时候特殊约束处理，视作矩阵B

综上可知

求解$AV_q=B$

其中 $A=\sum\frac{A^T_qA_q}{N_S}+\lambda_L\sum\frac{K_j^TC_j^TC_jK_j}{N_L}+\lambda_BQ$

简而言之，就是把每一个网络、每一个线段、每一个约束条件看作一个线性方程计算Energy