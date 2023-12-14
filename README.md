# Linearequation<br>

**개요**<br>
* [gaussElim( )]()<br>
* [LUdecomp( )]()<br>
* [solveLU() )]()<br>



<hr>

## gaussElim( )<br>

```c
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d)
```
>Parameter<br>
**A** : 계수가 저장된 배열 <br>
**b** : vector b <br>
**U** : 상삼각 행렬 <br>
**d** : 해를 저장할 배열

1. myMatrix.h 안에 선언되어 있다.
2. 다음과 같은 원리를 기반으로 작성되었다.<br>

$$f(x_k)=\sum_{i=0}^ny_i=\prod_{j=0,j\neq i}^n{(x-x_j)\over (x_i-x_j)}$$

* Step1 : First pivot

$$
\begin{equation}
  \begin{bmatrix}
  R_1 \rightarrow R_1\\
  R_2 \rightarrow R_2-({a_{21}\over a_{11}})R_1\\
  R_3 \rightarrow R_3-({a_{31}\over a_{11}})R_1\\
  R_4 \rightarrow R_4-({a_{41}\over a_{11}})R_1\\
  \end{bmatrix}\Rightarrow
  \begin{bmatrix}
  R_1^{'}\\
  R_2^{'} \\
  R_3^{'} \\
  R_4^{'}\\
  \end{bmatrix}
\end{equation}
$$

* Step2 : Second pivot

$$
\begin{equation}
  \begin{bmatrix}
  R_1^{'} \rightarrow R_1^{'}\\
  R_2^{'} \rightarrow R_2^{'}\\
  R_3^{'} \rightarrow R_3^{'}-({a_{32}^{'}\over a_{22}^{'}})R_2^{'}\\
  R_4^{'} \rightarrow R_4^{'}-({a_{42}^{'}\over a_{22}^{'}})R_2^{'}\\
  \end{bmatrix}\Rightarrow
  \begin{bmatrix}
  R_1^{''}\\
  R_2^{''} \\
  R_3^{''} \\
  R_4^{''}\\
  \end{bmatrix}
\end{equation}
$$


* Repeat until the system become a Upper triangular matrix.

* Backsubstitution

$$x_4={b_4 \over a_{44}}$$
$$x_3={b_3 - {a_{34}b_4 \over a_{44}} \over a_{33}}$$

$$x_2={{b_2-a_{23}x_3-a_{24}x_4} \over a_{22}}$$

$$x_1={{b_1a_{22}-a_{13}x_3-a_{14}x_4-a_{12}x_2}\over a_{11}}$$

<hr>

## Example <br>
```c++
int main()
{
	Matrix prob1_matA = txt2Mat("C:\\NP_Data\\Assignment5\\", "prob1_matA");
	printMat(prob1_matA, "Matrix A");
	Matrix prob1_vecb = txt2Mat("C:\\NP_Data\\Assignment5\\", "prob1_vecb");
	printMat(prob1_vecb, "Matrix b");

	
	Matrix U = createMat(prob1_matA.rows, prob1_matA.cols);
	Matrix d = createMat(prob1_vecb.rows, 1); 
	gaussElim(prob1_matA, prob1_vecb, U, d);

	printMat(U, "Matrix U");
	printMat(d, "Matrix d");

	Matrix x = createMat(prob1_vecb.rows, 1);
	backsub(U, d, x);
	printMat(x, "Matrix x");
}
```

## Output <br>
```c
Matrix A =
       9.500000       -2.500000        0.000000       -2.000000        0.000000
      -2.500000       11.000000       -3.500000        0.000000       -5.000000
       0.000000       -3.500000       15.500000        0.000000       -4.000000
      -2.000000        0.000000        0.000000        7.000000       -3.000000
       0.000000       -5.000000       -4.000000       -3.000000       12.000000

Matrix b =
      12.000000
     -16.000000
      14.000000
      10.000000
     -30.000000

Matrix U =
       9.500000       -2.500000        0.000000       -2.000000        0.000000
       0.000000       10.342105       -3.500000       -0.526316       -5.000000
       0.000000        0.000000       14.315522       -0.178117       -5.692112
       0.000000        0.000000        0.000000        6.549947       -3.325276
       0.000000        0.000000        0.000000        0.000000        5.631235
________________________________________________________
Matrix d =
      12.000000
     -12.842105
       9.653944
      11.992890
     -26.281520

Matrix x =
       1.408467
      -0.912619
       0.697151
       1.830991
      -4.667097
```

## Warning
>A의 행의 개수와 U의 열의 개수가 같아야 한다.<br>
A의 열의 개수와 b의 열의 개수가 같아야 한다.
## Error Handling
```c
if (A.rows != b.rows || A.cols != U.rows )
{
	printf("Error : A의 행의 개수와 b의 행의 개수가 같지 않거나 A의 행의 개수와 U의 열의 개수가 같지 않습니다.\n");

	return;
}
```
<br>


## LUdecomp( )<br>

```c
void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P) 
```
>Parameter<br>
**A** : LU분해할 행렬 <br>
**L** : lower triangular matrix <br>
**U** : 상삼각 행렬 <br>

1. myMatrix.h 안에 선언되어 있다.
2. 다음과 같은 원리를 기반으로 작성되었다.<br>

* k=1
  
$$A^{(1)}=M^{(1)}A^{(0)}=M^{(1)}A$$
$$A^{(2)}=M^{(2)}A^{(1)}=M^{(2)}M^{(1)}A$$
$$A^{(k)}=M^{(k)}A^{(k-1)}=M^{(k)}M^{(k-1)}\cdots M^{(1)}A$$
$$U=A^{(m-1)}=M^{(m-1)}A^{(m-2)}=M^{(m-1)}M^{(m-2)}\cdots M^{(1)}A$$

* Therefore, the final row-echeleon form U is express in terms of original A as

$$U=M^{(m-1)}M^{(m-2)}\cdots M^{(1)}A$$
$$A=(M^{(m-1)}M^{(m-2)}\cdots M^{(1)})^{-1}U$$
$$=LU$$

* where
$$m_{21}={a_{21}\over a_{11}},m_{31}={a_{31}\over a_{11}},m_{32}={a_{32}\over a_{22}}\cdots $$

<hr>

## Example <br>
```c++
int main()
{
	Matrix A = txt2Mat(path, "prob1_matK");
	printMat(A, "Matrix A");
	Matrix b = txt2Mat(path, "prob1_vecf");
	printMat(b, "Matrix b");


	


	Matrix L = eye(A.rows, A.cols);
	Matrix U = eye(A.rows, A.cols);
	Matrix P = eye(A.rows, A.cols);


	LUdecomp(A, L, U, P);

	Matrix x = zeros(b.rows, b.cols);


	solveLU(L, U, P, b, x);

	printMat(x, "Matirx x");
}
```

## Output <br>
```c
Matrix A =
      75.000000      -20.000000        0.000000
     -20.000000       35.000000      -15.000000
       0.000000      -15.000000       15.000000

Matrix b =
      19.620001
      29.430000
      14.715000

Matirx x =
       1.159364
       3.366614
       4.347614
```

## Warning
>A와 L,U의 행렬 크기는 같아야 한다.
A는 정방행렬이어야 한다.
## Error Handling
```c
if (A.rows != A.cols)
{
	printf("Error :  A is not nxn square matrix\n");

	return;
}

if (A.rows != L.rows || A.rows != U.rows || A.cols != L.cols || A.cols != U.cols)
{
	printf("Error : The size of A and L is not the ame or The size of A and U is not the same\n");

	return;
}
```
<br>

## solve LU( )<br>

```c
void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x);
```
>Parameter<br>
**L** : LU분해한 행렬 <br>
**U** : LU분해한 행렬 <br>
**b** : vector <br>
**x** : 구하고자하는 행렬 <br>

1. myMatrix.h 안에 선언되어 있다.
2. 다음과 같은 원리를 기반으로 작성되었다.<br>

  
$$PAx-Pb=LUx=b$$
$$PA=LU$$
$$L(Ux)=Ly=Pb=b$$

* step1 : Solve for y in Ly=b using Forward Sub

$$
\begin{equation}
  \begin{bmatrix}
  l_{11} & 0 & 0 & \cdots\\
  l_{21} & l_{22} & 0 & \cdots\\
  l_{31} & l_{32} & l_{33} & \cdots\\
  l_{41} & l_{42} & l_{43} & \cdots\\
  \end{bmatrix}
  \begin{bmatrix}
  y_1\\
  y_2\\
  \vdots
  \end{bmatrix}=
  \begin{bmatrix}
  b_1\\
  b_2\\
  \vdots
  \end{bmatrix}=
\end{equation}
$$

$$
y_1={b_1\over l_{11}}$$
$$y_2={(b_2-l_{21}y_1)\over l_{22}}$$
$$\vdots$$

* step2 : Solve for x, in Ux=y using Backward sub

$$
\begin{equation}
  \begin{bmatrix}
  U_{11} & U_{12} & U_{13} & \cdots\\
  0 & U_{22} & U_{23} & \cdots\\
  0 & 0 & U_{33} & \cdots\\
  0 & 0 & 0 & \cdots\\
  \end{bmatrix}
  \begin{bmatrix}
  x_1\\
  x_2\\
  \vdots
  \end{bmatrix}=
  \begin{bmatrix}
  y_1\\
  y_2\\
  \vdots
  \end{bmatrix}=
\end{equation}
$$

$$
x_n={y_n\over U_{nn}}$$
$$y_{n-1}={(y_{n-1}-U_{nn-1}y_n)\over U_{n-1n-1}}$$
$$\vdots$$

<hr>

## Example <br>
```c++
int main()
{
	Matrix A = txt2Mat(path, "prob1_matK");
	printMat(A, "Matrix A");
	Matrix b = txt2Mat(path, "prob1_vecf");
	printMat(b, "Matrix b");


	


	Matrix L = eye(A.rows, A.cols);
	Matrix U = eye(A.rows, A.cols);
	Matrix P = eye(A.rows, A.cols);


	LUdecomp(A, L, U, P);

	Matrix x = zeros(b.rows, b.cols);


	solveLU(L, U, P, b, x);

	printMat(x, "Matirx x");
}
```

## Output <br>
```c
Matrix A =
      75.000000      -20.000000        0.000000
     -20.000000       35.000000      -15.000000
       0.000000      -15.000000       15.000000

Matrix b =
      19.620001
      29.430000
      14.715000

Matirx x =
       1.159364
       3.366614
       4.347614
```

## Warning
>L과 U는 정방행렬이어야 하고 사이즈가 같아야 한다.
b의 열과 L의 열은 같아야 하고 b는 벡터이다 b의 열과 x의 열은 같아야 하고 x는 벡터이다.<br>
## Error Handling
```c
if (L.cols != L.rows || L.rows != U.cols || L.rows != U.rows || b.rows != L.rows || b.cols != 1 || b.rows != x.rows || b.cols != x.cols)
{
	printf("Error\n");

	return;
}
```
<br>
