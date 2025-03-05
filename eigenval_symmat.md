# Calculation for Eigenvalue and Eigenvector

## Eigen Value Calculation for 3X3 Symmetric Matrix

$$
\mathbf{A} = \begin{bmatrix} 
A_{11}  & A_{12} & A_{13} \\
A_{21} & A_{22}  & A_{23} \\
A_{31} & A_{32} & A_{33} 
\end{bmatrix}
$$

### Matrix Symmetry Check

The code first checks if the input matrix is symmetric within a tolerance(1e-6):

$$
|\mathbf{A}_{ij} - \mathbf{A}_{ji}|=0  \Rightarrow |\mathbf{A}_{ij} - \mathbf{A}_{ji}|<\text{tolerance}
$$

*Note: We will use a tolerance (a small value) when comparing numbers for equality or when checking if a number is close to zero. This is because floating-point arithmetic can introduce tiny rounding errors, making direct comparisons unreliable.*

```fortran
*GET, row_, PARM, %ARG1%, DIM, 1, 
*GET, col_, PARM, %ARG1%, DIM, 2, 
! Matrix symmetry check with TOLER_
*IF, row_, EQ, 3, AND, col_, EQ, 3, THEN
    *do,i,1,3
        *do,j,i+1,3
            xx_ = %ARG1%(i,j)
            YY_ = %ARG1%(j,i)
            *IF,abs(xx_-yy_),gt,TOLER_,then
                    Non_Sym_Flag = 1
                *exit
            *ENDIF
        *enddo
        *IF,Non_Sym_Flag,EQ,1,exit
    *enddo
*ELSE
    Non_Sym_Flag = 1
*ENDIF
```


### For a 3×3 Symmetric Matrix $\mathbf{A}$:

[Source: Use of Hyperbolic Cosines in Solving Cubic Polynomials](https://users.math.msu.edu/users/newhous7/math_235/lectures/cubic_gc_holmes.pdf)

The characteristic polynomial of a 3×3 matrix has the form:

$$
\lambda^3 - \text{Tr}(A)\lambda^2 + \text{pm_sum}\lambda - \det(A) = 0
$$

- `pm_sum`: principal minors
- `Tr(A)`: Trace of A
- `det(A)`: Determinant

1. **Matrix Properties Calculation**:
    - Trace: $\text{Tr}(\mathbf{A}) = \mathbf{A}_{11} + \mathbf{A}_{22} + \mathbf{A}_{33}$
    - Sum of principal minors:

$$
\text{pm_sum} = (\mathbf{A}_{22}\mathbf{A}_{33} - \mathbf{A}_{23}^2) + (\mathbf{A}_{11}\mathbf{A}_{33} - \mathbf{A}_{13}^2) + (\mathbf{A}_{11}\mathbf{A}_{22} - \mathbf{A}_{12}^2)
$$
- Determinant:
$$
\det(\mathbf{A}) = \mathbf{A}_{11}(\mathbf{A}_{22}\mathbf{A}_{33} - \mathbf{A}_{23}^2) - \mathbf{A}_{12}(\mathbf{A}_{12}\mathbf{A}_{33} - \mathbf{A}_{23}\mathbf{A}_{13}) + \mathbf{A}_{13}(\mathbf{A}_{12}\mathbf{A}_{23} - \mathbf{A}_{22}\mathbf{A}_{13})
$$
2. **Cubic Characteristic Equation Parameters**:
    - $p = \text{pm_sum} - \frac{\text{Tr}(\mathbf{A})^2}{3}$
    - $q = -\frac{2\text{Tr}(\mathbf{A})^3}{27} + \frac{\text{pm_sum} \cdot \text{Tr}(\mathbf{A})}{3} - \det(\mathbf{A})$
    - $\Delta = \left(\frac{q}{2}\right)^2 + \left(\frac{p}{3}\right)^3$
    - $r = \sqrt{-\frac{p}{3}}$
```fortran
    p_par = pm_sum - (Trace**2)/3
    q_par_ = (-2*Trace**3)/27 + (pm_sum*Trace)/3 - det
    sqrt_term_ = sqrt(-p_par/3)  !r
    delta = ((q_par_/2)**2) + (p_par/3)**3
```

3. **Eigenvalue Computation**:
    - If $r < \text{tolerance}$ (cubic equation has triple root):

$$
\lambda_1 = \lambda_2 = \lambda_3 = \frac{\text{Tr}(\mathbf{A})}{3}
$$
    - Otherwise, using trigonometric solution (when $\Delta \leq 0$):

$$
\cos\phi = \frac{-q}{2r^3} \quad \text{(clamped to } [-1,1] \text{)}
$$

$$
\phi = \arccos(\cos\phi)
$$

$$
\lambda_1 = 2r\cos\left(\frac{\phi}{3}\right) + \frac{\text{Tr}(\mathbf{A})}{3}
$$

$$
\lambda_2 = 2r\cos\left(\frac{\phi + 2\pi}{3}\right) + \frac{\text{Tr}(\mathbf{A})}{3}
$$

$$
\lambda_3 = 2r\cos\left(\frac{\phi + 4\pi}{3}\right) + \frac{\text{Tr}(\mathbf{A})}{3}
$$
```fortran
        *IF, sqrt_term_, LT, 1e-6, THEN
            ! Calculate eigenvalues
            %ARG2%(1) = Trace/3
            %ARG2%(2) = Trace/3
            %ARG2%(3) = Trace/3
        *ELSE
            cos_phi = -q_par_ /(2*sqrt_term_**3)
            ! Clamp cosine value
            *IF,cos_phi,GT,1,then
                cos_phi = 1
            *ELSEIF,cos_phi,LT,-1,then
                cos_phi = -1
            *ENDIF
            phi = ACOS(cos_phi)
            ! Calculate eigenvalues
            %ARG2%(1) = 2*sqrt_term_*cos(phi/3) + Trace/3
            %ARG2%(2) = 2*sqrt_term_*cos((phi + 2*pi)/3) + Trace/3
            %ARG2%(3) = 2*sqrt_term_*cos((phi + 4*pi)/3) + Trace/3
        *ENDIF
```

4. **Eigenvector Computation**:
For each eigenvalue $\lambda_i$:
    - Form the characteristic matrix: $\mathbf{C} = \mathbf{A} - \lambda_i \mathbf{I}$
    - Calculate cross products of rows to find the null space vector:
        - Try cross product of first two rows: $\vec{v} = \vec{r}_1 \times \vec{r}_2$
        - If $\|\vec{v}\| < \text{tolerance}$, try other row combinations
        - For degenerate cases, try cross products with unit vectors
    - Normalize the resulting vector: $\vec{e}_i = \frac{\vec{v}}{\|\vec{v}\|}$
```
flowchart TD
    A[Input 3×3 Matrix] --> B{Is Matrix Symmetric?}
    B -->|No| C[Error: Non-symmetric Matrix]
    B -->|Yes| D[Calculate Matrix Properties]
    D --> D1[Calculate Trace]
    D --> D2[Calculate Sum of Principal Minors]
    D --> D3[Calculate Determinant]
    D1 --> E[Solve Cubic Characteristic Equation]
    D2 --> E
    D3 --> E
    E --> E4[Calculate r = √-p/3]
    E4 --> F{Check Δ}
    F -->|Δ > 0| G[Error: Complex Roots]
    F -->|Δ ≤ 0| H{Is r < 1e-6}
    H -->|Yes| I[Triple Eigenvalue: λ₁=λ₂=λ₃=Trace/3]
    H -->|No| J[Calculate φ ]
    J --> K[Calculate Eigenvalues]
    K --> K1[λ₁]
    K --> K2[λ₂]
    K --> K3[λ₃]
    I --> L[Calculate Eigenvectors]
    K1 --> L
    K2 --> L
    K3 --> L
    L --> M[For each eigenvalue λ]
    M --> N[Form Matrix C = A - λI]
    N --> O[Calculate Cross Products of Rows]
    O --> P{Is Vector Norm > tolerance?}
    P -->|No| Q[Try Different Row Combinations]
    Q --> P
    P -->|Yes| R[Normalize Eigenvector]
    R --> S[Store Eigenvector]
    S --> T{All Eigenvalues Processed?}
    T -->|No| M
    T -->|Yes| U[Return Eigenvalues and Eigenvectors]
    V[end]
    U --> V
    C --> V
    G --> V
```


## Eigen Vector Calculation

[Source: method-of-finding-the-eigenvector](https://math.stackexchange.com/questions/989050/a-method-of-finding-the-eigenvector-that-i-dont-fully-understand)

For a symmetric 3×3 matrix $\mathbf{A}$, once we have computed its eigenvalues $\lambda_i$, we need to find the corresponding eigenvectors. The code implements a robust approach that handles various degenerate cases based on the rank of the characteristic matrix.

### Mathematical Foundation

An eigenvector $\mathbf{v}$ corresponding to eigenvalue $\lambda$ satisfies:

$$
\mathbf{A}\mathbf{v} = \lambda\mathbf{v}
$$

Rearranging:

$$
(\mathbf{A} - \lambda\mathbf{I})\mathbf{v} = \mathbf{0}
$$

This means that $\mathbf{v}$ lies in the null space (kernel) of the matrix $(\mathbf{A} - \lambda\mathbf{I})$. The null space can be computed by finding vectors orthogonal to the row space of $(\mathbf{A} - \lambda\mathbf{I})$.

### Eigenvector Computation Process

For each eigenvalue $\lambda_i$ (where $i = 1,2,3$), the code computes the corresponding eigenvector using the following steps:

1. **Form the Characteristic Matrix $\mathbf{C} = \mathbf{A} - \lambda_i\mathbf{I}$**

The characteristic matrix is created by subtracting the eigenvalue from the diagonal elements of matrix $\mathbf{A}$:

$$
\mathbf{C} = \begin{bmatrix} 
A_{11} - \lambda_i & A_{12} & A_{13} \\
A_{21} & A_{22} - \lambda_i & A_{23} \\
A_{31} & A_{32} & A_{33} - \lambda_i
\end{bmatrix}
$$

In the code, this is implemented as:

```fortran
char_mat_(1,1) = %ARG1%(1,1) - current_lambda
char_mat_(1,2) = %ARG1%(1,2)
char_mat_(1,3) = %ARG1%(1,3)
```

2. **Extract Rows from the Characteristic Matrix**

The code extracts the three rows from $\mathbf{C}$:

$$
\vec{r}_1 = [C_{11}, C_{12}, C_{13}]
$$

$$
\vec{r}_2 = [C_{21}, C_{22}, C_{23}]
$$

$$
\vec{r}_3 = [C_{31}, C_{32}, C_{33}]
$$

3. **Compute Eigenvector Based on Rank**

The algorithm determines eigenvectors using an approach that adapts to the rank of the characteristic matrix $\mathbf{C}$.

3.1. **Rank of $\mathbf{C}$ is 2**

When the rank is 2, the null space is one-dimensional. The eigenvector can be found by computing the cross product of any two linearly independent rows of $\mathbf{C}$.

**Attempt 1: Cross product of first two rows**

$$
\vec{v} = \vec{r}_1 \times \vec{r}_2
$$

In vector form:

$$
\vec{v} = \begin{bmatrix} 
r_{1,2}r_{2,3} - r_{1,3}r_{2,2} \\
r_{1,3}r_{2,1} - r_{1,1}r_{2,3} \\
r_{1,1}r_{2,2} - r_{1,2}r_{2,1}
\end{bmatrix}
$$

In the code:

```fortran
eig_vec(1) = row1(2)*row2(3) - row1(3)*row2(2)
eig_vec(2) = row1(3)*row2(1) - row1(1)*row2(3)
eig_vec(3) = row1(1)*row2(2) - row1(2)*row2(1)
```

**Attempt 2 \& 3: Try other row combinations**
If the norm of $\vec{v}$ is too small (< 1e-6), the code tries cross products of other row combinations:

- $\vec{r}_1 \times \vec{r}_3$
- $\vec{r}_2 \times \vec{r}_3$

3.2. **Rank of $\mathbf{C}$ is 1**

When the rank is 1, the null space is two-dimensional. The code creates two unit vectors and computes cross products with the remaining non-zero row.

Define unit vectors:

$$
\vec{u}_1 = [1, 0, 0]
$$

$$
\vec{u}_2 = [0, 1, 0]
$$

The code then computes:

- $\vec{r}_1 \times \vec{u}_1$
- $\vec{r}_1 \times \vec{u}_2$

This procedure finds vectors that are orthogonal to both $\vec{r}_1$ and one of the unit vectors, which span a subspace of the null space.

3.3. **Rank of $\mathbf{C}$ is 0**

In the extremely rare case where the rank is 0 (all entries of $\mathbf{C}$ are nearly zero, which can happen with repeated eigenvalues), the code uses the original matrix row:

$$
\vec{v} = [A_{i,1}, A_{i,2}, A_{i,3}]
$$

This is a fallback mechanism for the case where the matrix is a scalar multiple of the identity.

4. **Normalize the Eigenvector**

Once a non-zero eigenvector is found, it is normalized to have unit length:

$$
\|\vec{v}\| = \sqrt{v_1^2 + v_2^2 + v_3^2}
$$

$$
\vec{e} = \frac{\vec{v}}{\|\vec{v}\|}
$$

In the code:

```fortran
norm = SQRT(eig_vec(1)**2 + eig_vec(2)**2 + eig_vec(3)**2)
...
%ARG3%(j,i) = eig_vec(j)/norm
```


### Mathematical Justification

1. **Cross Product Method**: For a 3×3 matrix, the cross product of two rows produces a vector orthogonal to both rows. If these rows are linearly independent, the resulting vector lies in the null space of the matrix.
2. **Handling Degenerate Cases**: The algorithm systematically handles different rank cases:
    - Rank 2: One-dimensional null space, found via cross products
    - Rank 1: Two-dimensional null space, found via cross products with unit vectors
    - Rank 0: Three-dimensional null space, any vector is an eigenvector
3. **Normalization**: Ensures the eigenvectors have unit length, which is a standard convention.
```
flowchart TD
    A[Compute Eigenvalue λᵢ] --> B[Form Characteristic Matrix C = A - λᵢI]
    B --> C[Extract Rows from C]
    C --> D{Try Cross Product of Rows 1 & 2}
    D -->|Success: Norm ≥ 1e-6| J[Normalize Eigenvector]
    D -->|Failure: Norm < 1e-6| E{Try Cross Product of Rows 1 & 3}
    E -->|Success| J
    E -->|Failure| F{Try Cross Product of Rows 2 & 3}
    F -->|Success| J
    F -->|Failure| G[Handle Lower Rank Case]
    
    G --> H1[Create Unit Vectors]
    H1 --> H2{Try Cross Product of Row 1 & Unit Vector 1}
    H2 -->|Success| J
    H2 -->|Failure| H3{Try Cross Product of Row 1 & Unit Vector 2}
    H3 -->|Success| J
    H3 -->|Failure| I[Handle Null Rank Case]
    
    I --> I1[Use Original Matrix Row as Eigenvector]
    I1 --> J
    
    J --> K[Store Normalized Eigenvector]
    K --> L[Proceed to Next Eigenvalue]
```


### Conclusion

This approach for eigenvector computation is robust and handles various special cases that can arise, particularly with symmetric matrices that may have repeated eigenvalues. The algorithm intelligently adapts to the rank of the characteristic matrix to find appropriate eigenvectors in all situations.

---

## APDL Code Test

### **Usage**

#### **Syntax:**

```fortran
*use, <macro_file_path>, '<InputMatrix>', '<EigenValuesArray>', '<EigenVectorsArray>'
```


#### **Parameters:**

- `<InputMatrix>`: A **3×3 symmetric matrix**.
- `<EigenValuesArray>`: Output **array of eigenvalues**.
- `<EigenVectorsArray>`: Output **matrix where each column represents an eigenvector**.


#### **Example:1**

**APDL**

```fortran
*DEL,A,,NOPR
*dim,A,,3,3         
A(1,1) = 4 , 1,  2  
A(1,2) = 1 , 3 , 1  
A(1,3) = 2 , 1 , 5  
*use, C:\Users\pramod.kumar\Documents\TOOL_SymMat_Eigen.MAC,'A','EIGEN_VAL','EIGEN_VEC'

 PARAMETER STATUS- EIGEN_VAL
      LOCATION                VALUE
        1       1       1    7.04891734     
        2       1       1    2.30797853     
        3       1       1    2.64310413 
 PARAMETER STATUS- EIGEN_VEC
      LOCATION                VALUE
        1       1       1   0.591009049     
        2       1       1   0.327985278     
        3       1       1   0.736976229     
        1       2       1  -0.736976229     
        2       2       1   0.591009049     
        3       2       1   0.327985278     
        1       3       1   0.327985278     
        2       3       1   0.736976229     
        3       3       1  -0.591009049
```

**PYTHON**

```
eigenvalues: [7.04891734 2.30797853 2.64310413], 
eigenvectors:
[[ 0.59100905  0.73697623 -0.32798528]
 [ 0.32798528 -0.59100905 -0.73697623]
 [ 0.73697623 -0.32798528  0.59100905]]
```


#### **Example:2**

**APDL**

```fortran
 PARAMETER A(1,1) =    0.3065029800    -0.2826874000E-01-0.1374234000E-01
 PARAMETER A(1,2) =   -0.2826874000E-01 0.3047784800     0.2714881000E-01
 PARAMETER A(1,3) =   -0.1374234000E-01 0.2714881000E-01 0.2887970400

 USE MACRO FILE  C:\Users\pramod.kumar\Documents\TOOL_SymMat_Eigen.MAC
 ARG  1 = A
 ARG  2 = EIGEN_VAL
 ARG  3 = EIGEN_VEC

 PARAMETER STATUS- EIGEN_VAL
      LOCATION                VALUE
        1       1       1   0.348189320
        2       1       1   0.266968071
        3       1       1   0.284921109
 PARAMETER STATUS- EIGEN_VEC
      LOCATION                VALUE
        1       1       1  -0.598582695
        2       1       1   0.667121575
        3       1       1   0.443449616
        1       2       1  -0.246086187
        2       2       1  -0.679948281
        3       2       1   0.690732889
        1       3       1  -0.762325617
        2       3       1  -0.304333929
        3       3       1  -0.571174679
```

**PYTHON**

```python
eigenvalues: [0.34818932 0.28492111 0.26696807], 
eigenvectors:
[[ 0.59858269 -0.76232562 -0.24608619]
 [-0.66712157 -0.30433393 -0.67994828]
 [-0.44344962 -0.57117468  0.69073289]]
```

---

# Limiting Cases

### **1. All Data Lies on a Straight Line**

#### **Covariance Matrix**:

If the data lies on a straight line (e.g., along a direction vector **v**), the covariance matrix has **rank 1**. For example, if the line is along the direction $[1, 1, 1]^T$:

$$
\mathbf{\Sigma} = a \begin{bmatrix}
1 & 1 & 1 \\
1 & 1 & 1 \\
1 & 1 & 1
\end{bmatrix},
$$

where $a$ is the variance along the line.

#### **Eigenvalues**:

- **Non-zero eigenvalue**: $3a$ (variance along the line).
- **Zero eigenvalues**: $0, 0$ (no variance perpendicular to the line).


#### **PCA Significance**:

- The **first principal component (PC1)** aligns with the line’s direction.
- All data is perfectly collinear—**dimensionality reduces to 1D**.

```fortran
 PARAMETER A(1,1) =     1.000000000      1.000000000      1.000000000
 PARAMETER A(1,2) =     1.000000000      1.000000000      1.000000000
 PARAMETER A(1,3) =     1.000000000      1.000000000      1.000000000

 USE MACRO FILE  C:\Users\pramod.kumar\Documents\TOOL_SymMat_Eigen.MAC
 ARG  1 = A
 ARG  2 = EIGEN_VAL
 ARG  3 = EIGEN_VEC

 PARAMETER STATUS- EIGEN_VAL
      LOCATION                VALUE
        1       1       1    3.00000000
        2       1       1   4.440892099E-016
        3       1       1  -8.881784197E-016
 PARAMETER STATUS- EIGEN_VEC 
      LOCATION                VALUE
        1       1       1   0.577350269
        2       1       1   0.577350269
        3       1       1   0.577350269
        1       2       1    0.00000000
        2       2       1   0.707106781
        3       2       1  -0.707106781
        1       3       1    0.00000000
        2       3       1   0.707106781
        3       3       1  -0.707106781
```

---

### **2. All Data Lies on a Perfect Sphere (Evenly Distributed)**

#### **Covariance Matrix**:

For a sphere centered at the origin with radius $r$, the covariance matrix is **isotropic** (diagonal with equal variances):

$$
\mathbf{\Sigma} = \sigma^2 \begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix},
$$

where $\sigma^2 = \frac{r^2}{3}$ (variance in each axis direction).

#### **Eigenvalues**:

- All eigenvalues are **equal**: $\lambda_1 = \lambda_2 = \lambda_3 = \sigma^2$.


#### **PCA Significance**:

- No preferred direction—**all principal components are arbitrary**.
- The data is **isotropic**, so PCA does not reduce conditionality.

```fortran
 PARAMETER A(1,1) =     1.000000000      0.000000000      0.000000000
 PARAMETER A(1,2) =     0.000000000      1.000000000      0.000000000
 PARAMETER A(1,3) =     0.000000000      0.000000000      1.000000000

 USE MACRO FILE  C:\Users\pramod.kumar\Documents\TOOL_SymMat_Eigen.MAC
 ARG  1 = A
 ARG  2 = EIGEN_VAL
 ARG  3 = EIGEN_VEC
 
 PARAMETER STATUS- EIGEN_VAL
      LOCATION                VALUE
        1       1       1    1.00000000
        2       1       1    1.00000000
        3       1       1    1.00000000
        
 PARAMETER STATUS- EIGEN_VEC
      LOCATION                VALUE
        1       1       1    1.00000000
        2       1       1    0.00000000
        3       1       1    0.00000000
        1       2       1    0.00000000
        2       2       1    1.00000000
        3       2       1    0.00000000
        1       3       1    0.00000000
        2       3       1    0.00000000
        3       3       1    1.00000000

```

---

### **3. All Data Lies on a Cylinder Surface (Evenly Distributed)**

Assume the cylinder is aligned along the $z$-axis with radius $r$ and height $h$.

#### **Covariance Matrix**:

$$
\mathbf{\Sigma} = \begin{bmatrix}
\frac{r^2}{2} & 0 & 0 \\
0 & \frac{r^2}{2} & 0 \\
0 & 0 & \frac{h^2}{12}
\end{bmatrix},
$$

- **Radial variance** ($x,y$): $\frac{r^2}{2}$ (uniform distribution around the circle).
- **Axial variance** ($z$): $\frac{h^2}{12}$ (uniform distribution along the height).


#### **Eigenvalues**:

- **Two equal eigenvalues** (radial directions): $\lambda_1 = \lambda_2 = \frac{r^2}{2}$.
- **Third eigenvalue** (axial direction): $\lambda_3 = \frac{h^2}{12}$.


#### **PCA Significance**:

- **PC1 and PC2**: Capture radial variance (directions in the $x$-$y$ plane).
- **PC3**: Captures variance along the cylinder’s axis ($z$).
- Dimensionality depends on the task:
    - If $h \gg r$, PC3 dominates (dimensionality reduces to 1D along $z$).
    - If $r \gg h$, PC1/PC2 dominate (dimensionality reduces to 2D radial plane).

```fortran
 PARAMETER A(1,1) =     1.000000000      0.000000000      0.000000000
 PARAMETER A(1,2) =     0.000000000      1.000000000      0.000000000
 PARAMETER A(1,3) =     0.000000000      0.000000000      5.000000000

 USE MACRO FILE  C:\Users\pramod.kumar\Documents\TOOL_SymMat_Eigen.MAC
 ARG  1 = A
 ARG  2 = EIGEN_VAL
 ARG  3 = EIGEN_VEC

 PARAMETER STATUS- EIGEN_VAL
      LOCATION                VALUE
        1       1       1    5.00000000
        2       1       1    1.00000000
        3       1       1    1.00000000

 PARAMETER STATUS- EIGEN_VEC
      LOCATION                VALUE
        1       1       1    0.00000000
        2       1       1    0.00000000
        3       1       1    1.00000000
        1       2       1    0.00000000
        2       2       1    1.00000000
        3       2       1    0.00000000
        1       3       1    0.00000000
        2       3       1    0.00000000
        3       3       1    1.00000000
```

---

### **Key Takeaways**

| **Shape** | **Covariance Matrix** | **Eigenvalues** | **PCA Interpretation** |
| :-- | :-- | :-- | :-- |
| Straight Line | Rank 1 (all rows/columns identical) | `(3a, 0, 0)` | 1D structure; perfect collinearity |
| Sphere | Diagonal with equal entries | `(sigma^2, sigma^2, sigma^2)` | No structure; isotropic variance |
| Cylinder Surface | Diagonal with two equal entries + one distinct | `(r^2/2, r^2/2, h^2/12)` | 2D radial + 1D axial structure |

These results highlight how geometric structure in data directly determines the covariance matrix’s eigenvalues and the principal components in PCA.

---

## **Why Are Eigenvalues of a Symmetric Matrix Always Real?**

[Source: Symmetric matrices have real eigenvalues](https://pi.math.cornell.edu/~jerison/math2940/real-eigenvalues.pdf)

The Spectral Theorem states that if $A$ is an $n \times n$ symmetric matrix with real entries, then it has $n$ orthogonal eigenvectors. The first step of the proof is to show that all the roots of the characteristic polynomial of $A$ (i.e. the eigenvalues of $A$ ) are real numbers.

Recall that if $z=a+b i$ is a complex number, its complex conjugate is defined by $\bar{z}=a-b i$. We have $z \bar{z}=(a+b i)(a-b i)=a^{2}+b^{2}$, so $z \bar{z}$ is always a nonnegative real number (and equals 0 only when $z=0$ ). It is also true that if $w, z$ are complex numbers, then $\overline{w z}=\bar{w} \bar{z}$.

Let $\mathbf{v}$ be a vector whose entries are allowed to be complex. It is no longer true that $\mathbf{v} \cdot \mathbf{v} \geq 0$ with equality only when $\mathbf{v}=\mathbf{0}$. For example,

$$
\left[\begin{array}{l}
1 \\
i
\end{array}\right] \cdot\left[\begin{array}{l}
1 \\
i
\end{array}\right]=1+i^{2}=0
$$

However, if $\overline{\mathbf{v}}$ is the complex conjugate of $\mathbf{v}$, it is true that $\overline{\mathbf{v}} \cdot \mathbf{v} \geq 0$ with equality only when $\mathbf{v}=\mathbf{0}$. Indeed,

$$
\left[\begin{array}{c}
a_{1}-b_{1} i \\
a_{2}-b_{2} i \\
\vdots \\
a_{n}-b_{n} i
\end{array}\right] \cdot\left[\begin{array}{c}
a_{1}+b_{1} i \\
a_{2}+b_{2} i \\
\vdots \\
a_{n}+b_{n} i
\end{array}\right]=\left(a_{1}^{2}+b_{1}^{2}\right)+\left(a_{2}^{2}+b_{2}^{2}\right)+\cdots+\left(a_{n}^{2}+b_{n}^{2}\right)
$$

which is always nonnegative and equals zero only when all the entries $a_{i}$ and $b_{i}$ are zero.
With this in mind, suppose that $\lambda$ is a (possibly complex) eigenvalue of the real symmetric matrix $A$. Thus there is a nonzero vector $\mathbf{v}$, also with complex entries, such that $A \mathbf{v}=\lambda \mathbf{v}$. By taking the complex conjugate of both sides, and noting that $\bar{A}=A$ since $A$ has real entries, we get $\overline{A \mathbf{v}}=\overline{\lambda \mathbf{v}} \Rightarrow A \overline{\mathbf{v}}=\bar{\lambda} \overline{\mathbf{v}}$. Then, using that $A^{T}=A$,

$$
\begin{aligned}
\overline{\mathbf{v}}^{T} A \mathbf{v} & =\overline{\mathbf{v}}^{T}(A \mathbf{v})=\overline{\mathbf{v}}^{T}(\lambda \mathbf{v})=\lambda(\overline{\mathbf{v}} \cdot \mathbf{v}) \\
\overline{\mathbf{v}}^{T} A \mathbf{v} & =(A \overline{\mathbf{v}})^{T} \mathbf{v}=(\bar{\lambda} \overline{\mathbf{v}})^{T} \mathbf{v}=\bar{\lambda}(\overline{\mathbf{v}} \cdot \mathbf{v})
\end{aligned}
$$

Since $\mathbf{v} \neq \mathbf{0}$, we have $\overline{\mathbf{v}} \cdot \mathbf{v} \neq 0$. Thus $\lambda=\bar{\lambda}$, which means $\lambda \in \mathbf{R}$.

### **Conclusion**

A **real symmetric matrix always has real eigenvalues**. So, there is **no possibility of an imaginary eigenvalue** for a symmetric matrix. If a matrix has complex eigenvalues, then it must be **non-symmetric**.

---

## Plane Fit Macro Test

Here we will test the plane fit macro:

* Select node cloud from an existing database
* Fetch nodal location coordinates into an array
* Export array data to a text file for external validation
* Use the same data as input for the plane fit macro


### Test-1

- Here we selected component TA (DDTT2 TR02 GR04)

**APDL result**

```fortran
 NAME                              VALUE                        TYPE  DIMENSIONS
 NUMNODE                        40496.0000                  SCALAR

 USE MACRO FILE  C:\Users\pramod.kumar\Documents\NewFolder\Tool_FitPlanePCA.mac
 ARG  1 = COORD_COMP_BF
 ARG  2 =    1.0000
 ARG  3 = out1

 PARAMETER STATUS- OUT1_CENT
      LOCATION                VALUE
        1       1       1   -4.80176175
        1       2       1   -2.65859686
        1       3       1   3.413642991E-003

 PARAMETER STATUS- OUT1_NORMAL
      LOCATION                VALUE
        1       1       1   0.103283011
        1       2       1   0.994649657
        1       3       1  -2.163258547E-003

 PARAMETER STATUS- OUT1_DIST
 NAME                              VALUE                        TYPE  DIMENSIONS
 OUT1_DIST                         3.14032026                  SCALAR
```

**Python result**

```python
Number of points: 40496
Centr: [-4.80176175e+00 -2.65859686e+00  3.41364298e-03]
normal:[-0.10328301 -0.99464966  0.00216326]
Distance: -3.1403202555184895
```

**Plot**

### Test-2

**APDL result**

```fortran
 USE MACRO FILE  C:\Users\pramod.kumar\Documents\NewFolder\Tool_FitPlanePCA.mac
 ARG  1 = COORD_COMP_BF
 ARG  2 =    1.0000
 ARG  3 = out1

 SG5X DRB+ TR00 GR00 GI24

 PARAMETER STATUS- OUT1_CENT
      LOCATION                VALUE
        1       1       1   -9.53569460
        1       2       1    1.42200000
        1       3       1   0.749431723

 PARAMETER STATUS- OUT1_NORMAL 
      LOCATION                VALUE
        1       1       1  -2.205791392E-015
        1       2       1    1.00000000
        1       3       1  -1.771444401E-018

 PARAMETER STATUS- OUT1_DIST 
 NAME                              VALUE                        TYPE  DIMENSIONS
 OUT1_DIST                        -1.42200000                  SCALAR
```

**Python result**

```
number of points: 1333
Centro: [-9.5356946   1.422       0.74943172]
normal:[0. 1. 0.]
Distance: -1.4220000000000002
```

**Plot**

### Test-3

**APDL result**

```
 USE MACRO FILE  C:\Users\pramod.kumar\Documents\NewFolder\Tool_FitPlanePCA.mac
 ARG  1 = COORD_COMP_BF
 ARG  2 =    1.0000
 ARG  3 = out1
 
 PARAMETER STATUS- OUT1_CENT
      LOCATION                VALUE
        1       1       1   -6.25710923
        1       2       1   3.381026923E-005
        1       3       1  -1.178140284E-004

 PARAMETER STATUS- OUT1_NORMAL
      LOCATION                VALUE
        1       1       1  -0.999999929
        1       2       1   3.513700982E-004
        1       3       1  -1.339453027E-004

 PARAMETER STATUS- OUT1_DIST
 NAME                              VALUE                        TYPE  DIMENSIONS
 OUT1_DIST                        -6.25710882                  SCALAR
```

**Python result**

```python
number of points: 1018
Centro: [-6.25710923e+00  3.38102750e-05 -1.17813969e-04]
normal:[ 9.99999929e-01 -3.51370287e-04  1.33945324e-04]
Distance: 6.257108818143692
```

**Plot**

### Test-4

**APDL result**

```
 USE MACRO FILE  C:\Users\pramod.kumar\Documents\NewFolder\Tool_FitPlanePCA.mac
 ARG  1 = COORD_COMP_BF
 ARG  2 =    1.0000
 ARG  3 = out1

 PARAMETER STATUS- OUT1_CENT  (   9428 PARAMETERS DEFINED)
                  (INCLUDING       22 INTERNAL PARAMETERS)

      LOCATION                VALUE
        1       1       1   -8.08733348
        1       2       1   0.381147282
        1       3       1   4.987598274E-002

 PARAMETER STATUS- OUT1_NORMAL  (   9428 PARAMETERS DEFINED)
                  (INCLUDING       22 INTERNAL PARAMETERS)

      LOCATION                VALUE
        1       1       1  -1.368165755E-002
        1       2       1  -0.997867922
        1       3       1  -6.381553241E-002

 PARAMETER STATUS- OUT1_DIST
 NAME                              VALUE                        TYPE  DIMENSIONS
 OUT1_DIST                        0.272869382                  SCALAR
```

**Python result**

```
number of points: 77102
Centro: [-8.08733348  0.38114728  0.04987598]
normal:[-0.01368166 -0.99786792 -0.06381553]
Distance: 0.272869381882859
```

