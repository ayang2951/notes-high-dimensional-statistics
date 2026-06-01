_Overview: We carefully derive the high-dimensional asymptotic results for the single-component spiked matrix model._

## Single-Component Spiked Matrix Model

We'll spend a lot of time here, on the single-component spiked matrix model. We'll discuss a few key interesting properties of this problem, especially the phase transition, and highlight the application of some important tools.

### Problem Setup

Suppose we have
$$
X = \mu u v_*^\top + Z,
$$
where $X \in \mathbb R^{n \times p}$, and $Z \in \mathbb R^{n \times p}$ is the noise matrix, $Z_{ij} \overset{iid}\sim \mathcal N \left(0, \frac{1}{n}\right)$, and the signal-to-noise ratio is $\mu \in \mathbb R$. Suppose that we have a deterministic norm-1 vector $v_* \in \mathbb R^{p}$, and a random vector $u \in \mathbb R^{n}$ where $u_i \overset{iid}\sim \mathcal N\left(0, \frac{1}{n}\right)$. 

Our aim ist to estimate the direction of $v_*$, which is norm-1. This boils down to principal component analysis: finding the "important direction" of the data.

### Computing the Maximum Likelihood Estimate

Let's now think about this problem in the high-dimensional regime: where $p$, the dimension of $v_*$, grows with $n$. In classical asymptotics, we know that the MLE is the best estimator, and that it is consistent, i.e. $\hat v_{MLE} \overset{p}\rightarrow v_*$. We also know that, by the asymptotics of the MLE, that $\sqrt{n} (\hat v_{MLE} - v_*) \overset{d}\rightarrow \mathcal N \left(0, I^{-1}(v_*)\right)$.

Let's compute the MLE for this problem and do an asymptotic analysis. 

Consider the rows of $X$: we see from
$$
X = \begin{bmatrix}
\ X_1^\top \ \\\\
\ X_2^\top \ \\\\
\ \vdots \ \\\\
\ X_n^\top \ 
\end{bmatrix}
= \begin{bmatrix}
\ \mu u_1 v_*^\top + Z_1^\top \ \\\\
\ \mu u_2 v_*^\top + Z_2^\top \ \\\\
\ \vdots \ \\\\
\ \mu u_n v_*^\top + Z_n^\top \ 
\end{bmatrix}
$$
that they are independent, as $u$ and $Z$ have independent components. The distribution of each row is $X_i \sim \mathcal N\left(0, \frac{\mu^2}{n} v_* v_*^\top + \frac{1}{n} \mathbf I_p\right)$. Hence, we have that the distribution of $X$ is
$$
\mathbb P_{v_*}(X) = \prod_{i = 1}^n \mathbb P_{v_*}(X_i) = \text{det} \left(\frac{\mu^2}{n} v_* v_*^\top + \frac{1}{n} \mathbf I_p\right)^{-\frac{1}{2}} \exp \left\\{-\frac{n}{2} X_i^\top (\mu^2 v_* v_*^\top + \mathbf I_p)^{-1} X_i\right\\}.
$$
Note that, for the "rank-one updated" matrix $\mu^2 v_* v_*^\top + \mathbf I_p$, the eigenvalues are 1, which has multiplicity $p - 1$, and $\mu^2 + 1$. The eigenvector for $\mu^2 + 1$ has associated eigenvector $v$, and the rest are the $p - 1$ orthogonal directions. Returning to our maximum likelihood calculation, we have
$$
\hat v_{MLE} \overset{(1)}= \arg\min_{\Vert v \Vert_2 = 1} \sum_{i = 1}^n X_i^\top (\mu^2 v_* v_*^\top + \mathbf I_p)^{-1} X_i \overset{(2)}= \arg\min_{\Vert v \Vert_2 = 1} \sum_{i = 1}^n X_i^\top \left(\mathbf I_p - \frac{\mu^2 v v^\top}{\mu^2 v^\top v + 1}\right) X_i \overset{(3)}= \arg\max_{\Vert v \Vert_2 = 1} \sum_{i = 1}^n X_i^\top v v^\top X_i \overset{(4)}= \arg\max_{\Vert v \Vert_2 = 1} v^\top X^\top X v.
$$
Equality $(1)$ we have because maximizing the likelihood means is equivalent to minimizing the negative of the exponential (dropping the positive constant), $(2)$ by using Sherman-Morrison on the matrix inversion, $(3)$ by dropping the constant and because $v$ is norm-1, and the last step $(4)$ by the trace trick.

We see that the last expression is equivalent to PCA: finding the vector that maximizes the quadratic form with the covariance matrix. Hence, under classical asymptotics, the PCA of this particular problem is completely understood through the "standard" MLE asymptotics.

### High-Dimensional MLE Asymptotics

Let's now consider PCA in the high-dimensional case. Suppose that the estimand $v_*$ has dimension $p$ that is allowed to scale with $n$ in some prescribed way. Let's subscript everything with a dummy index $\ell = 1, 2 \ldots $: we have $\{p_1, p_2 \ldots\}$, $\{n_1, n_2 \ldots\}$, even $\{\mu_1, \mu_2 \ldots\}$, and the estimand $\{{v_*}\_1, {v_*}\_2 \ldots\}$. As $\ell \rightarrow \infty$, for each indexed problem, we can use the principal component of $X_\ell^\top X_\ell$ to estimate ${v_*}_\ell$.

If we do that, how should we evaluate the estimator we obtain? For example, we can consider
$$
\Vert \hat v(X_\ell^\top X_\ell) - v_* \Vert_2^2 \overset{p}\rightarrow \fbox{ ? }
$$
and
$$
\langle \hat v, v_* \rangle \overset{p}\rightarrow \fbox{ ? }
$$
to see how well our estimator can do.

For the entire remainder of this section, we will assume, without loss of generality, that $v_* = e_1 = [1, 0, 0 \ldots] \in \mathbb R^p$ for ease of analysis.

Let's do some heuristics first. Without a heuristic analysis, it can be difficult to find a target to work toward. If we do some heuristics first, even if highly non-rigorous, we can get an understanding for what we should mathematically expect.

For some $\hat v \approx v_* = e_1$ (by some metric of comparison), we should have that $X^\top X e_1 \approx \lambda e_1$, because $\hat v$ is an eigenvector, and the eigenvectors should be somewhat robust to small perturbations.

With $v_* = e_1$, we have that
$$
X = \mu u v_*^\top + Z = \begin{bmatrix}
\ \mu u_1 + Z_{11}^\top & Z_{12} & \ldots & Z_{1p} \ \\\\
\ \mu u_1 + Z_{21}^\top & Z_{22} & \ldots & Z_{2p} \ \\\\
\ \vdots & & & \vdots \\\\
\ \mu u_1 + Z_{n1}^\top & Z_{n2} & \ldots & Z_{np} \
\end{bmatrix}.
$$
Denote by $X_1$ the first column of $X$, and $X_{-1}$ the remaining columns, and we have
$$
X^\top X = \begin{bmatrix}
\ X_1 ^\top X_1 & X_1^\top X_{-1} \ \\\\
\ X_{-1}^\top X_1 & X_{-1}^\top X_{-1} \ \\\\
\end{bmatrix} .
$$
We can check if $e_1$ is close to an eigenvector by calculating $X^\top X e_1 = [X_1^\top X_1 \ldots X_p^\top X_1]^\top$ and hoping that it is "close" (by some metric of comparison) to $\lambda e_1$. First computing $X_1^\top X_1$, for example, we have
$$
X_1^\top X_1 = \sum_{i = 1}^n (\mu u_i + z_{i1})^2 \overset{p}\rightarrow \mathbb E[(\mu u_i + z_{i1})^2] = \mu^2 + 1
$$
by the weak law of large nunmbers and independence.

Let's next compute $X_1^\top X_2$: we have
$$
X_1^\top X_2 = \sum_{i = 1}^n (\mu u_i + z_{i1})z_{i2} \rightarrow 0.
$$
Next, we want to know how <i>fast</i> this converges to 0. We do this by concentration, and calculate the variance:
$$
\mathbb E[(X_1^\top X_2)^2] = \mathbb E[X_1^\top X_2 \cdot X_1^\top X_2] = \mathbb E[X_1^\top \mathbb E[X_2 X_2^\top] X_1] = \mathbb E[X_1^\top \frac{1}{n} \mathbf I_n X_1] = \frac{1}{n} \mathbb E[X_1^\top X_1] = \frac{1}{n}(\mu^2 + 1).
$$
Hence, we have $X_1 ^\top X_2 = O_p\left(\sqrt{\frac{\mu^2 + 1}{n}}\right)$ as the rate.

Suppose, for the moment, that $\mu$ is fixed, i.e. that $\mu_\ell = \mu$ for all $\ell$. We then have that
$$
X^\top X e_1 = \begin{bmatrix}
\ X_1^\top X_1 \ \\\\
\ X_1^\top X_2 \ \\\\
\ \vdots \ \\\\
\ X_1^\top X_p \ \\\\
\end{bmatrix} = \begin{bmatrix}
\ 1 + \mu^2 \ \\\\
\ \sqrt{\frac{1 + \mu^2}{n}} \\\\
\ \vdots \ \\\\
\ \sqrt{\frac{1 + \mu^2}{n}}
\end{bmatrix} .
$$
We can see that this should converge to something roughly proportional to $e_1$. 

Recall, however, that we are in a high-dimensional setting. The dimension of $X^\top X e_1$ also grows... and small errors can add up & overwhelm the signal with noise accumulation. Let's check the errors (the accumulation of positive values that are 0 in the estimand):
$$
\sqrt{\sum_{j = 2}^p (X_j^\top X_1)^2} \approx \sqrt{\frac{p (1 + \mu^2)}{n}},
$$
so we want that $1 + \mu^2 \gg \sqrt{\frac{p}{n}(1 + \mu^2)}$. We have $1 + \mu^2 \gg \frac{p}{n}$, and because $\mu$ is a constant, this is equivalent to $n = \omega(p)$.

Heuristically, then, if $\frac{p}{n} \rightarrow 0$, then $\langle \hat v, v_* \rangle \overset{p}\rightarrow 1$, and $\hat \lambda_1 \overset{p}\rightarrow \mu^2 + 1$. On the other hand, if $n = o(p)$, then $\frac{p}{n} \rightarrow \infty$, and it may be true that $\langle \hat v, v_* \rangle \rightarrow 0$: the first principal component would be completely orthogonal to the true estimand $v_* = e_1$, and the true direction is drowned out by the noise.

An edge case, termed endearingly by physicists as the only "nontrivial" case, is when $p = \Theta(n)$. 

### Asymptotics of PCA

We now formalize these heuristic arguments. 

<div class="callout theorem"><span class="label">Theorem: Asymptotics of PCA for Rank-One Spiked Matrix</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Consider the rank-one (single spike) spiked matrix,
$$
X = \mu u v_*^\top + Z,
$$
where $u_i \overset{iid}\sim \mathcal N \left(0, \frac{1}{n}\right)$, $Z_{ij} \overset{iid}\sim \mathcal N \left(0, \frac{1}{n}\right)$, and $v_*$ is a norm-1 deterministic vector we aim to estimate.

Then, in the high-dimensional limit where $p, n \rightarrow \infty$ and $\frac{p}{n} \rightarrow 0$ and $\mu$ is fixed, then
<ol type="i">
  <li>$\lambda_1(X^\top X) \overset{p}\rightarrow \mu^2 + 1$.</li>
  <li>$| \langle \hat v, v_* \rangle | \overset{p}\rightarrow 1$.</li>
</ol>
Hence, $\hat v$ is consistent.
</div>

To prove this result, we need to introduce a bit of background in linear algebra. 

<div class="callout proposition"><span class="label">Proposition: Schur's Formula</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
For a block matrix 
$$
S = \begin{bmatrix}
\ A & B \ \\
\ C & D \
\end{bmatrix},
$$
we have that
$$
\text{det}(S) = \text{det}(A) \cdot \text{det}(D - CA^{-1}B).
$$
</div>

<div class="callout proposition"><span class="label">Proposition: Bounds on the Rayleigh Quotient</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $v \in \mathbb R^n$ be an arbitrary vector and let $A \in \mathbb R^{n \times n}$ be a symmetric matrix. Then
$$
\lambda_{\min}(A) \leq \frac{v^\top A v}{v^\top v } \leq \lambda_{\max}(A).
$$
The quantity $\frac{v^\top A v}{v^\top v}$ is called the <i><strong>Rayleigh quotient</strong></i>.
</div>

<div class="callout definition"><span class="label">Definition: Principal Submatrix</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $A \in \mathbb R^{ \times n}$ be a symmetric matrix, and let $\mathcal I \subseteq \{1, 2 \ldots n\}$ be a subset of indices. Then
$A[I] \triangleq (a_{ij})_{i, j \in \mathcal I}$
is called a <i><strong>principal submatrix</strong></i> of $A$.
</div>

<div class="callout theorem"><span class="label">Theorem: Cauchy's Interlacing Theorem</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $A \in \mathbb R^{n \times n}$ be a symmetric matrix with eigenvalues $\lambda_1(A) \geq \lambda_2(A) \ldots \geq \lambda_n(A)$. Let $B \in \mathbb R^{m \times m}$ be a principal submatrix of $A$. Then, for all $j \in \{1, 2 \ldots m\}$, we have that
$$
\lambda_j(A) \geq \lambda_j(B) \geq \lambda_{j + n - m}(A).
$$
This result is called <i><strong>Cauchy's interlacing theorem</strong></i>.
</div>

A variant of Cauchy's interlacing theorem is the one where the $n \times n$ principal submatrix is obtained by deleting only the $i^{\text{th}}$ row and column.

<div class="callout theorem"><span class="label">Theorem: Cauchy's Interlacing Theorem, Special $n - 1$ Case</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $A \in \mathbb R^{n \times n}$ be a symmetric matrix with eigenvalues $\lambda_1(A) \geq \lambda_2(A) \ldots \geq \lambda_n(A)$. Let $B \in \mathbb R^{n-1 \times n-1}$ be a principal submatrix of $A$ with eigenvalues $\lambda_1(B) \geq \lambda_2(B) \ldots \geq \lambda_{n - 1}(B)$. Then we have that
$$
\lambda_j(A) \geq \lambda_j(B) \geq \lambda_{j +1 }(A)
$$
for $j \in [n - 1]$, or, equivalently, 
$$
\lambda_1(A) \geq \lambda_1(B) \geq \lambda_2(A) \geq \lambda_2(B) \ldots \geq \lambda_{n - 1}(A) \geq \lambda_{n - 1}(B) \geq \lambda_n(A),
$$
where the eigenvalues of $B$ are interlaced with those of $A$.
</div>

<details class="collapsible">
<summary>Proof</summary>
<div class="collapsible__content">
The proof of Cauchy's interlacing theorem requires the <i><strong>Courant-Fischer</i></strong> theorem, which we state here.

<div class="callout theorem"><span class="label">Theorem: Courant-Fischer</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $M \in \mathbb R^{m \times m}$ be a symmetric matrix, and let $\theta_1 \geq \theta_2 \ldots \geq \theta_m$ be its eigenvalues. Then, for $j \in [m]$, where $S_k^m$ is a subspace of $\mathbb R^m$ with dimension $k$, we have that
$$
\theta_j = \max_{S^m_j} \min_{x \in S_j^m\backslash\{0\}} R_M(x) = \min_{S^m_{m - j + 1}} \max_{x \in S^m_{m - j + 1}\backslash\{0\}} R_M(x),
$$
where $R_M(x)$ is the Rayleigh quotient of $M$ and $x$.
</div>
  
</div>
</details>



With the special case, we see that only the first and last (largest and smallest) eigenvalues are not "controlled"" by the eigenvalues of the principal submatrix. This has an important implication for us. As we discussed before, we have that
$$
X^\top X = \begin{bmatrix}
\ X_1 ^\top X_1 & X_1^\top X_{-1} \ \\\\
\ X_{-1}^\top X_1 & X_{-1}^\top X_{-1} \ \\\\
\end{bmatrix},
$$
so, $X_{-1}^\top X_{-1}$ being a principal submatrix of $X^\top X$ of size $p - 1$ and hence being the "bulk" of the whole matrix, should also mostly "determine" the properties of the whole matrix.

Let us denote by $\tilde \lambda_1, \tilde \lambda_2 \ldots \tilde \lambda_{p - 1}$ the eigenvalues of the matrix $X_{-1}^\top X_{-1}$. Then we have by Cauchy's interlacing theorem that
$$
\hat \lambda_1 \geq \tilde \lambda_1 \geq \hat\lambda_2 \geq \ldots \tilde \lambda_{p - 1} \geq \hat\lambda_p\dots
$$
As $\{\tilde \lambda\}$ should be nearly 1, we should have also that $\hat\lambda_2 \ldots \hat\lambda_{p - 1}$ are close to 1 with high probability, i.e. for any $\epsilon > 0$, with probability converging to 1, $\hat\lambda_2 \hat\lambda_3 \ldots \hat\lambda_{p - 1} \in (1 - \epsilon, 1 + \epsilon)$.

Knowing what happens for the $p - 2$ eigenvalues, we now need to think about the last two extremes. Solving $\text{det}(\lambda \mathbf I_p - X^\top X)$ seems a bit difficult, so let's break it down into smaller pieces using Schur's formula above. We have
$$
\text{det}(\lambda \mathbf I_p - X^\top X) = \text{det}(X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1}) \cdot \text{det}\bigg(X_1^\top X_1 - \lambda - X_1^\top X_{-1} (X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1} X_{-1}^\top X_1\bigg) = 0.
$$
From now on, let $\epsilon > 0$ be arbitrary, and consider the event that the $p - 2$ eigenvalues of $X_{-1}^\top X_{-1}$ are between $1 - \epsilon$ and $1 + \epsilon$.

Conditioned on this event, we know that for a solution to the determinant expression above where $\lambda > 1 + \epsilon$, the solution must be to $X_1^\top X_1 - \lambda - X_1^\top X_{-1} (X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1} X_{-1}^\top X_1$ (we omit the determinant because the expression is a scalar).

This is where our heuristic analysis comes in handy. We "know" (guess) that $\lambda_\max(X^\top X)$ should converge to $1 + \mu^2$. Hence, we should expect that 
$$
X_1^\top X_1 - \lambda - X_1^\top X_{-1} (X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1} X_{-1}^\top X_1 = 0
$$ 
should have the solution $\lambda = 1 + \mu^2$. As we calculated before, $X_1^\top X_1 \rightarrow 1 + \mu^2$, so the other term, $X_1^\top X_{-1} (X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1} X_{-1}^\top X_1$, should converge to 0. This is our new immediate goal.

<div class="callout remark"><span class="label">Remark: Consistency of the MLE</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Note that proving that $X_1^\top X_1 - (1 + \mu^2) - X_1^\top X_{-1} (X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1} X_{-1}^\top X_1 \rightarrow 0$ is insufficient to show that $\lambda_{\max}(X^\top X) \rightarrow 1 + \mu^2$. Recall that the proof of consistency for the maximum likelihood estimator required uniform convergence. We would also need to show that a solution exists at $\lambda = 1 + \mu^2$ with high probability using a uniform convergence argument.
</div>

Returning to our goal, we use the proposition above bounding the Rayleigh quotient and obtain that
$$
\lambda_{\min} \bigg((X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1}\bigg) \cdot (X_1^\top X_{-1}) (X_1^\top X_{-1})^\top \leq (X_1^\top X_{-1})\bigg((X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1}\bigg)(X_1^\top X_{-1})^\top \leq \lambda_{\max} \bigg((X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1}\bigg) \cdot (X_1^\top X_{-1})(X_1^\top X_{-1})^\top,
$$ 
which, when we take the negative and the inverse of the matrix, implies that 
$$
\frac{1}{\lambda_{\max}(\lambda \mathbf I - X_{-1}^\top X_{-1})}{(X_1^\top X_{-1})(X_1^\top X_{-1})}^\top \leq - (X_1^\top X_{-1}) (X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1} (X_1^\top X_{-1})^\top \leq \frac{1}{\lambda_{\min}(\lambda \mathbf I - X_{-1}^\top X_{-1})}{(X_1^\top X_{-1})(X_1^\top X_{-1})}^\top.
$$
We see that we should utilize the squeeze theorem and show that the upper and lower bounds converge to 0. This will be split into two parts: showing that $\lambda_\max$ and $\lambda_\min$ are bounded away from 0, and that $X_1^\top X_{-1} X_{-1}^\top X_1$ converges to 0 in probability.

To deal with the former, note that $X_{-1}^\top X_{-1}$ has eigenvalues that, with high probability, are in the interval $(1 - \epsilon, 1 + \epsilon)$. Hence, all the eigenvalues of $\lambda \mathbf I_{p - 1} - X_{-1}^\top X_{-1}$ are, with high probability, in the interval $(\lambda - 1 - \epsilon, \lambda - 1 + \epsilon)$. For any fixed epsilon, we could consider $\lambda > 1 + 2\epsilon$, for example (since we hope that the largest eigenvalue is distinguishably away from 1), and we'd have that the eigenvalues lie, with high probability, in the interval $(\epsilon, 3\epsilon)$. We'd then have that
$$
\frac{1}{3\epsilon} X_1^\top X_{-1} X_{-1}^\top X_1 \ \leq \ (X_1^\top X_{-1})\bigg((X_{-1}^\top X_{-1} - \lambda \mathbf I_{p - 1})^{-1}\bigg)(X_1^\top X_{-1})^\top \ \leq \ \frac{1}{\epsilon} X_{-1} X_{-1}^\top X_1,
$$
where we'd hope to show that $\frac{1}{\epsilon}|X_{-1} X_{-1}^\top X_1| \rightarrow 0$.

If we want to use a $\lambda$ that has no dependence on $\epsilon$, we could also use our conjectured solution $\lambda = 1 + \mu^2$, wwhere recall that $\mu$ is fixed across the dummy index $\ell$ and hence $\mu \gg \epsilon$, the eigenvalues should lie within $(\mu^2 - \epsilon, \mu^2 + \epsilon)$ with high probability, and are hence bounded away from 0.

Let's now analyze the quantity $X_1^\top X_{-1} X_{-1}^\top X_1$. Note that $X_1$ is independent from $X_{-1}$, so using our concentration results, we should have that the expression concentrates around its mean, which is $\mathbb E \left[\left| X_1^\top X_{-1} X_{-1}^\top X_1 - \mathbb E[\text{Tr}()]\right|\right]$

<!-- 

<div class="callout remark"><span class="label">Remark: Remark for Above Example</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Here is a remark.
</div>


<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">

</div>
</details>


<details class="collapsible">
  <summary>Proof</summary>
  <div class="collapsible__content">
  Here is the proof of the above proposition.
  <details class="collapsible">
      <summary>Proof of the sub-proposition.</summary>
      <div class="collapsible__content">
        Here is the sub-proof.
      </div>
    </details>
  </div>
</details>


<div class="callout definition"><span class="label">Definition: Object to Define</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Here is the definition. Here are the list of required properties:
<ol type="i">
  <li>property 1.</li>
  <li>property 2.</li>
  <li>property 3.</li>
</ol>
</div>

 -->