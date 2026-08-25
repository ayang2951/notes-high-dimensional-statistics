## Introduction to the Replica Method

The replica method comes from statistical physics. It is a highly applicable and useful tool for analyzing high dimensional problems. 

Why important?

- It is a standard piece of machinery that is highly applicable.
- The brute-force method is often very messy.
- It is especially useful now because AI tools can help with some of the ugly algebraic calculations.
- Homework 3 idea: do a problem that is not already extremely studied and use the replica method.

> **Added intuition.** The statistical-physics viewpoint is to turn an inference problem into a problem about a **partition function**. The partition function adds up the weight of every possible configuration. Once we understand its logarithm, many quantities we care about can be recovered by differentiating it with respect to carefully chosen parameters.
>
> **Return to handwritten notes below.**

To get a more concrete feel for the method, we use our previously-studied high-dimensional PCA as the example to work through.

## Replica Method on High-Dimensional PCA

We have the rank-one spiked matrix model
$$
X=\mu u_*v_*^\top+Z,
$$
where
$$
u_*\in\mathbb R^n,\quad \Vert u_*\Vert=1,\quad v_*\in\mathbb R^p,\quad \Vert v_*\Vert=1,\quad Z_{ij}\stackrel{\mathrm{iid}}{\sim}N\left(0,\frac{1}{n}\right).
$$
Consider for now the signal strength $\mu$ to be fixed, and let $p$ and $n$ scale in the typical way for proportional high dimensional asymptotics, i.e.
$$
n,p\to\infty,\qquad \frac pn\to\gamma\in(0,\infty).
$$
Our objective is to study the asymptotic correlation between the estimator $\hat v$ and the estimand $v_*$:
$$\langle \hat v,v_*\rangle^2\overset{p}\rightarrow ?$$
The likelihood of the observed matrix $X$ given $u$ and $v$ is, up to constants independent of $u,v$,
$$
f(X\mid u,v)\propto \exp\left\{-\frac n2\operatorname{Tr}[(X-\mu uv^\top)^\top(X-\mu uv^\top)]\right\}\propto \exp\left\{-\frac n2\operatorname{Tr}(X^\top X)+n\mu u^\top Xv\right\}.
$$
Note that because we use priors on $u, v$ such that $\Vert u \Vert_2 =\Vert v \Vert_2=1$, the term $-\frac{1}{2}n\mu^2\Vert uv^\top\Vert_F^2=-\frac{1}{2}n\mu^2$ is constant in $u,v$.

Suppose we have some priors on $u,v$: for example, uniform priors on the spheres in $\mathbb R^n$ and $\mathbb R^p$. Let us write these priors as $d\sigma(u)$ and $d\lambda(v)$. 

Typically, the quantity we want in these Bayesian problems is the posterior distribution $p_n(u, v \mid X)$, which is
$$
p_n(u,v\mid X)=\frac{1}{Z_n}\exp\{n\mu u^\top Xv\}d\sigma(u)d\lambda(v).
$$
The normalizing factor is
$$
Z_n=\int\exp\{n\mu u^\top Xv\}d\sigma(u)d\lambda(v).
$$
This $Z_n$ is called the <i>partition function</i>. 

<div class="callout remark"><span class="label">Remark: Remark for Above Example</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
This partition function, the normalizing factor, is very important: some go so far as to say that calculating the partition function is the <i>core</i> of statistical physics.
</div>

This next step will be quite confusing at first: we shall add parameters $\zeta, h$ and the true estimand $v_*$ into the posterior density function. Adding these parameters, we have the modified function (Arian called this the "drift"?)
$$
\tilde p_n(u,v\mid X)=\frac{1}{\tilde Z_n(h)}\exp\{\zeta\mu n u^\top Xv+nhv_*^\top v\}d\sigma(u)d\lambda(v),
$$
where we would obtain the normalizing factor $\tilde Z_n(h)$
$$
\tilde Z_n(h)=\int\exp\{\zeta\mu n u^\top Xv+nhv_0^\top v\}d\sigma(u)d\lambda(v).
$$
The parameter $\zeta$ controls how much we use the data relative to the prior. The parameter $h$ is not part of the original model; it is added so that differentiating the log partition function gives the overlap with $v_*$.

> **Added intuition.** The $h$-term is like adding a small external field pointing in the true direction $v_0$. If increasing $h$ makes the partition function grow quickly, then the posterior is putting a lot of mass on vectors $v$ aligned with $v_0$. This is why a derivative with respect to $h$ measures alignment.
> 
> (ChatGPT calls this a temperature/data-weight parameter and a source field + source direction, resp.) 
> 
> **Return to handwritten notes below.**

Getting intuition from the Laplace method, our goal will be to focus on the exponent of this function and calculate the normalized (scaled) log of the partition function
$$
\Phi_n(h)=\frac{1}{n}\log\tilde Z_n(h).
$$
Taking $h\to0$ recovers the original posterior. More precisely (assuming the exchangability of integrataion and differentiation),
$$
\frac{\partial\Phi_n(h)}{\partial h}=\frac{\int(v_*^\top v)\exp\{\zeta\mu n u^\top Xv+nhv_*^\top v\}d\sigma(u)d\lambda(v)}{\int\exp\{\zeta\mu n u^\top Xv+nhv_*^\top v\}d\sigma(u)d\lambda(v)}.
$$
Therefore
$$
\lim_{h\to0}\frac{\partial\Phi_n(h)}{\partial h}= \lim_{h\to 0} \frac{1}{n} \cdot \frac{\frac{\partial}{\partial h} \tilde Z_n(h)}{\tilde Z_n(h)} = v_*^\top \frac{1}{\tilde Z_n(0)} \cdot \int v \exp\{\zeta n\mu u^\top Xv\} d\sigma(u) d\lambda(v) = \left\langle v_*,\mathbb E_{\tilde p_n(\cdot\mid X,h=0)}[v]\right\rangle.
$$
Hence the first intuition is that we should modify the partition function so that derivatives of its log give quantities we want.

> **Added clarification.** The derivative above gives a one-replica overlap with $v_0$. If the target is a squared overlap, one often either differentiates a slightly different source term, or uses two replicas so that a product of overlaps appears. The main idea is the same: add a parameter, compute the log partition function as a function of that parameter, then differentiate.
>
> **Return to handwritten notes below.**

### Replica Trick #1

The main asymptotic target is
$$
\lim_{n\to\infty}\tilde\Phi_n(h)=\lim_{n\to\infty}\frac{1}{2}\log\tilde Z_n(h).
$$
Our hope is that we will be able to use concentration here via self-averaging of the function: that
$$
\frac{1}{n}\log\tilde Z_n(h)-\frac{1}{n}\mathbb E\log\tilde Z_n(h)\to0.
$$
If this holds, then instead of computing the quantity direction, we can instead calculate the simpler
$$
\frac{1}{n}\mathbb E[\log\tilde Z_n(h)].
$$
This logarithm, however, is bothersome, and we cannot compute the quantity easily because of its expression. The next replica trick is to omit the necessity of computing this expectation.

<!-- (We also want an exact computation at least a lot of the time, so Jensen-type bounds are not enough.) -->

### Replica Trick #2

A second trick used in the replica method is the identity
$$
\mathbb E[\log\tilde Z_n(h)]=\lim_{r\to0}\frac{\partial}{\partial r}\mathbb E[\tilde Z_n(h)^r],
$$
so the next goal is to calculate $\mathbb E[\tilde Z_n(h)^r]$ for $r$ in a neighborhood of $0$.

### Replica Trick #3

The question now is how to compute this expectation for small $r$.

The main heuristic flow here will be:

1. Calculate $\mathbb E[\tilde Z_n(h)]$, $\mathbb E[\tilde Z_n(h)^2]$, and so on.
2. More generally, calculate $\mathbb E[\tilde Z_n(h)^r]$ for $r\in\mathbb Z_+$.
3. Call the resulting expression $f(r)$, and then pretend that the domain of $f$ includes $r\in(0,1)$. In other words, pretend that the integer formula also applies to values between $0$ and $1$.

Tricks 1 and 2 can often be made rigorous. Trick 3 is the genuinely heuristic part. 

Remark: when the replica method fails, it is not often because of trick #3 anyway, it's because of something else (more discussion later).

> **Added intuition.** The word **replica** comes from the integer moment $\tilde Z_n^r$. When $r$ is a positive integer, multiplying $r$ copies of the same integral creates $r$ copies of the variables $(u,v)$. These copies are the replicas. The data matrix $X$ is shared across the copies, which is exactly why the replicas interact after we average over the randomness.
>
> **Return to handwritten notes below.**

From now on, until almost the end of the replica analysis, assume $r\in\mathbb Z_+$.

### Creating the Replicas

Let us return to our high-dimensional PCA problem. Using the replica method on the problem, we invoke trick #3 and proceed to calculate the $r^{\text{th}}$ moment of $\tilde Z_n(h)$:
$$
\mathbb E[\tilde Z_n(h)^r]=\mathbb E\left[\left(\int\exp\{\zeta\mu n u^\top Xv+nhv_*^\top v\}d\sigma(u)d\lambda(v)\right)^r \ \right].
$$
We shall do something clever here, where the replica method gets its name: instead of the power $r$, we shall create $r$ copies of the densities (split the $r$ powers of $u$ and $v$ into $r$ separate versions, keeping $X$ the same), also called <i>replicas</i>:
$$
\mathbb E[\tilde Z_n(h)^r]=\mathbb E\left[\prod_{i=1}^r\int\exp\{\zeta\mu n u_i^\top Xv_i+nhv_*^\top v_i\}d\sigma(u_i)d\lambda(v_i)\right].
$$
Now return to the model of the observable $X=\mu u_*v_*^\top+Z$ to remove most of the randomness from the system: if we can compute out the expectation over $Z$, we are left only with the rank-one matrix involving $u_*$ and $v_*$. 
$$
\mathbb E[\tilde Z_n(h)^r]=\mathbb E\left[\prod_{i=1}^r\int\exp\{\zeta\mu n u_i^\top Zv_i+\zeta\mu^2 n(u_i^\top u_*)(v_*^\top v_i)+nhv_*^\top v_i\}d\sigma(u_i)d\lambda(v_i)\right].
$$
We use the tower law and Tonelli's theorem to simplify:
$$
\mathbb E[\tilde Z_n(h)^r] = \mathbb E_{u_*, v_*, Z} \left[\int \prod_{i = 1}^r \exp\{\zeta n \mu u_i^\top Z v_i\}\prod_{i = 1}^r \exp\{\zeta n\mu^2 (u_i^\top u_*)(v_*^\top v_i) + n h v_*^\top v_i\}d\sigma(u_i) d\lambda(v_i)\right] = \mathbb E_{u_*, v_*} \left[\int \mathbb E_Z\left[\prod_{i = 1}^r \exp \{\zeta n \mu u_i^\top Z v_i\}\right] \prod_{i = 1}^r \exp\{\zeta n \mu^2 u_i^\top u_* v_*^\top v_i + n h v_*^\top v_i\}d\sigma(u_i) d\lambda(v_i)\right].
$$
From this expression, we can identify our next objective; we want to compute the expectation over $Z$ first:
$$
\mathbb E_Z\left[\prod_{i=1}^r\exp\{\zeta\mu n u_i^\top Zv_i\}\right].
$$
Let us define the quantities
$$
w_i=u_i^\top Zv_i,\qquad w=(w_1,\ldots,w_r)^\top.
$$
We hence have that $w \sim \mathcal N(0, \Sigma)$ is Gaussian with covariance matrix
$$
\Sigma_{ij} = \mathbb E_Z[(u_i^\top Zv_i)(u_j^\top Zv_j)]=\frac{1}{n}(u_i^\top u_j)(v_i^\top v_j).
$$
From the MGF of a zero-mean Gaussian, we have that the expectation is therefore
$$
\mathbb E_Z[\exp\{\zeta\mu n\mathbf 1^\top w\}]=\exp\left\{\frac{1}{2}\zeta^2\mu^2 n\sum_{i=1}^r\sum_{j=1}^r(u_i^\top u_j)(v_i^\top v_j)\right\}.
$$
These products, the terms $u_i^\top u_j$ and $v_i^\top v_j$, are called <i>overlaps</i>, and getting to the overlaps is one of the main intermediate goals. In fact, with replica method applications, what is often necessary is to simplify <i>until</i> the overlaps appear. After that, the rest will be to investigate these overlaps.


### Overlap Matrices

At this point, let us redefine some notation: rename the true signals $u_*$ and $v_*$ as $u_0$ and $v_0$, respectively, to make them more comparable to one of the replicas (but keep in mind that $u_0,v_0$ are philosophically distinct from the genuine replicas $u_i,v_i$, even though they are included as index $0$ in the same overlap matrices).

What we have now, with the subtitution in notation above and simplifying the expectation over $Z$, integrals numbering $r + 1$:
$$
\int_{u_0, v_0} \int_{u_1, v_1} \ldots \int_{u_r, v_r} \exp\{n\mu^2 \zeta^2 \mathbf 1^\top \Sigma \mathbf 1\} \prod_{i = 1}^r \exp\{n\mu^2 u_i^\top u_0 v_0^\top v_i + nhv_0^\top v_i\} d\sigma(u_r)d\lambda(v_r) \ldots d\sigma(u_1)d\lambda(v_1)d\sigma(u_0)d\lambda(v_0)
$$
Define the overlap matrices
$$
Q_{ij}=v_i^\top v_j,\qquad R_{ij}=u_i^\top u_j,\qquad 0\le i,j\le r.
$$
Let $f_n(Q)$ and $g_n(R)$ denote the induced distributions of the overlap matrices. Everything is now written in terms of overlaps:
$$
\mathbb E[\tilde Z_n(h)^r]=\iint\exp\left\{n\left[\frac{1}{2}\zeta^2\mu^2\sum_{i=1}^r\sum_{j=1}^r Q_{ij}R_{ij}+\zeta\mu^2\sum_{i=1}^rQ_{0i}R_{0i}+h\sum_{i=1}^rQ_{0i}\right]\right\}df_p(Q)dg_n(R).
$$
Let us take a brief look back to see what these 3 replica tricks have gotten us so far...

> **Added intuition.** This is the key compression step. The original integral ranges over high-dimensional spheres, but after averaging over the Gaussian noise, the exponent only depends on the pairwise angles between replicas. Those pairwise angles are exactly the entries of $Q$ and $R$.
>
> **Return to handwritten notes below.**

Now that we have simplified our problem to a double integration with $r \times r$ matrices, we shall use our large deviation and Laplace method results to simplify the densities of the matrices $Q$ and $R$ and make them workable.

To simplify this integral, we want a large-deviation principle for the overlap distributions. Informally, suppose
$$
f_n(Q)\asymp\exp\{-pI_v(Q)\},\qquad g_n(R)\asymp\exp\{-nI_u(R)\}.
$$
Recall Varadhan's theorem in the form needed here:
$$
\frac{1}{n}\log\int e^{nG(x)}\,d\mu_n(x)\to\sup_x\{G(x)-I(x)\}.
$$
Applying this with $\frac{p}{n}\to\gamma$ gives the variational form
$$
\psi(r,h)=\lim_{n\to\infty}\frac{1}{n}\log\mathbb E[\tilde Z_n(h)^r]=\sup_{Q,R}\left\{\frac{1}{2}\zeta^2\mu^2\sum_{a=1}^r\sum_{b=1}^r Q_{ab}R_{ab}+\zeta\mu^2\sum_{a=1}^rQ_{0a}R_{0a}+h\sum_{a=1}^rQ_{0a}-\gamma I_v(Q)-I_u(R)\right\}.
$$



> **Added clarification.** 
> 
> Remark: note that the $\gamma$ term appears because the rate functions of $Q$ and $R$ are the same; the only difference is in the dimension.
>
> After computing $\psi(r,h)$ for integer $r$, the replica prediction for the limiting averaged free energy is obtained by differentiating at $r=0$. In shorthand,
> 
>
> $$\lim_{n\to\infty}\frac{1}{n}\mathbb E[\log\tilde Z_n(h)]\approx \left.\frac{\partial}{\partial r}\psi(r,h)\right|_{r=0}.$$
>
> **Return to handwritten notes below.**

To use this, we need to show that the overlap matrices satisfy a large-deviation principle and identify the rate functions $I_v$ and $I_u$.

### Large deviations for overlaps

For the $v$-overlaps, represent uniform sphere vectors by normalized Gaussians:
$$
v_a=\frac{z_a}{\Vert z_a\Vert_2},\qquad z_a\sim N(0,I_p).
$$
Then
$$
v_a^\top v_b=\frac{z_a^\top z_b}{\Vert z_a\Vert_2\Vert z_b\Vert_2}.
$$
Start with the unnormalized Gram matrix
$$
\tilde Q=\frac{1}{p}\begin{bmatrix}
\ z_0^\top z_0 & z_0^\top z_1 & \cdots&z_0^\top z_r \
\\ 
\ z_1^\top z_0 & z_1^\top z_1 & \cdots & z_1^\top z_r \
\\ 
\ \vdots & \vdots &\ddots & \vdots \
\\ 
\ z_r^\top z_0 & z_r^\top z_1 & \cdots & z_r^\top z_r \
\end{bmatrix}.
$$
Equivalently, if $b_j=(z_{0j},z_{1j},\ldots,z_{rj})^\top\in\mathbb R^{r+1}$, then
$$
\tilde Q=\frac{1}{p}\sum_{j=1}^p b_jb_j^\top,\qquad b_j\stackrel{\mathrm{iid}}{\sim}N(0,I_{r+1}).
$$
So $\tilde Q$ is a sample covariance matrix. This puts us in the setting of Cramer's theorem.

Recall Cramer's theorem: if $X_1,\ldots,X_n$ are iid and $\Lambda(t)=\log\mathbb E[e^{tX_1}]$, then
$$
I(z)=\sup_t\{tz-\Lambda(t)\}.
$$
For the vector version, if $X_i\in\mathbb R^d$, then
$$
I(z)=\sup_{t\in\mathbb R^d}\{\langle z,t\rangle-\Lambda(t)\},\qquad \Lambda(t)=\log\mathbb E[e^{\langle t,X_1\rangle}].
$$
For the matrix version relevant here, the rate function for $\tilde Q$ is
$$
I(\tilde Q)=\sup_{\Lambda}\left\{\operatorname{Tr}(\Lambda\tilde Q)-\log\mathbb E\left[e^{\operatorname{Tr}(\Lambda bb^\top)}\right]\right\},\qquad b\sim N(0,I_{r+1}).
$$
Issue: $\tilde Q$ is symmetric, so we should use only the upper-triangular part, or equivalently optimize over symmetric matrices without double-counting off-diagonal entries.

> **Added clarification.** The clean way to avoid double-counting is to write the Gaussian quadratic form with a symmetric matrix $\Gamma$ and a factor $1/2$. This is the same calculation as the upper-triangular version, just written in matrix notation.
>
> **Return to handwritten notes below.**

Let $d=r+1$. Using a symmetric parameter $\Gamma$,
$$
I(\tilde Q)=\sup_{\Gamma=\Gamma^\top,\,I-\Gamma\succ0}\left\{\frac{1}{2}\operatorname{Tr}(\Gamma\tilde Q)-\log\mathbb E\left[\exp\left(\frac{1}{2} b^\top\Gamma b\right)\right]\right\}.
$$
The Gaussian integral gives
$$
\mathbb E\left[\exp\left(\frac{1}{2} b^\top\Gamma b\right)\right]=\det(I-\Gamma)^{-1/2}.
$$
Therefore, we have that
$$
I(\tilde Q)=\sup_{\Gamma=\Gamma^\top,\,I-\Gamma\succ0}\left\{\frac{1}{2}\operatorname{Tr}(\Gamma\tilde Q)+\frac{1}{2}\log\det(I-\Gamma)\right\}.
$$
Differentiate the objective with respect to $\Gamma$:
$$
\nabla_\Gamma\left[\frac{1}{2}\operatorname{Tr}(\Gamma\tilde Q)+\frac{1}{2}\log\det(I-\Gamma)\right]=\frac{1}{2}\tilde Q-\frac{1}{2}(I-\Gamma)^{-1}.
$$
Setting this equal to zero gives
$$
\Gamma^*=I-\tilde Q^{-1}.
$$
Plugging back in,
$$
I(\tilde Q)=\frac{1}{2}\operatorname{Tr}(\tilde Q)-\frac{1}{2}d-\frac{1}{2}\log\det(\tilde Q),\qquad d=r+1.
$$

> **Added cleanup.** The standard Wishart rate has the sign $-\frac{1}{2}\log\det(\tilde Q)$. With a $+\frac{1}{2}\log\det(\tilde Q)$ sign, the expression would not be a nonnegative rate function. At $\tilde Q=I$, the rate is $0$, as it should be.
>
> **Return to handwritten notes below.**

For the normalized sphere-overlap matrix $Q$, the diagonal entries are fixed at $1$. On this constraint, the trace term cancels the constant, so the rate function is, up to an additive constant,
$$
I_v(Q)=-\frac{1}{2}\log\det(Q),\qquad Q\succeq0,\quad Q_{aa}=1.
$$
The same calculation applies to $R$:
$$
I_u(R)=-\frac{1}{2}\log\det(R),\qquad R\succeq0,\quad R_{aa}=1.
$$
So the overlap variational problem can be written as
$$
\psi(r,h)=\sup_{Q,R}\left\{\frac{1}{2}\zeta^2\mu^2\sum_{a=1}^r\sum_{b=1}^r Q_{ab}R_{ab}+\zeta\mu^2\sum_{a=1}^rQ_{0a}R_{0a}+h\sum_{a=1}^rQ_{0a}+\frac{\gamma}{2}\log\det Q+\frac{1}{2}\log\det R\right\}.
$$
The remaining replica-analysis step would be to evaluate this finite-dimensional variational problem for integer $r$, analytically continue the answer, and then take $r\to0$.

> **Added intuition.** The log-determinant terms are entropy terms. They penalize overlap matrices that are too rigid or too aligned. The energy terms are the terms involving $\mu,\zeta,h$, which reward overlaps that fit the spike and the source field. The variational problem is therefore an energy-entropy tradeoff.
>
> **Return to handwritten notes below.**

