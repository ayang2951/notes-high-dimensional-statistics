_Overview: This section covers basic random matrix theory, including important tools for analysis and a canonical example._

## An Introduction to Random Matrix Theory 

A common and critically important issue in high dimensional statistics is the asymptotic distribution of functions of random matrices, such as its eigenvalues. 

In this section, we introduce some motivating examples and some critical tools, most notably the <i>Stieltjes transform</i>.

We begin with a relatively simple case, a random matrix $X \in \mathbb R^{n \times p}$ where its entries are independently and identically distributed $\mathcal N(0, \frac{1}{n})$. Let us consider the eigenvalues $\lambda_1, \lambda_2 \ldots \lambda_p$ of $X^\top X$, in the case where $p, n \rightarrow \infty$, and $\frac{p}{n} \rightarrow \gamma$. We are interested in the empirical spectral distribution of $X^\top X$, i.e.
$$
\mu_n \triangleq \frac{1}{p} \sum_{i = 1}^p \delta_{\lambda_i}.
$$
We are interested in the asymptotics of this object: what does $\mu_n$ converge to, and in what sense? For example, we could have
<ol type="i">
  <li>Weak convergence, i.e. $\mu_n \overset{w}\rightarrow \mu$, by definition $\mathbb E_{\mu_n}[f(X)] \rightarrow \mathbb E_\mu[f(X)]$ for all bounded and continuous functions $f$,</li>
  <li>Weak convergence <i>in probability</i>, i.e. $\mathbb E_{\mu_n}[f(X)] \overset{p}\rightarrow \mathbb E_{\mu}[f(X)]$, as $\mu_n$ themselves are random measures,</li>
  <li>Weak convergence almost surely, the anology of the above with almost sure convergence instead of in probability.</li>
</ol>

There are a few standard, general tools used for proving weak convergence, including:
<ol type="i">
  <li>Using the definition.</li>
  <li>Using the method of moments, if the distribution comes from a suitable class (finite moments that uniquely characterize the distribution): proving that $\mathbb E_{\mu_n}[X^t] = \mathbb E_\mu[X^t]$.</li>
  <li>Using characteristic functions: proving that $\mathbb E_{\mu_n}[e^{itx}] \rightarrow \mathbb E_\mu[e^{itx}]$ for all $t$. This is a neat, fast tool when there's independence in the problem&mdash;but this is precisely what does <i>not</i> exist for random matrices.</li>
</ol>

These tools happen not to be the best one for our problem, so we introduce a powerful new one: the Stieltjes transform.

## The Stieltjes Transform

<div class="callout definition"><span class="label">Definition: Stieltjes Transformn</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\mu$ denote a measure on $\mathbb R$. The <strong><i>Stieltjes transform</i></strong> of $\mu$ is given by
$$
\mathcal S_\mu(z) = \int \frac{1}{t - z} d\mu(t)
$$
for all $z \in \mathbb C \backslash \mathbb R$.
</div>

We now introduce a few properties of the Stieltjes transform, of which we highlight two important ones here to pay particular attention to:
<ul>
  <li>Given the Stieltjes transform $\mathcal S_\mu(z)$, we can recover $\mu$.</li>
  <li>If the Stieltjes transform converges for a sequence of measures for all $z$, so does the sequence of measures, weakly, i.e. if $\mathcal S_{\mu_n}(z) \rightarrow \mathcal S_\mu(z)$ for all $z \in \mathbb C \backslash \mathbb R$, then $\mu_n \overset{w}\rightarrow \mu$.</li>
</ul>

With these goals in mind, we proceed to formally presenting some results on the Stieltjes transform that will be useful to us later on.

<div class="callout proposition"><span class="label">Proposition: Stieltjes Transform for Compactly Supported Probability Measures</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\mu$ be a compactly supported probability measure, i.e. let there exist an $M$ such that $\mu([-M, M]) = 1$.

Then for all $z \in \mathbb C \backslash \mathbb R$ where $|z| > M$, 
$$
\mathcal S_\mu(z) = - \sum_{m = 0}^\infty \frac{m_n}{z^{n + 1}},
$$
where $m_n \triangleq \int t^n d\mu(t)$.
</div>


<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
We rewrite the integrand for the Stieltjes transform and utilize that $\mu$ is compactly supported:
$$
\mathcal S_\mu(z) = \int \frac{1}{t - z} d\mu(t) = -\frac{1}{z} \int \frac{1}{1 - \frac{t}{z}} d\mu(t) =  -\frac{1}{z} \int_{-M}^M \frac{1}{1 - \frac{t}{z}} d\mu(t).
$$
We have that $|z| > m$ and that, on the domain of integration, $|t| \leq M$, so we know that $|\frac{t}{z}| < 1$. Hence, we may write the expression above as the infinite geometric series, which is absolutely summable:
$$
\mathcal S_\mu(z) = -\frac{1}{z} \int_{-M}^M \sum_{n = 0}^\infty \left(\frac{t}{z}\right)^n d\mu(t) = - \int_{-M}^M \sum_{n = 0}^\infty \left(\frac{t^n}{z^{n + 1}}\right) d\mu(t).
$$
Because the sum is finite, we can use dominated convergence to obtain
$$
\mathcal S_\mu(z) = - \sum_{n = 0}^\infty \frac{1}{z^{n + 1}} \int_{-M}^M t^n d\mu(t) = - \sum_{n = 0}^\infty \frac{1}{z^{n + 1}} \int t^n d\mu(t) = - \sum_{n = 0}^\infty \frac{m_n}{z^{n + 1}}.
$$
</div>
</details>

The following theorem essentially states that, given the Stieltjes transform of a measure $\mu$, we can recover $\mu$.

<div class="callout theorem"><span class="label">Theorem: Stieltjes Inversion Formula</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\mu$ be a a measure on $\mathbb R$, and let $a < b$ be two continuity points of $\mu$. Then
$$
\lim_{\epsilon \rightarrow 0^+} \frac{1}{\pi} \int_{a}^b \text{Im} \left(\mathcal S_\mu(x + i\epsilon)\right) d\lambda(x) = \mu((a, b)),
$$
where $\lambda$ is the Lebesgue measure.
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
We write out the Stieltjes transform:
$$
\lim_{\epsilon \rightarrow 0^+} \frac{1}{\pi} \int_{a}^b \text{Im} \left(\int \frac{1}{t - (x + i \epsilon)}  d\mu(t)\right) d\lambda(x) = \lim_{\epsilon \rightarrow 0^+} \frac{1}{\pi} \int_{a}^b \text{Im} \left( \int \frac{t - x + i \epsilon}{(t - x)^2 + \epsilon^2}d\mu(t)\right)  d\lambda(x) = \lim_{\epsilon \rightarrow 0^+} \frac{1}{\pi} \int_{a}^b \int \left(\frac{\epsilon}{(t - x)^2 + \epsilon^2}\right) d\mu(t) d\lambda(x).
$$
We recognize the integrand as being the derivative of arctan, and, using Tonelli and bounded convergence, obtain
$$
\lim_{\epsilon \rightarrow 0^+} \frac{1}{\pi} \int_{a}^b \int_{\mathbb R} \frac{1}{\epsilon} \left({\left(\frac{x - t}{\epsilon}\right)^2 + 1}\right)^{-1} d\mu(t) d\lambda(x) = \frac{1}{\pi} \int_{\mathbb R} \lim_{\epsilon \rightarrow 0^+}  \left[\tan^{-1}\left(\frac{x - t}{\epsilon}\right) \bigg |_a^b \right]d\lambda(x) d\mu(t).
$$
We can now split the integral across the domain $\mathbb R$, obtaining integrals over $(- \infty, a)$, $\{a\}$, $(a, b)$, $\{b\}$, and $(b, \infty)$. Given that $a, b$ are continuity points, only three terms remain. We utilize that $x < a$, $x \in (a, b)$, and $x > b$ on each domain of integration to obtain 
$$
\frac{1}{\pi} \int_{(-\infty, a)} \lim_{\epsilon \rightarrow 0^+}  \left[\tan^{-1}\left(\frac{x - t}{\epsilon}\right) \bigg |_a^b \right] d\mu(t) + \frac{1}{\pi} \int_{(a, b)} \lim_{\epsilon \rightarrow 0^+}  \left[\tan^{-1}\left(\frac{x - t}{\epsilon}\right) \bigg |_a^b \right]d\mu(t) + \frac{1}{\pi} \int_{(b, \infty)} \lim_{\epsilon \rightarrow 0^+}  \left[\tan^{-1}\left(\frac{x - t}{\epsilon}\right) \bigg |_a^b \right]d\mu(t) = \frac{1}{\pi} \int_{(a, b)} (\frac{\pi}{2} - \left(- \frac{\pi}{2}\right)) d\mu(t) = \mu((a, b)).
$$
</div>
</details>

The next theorem is the most important of this set of results. It gives us equivalence in convergence of a sequence of measures and their Stieltjes transforms.

<div class="callout theorem"><span class="label">Theorem: Continuity of Stieltjes Transforms</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\mu_n$ be a sequence of measures and let $\mathcal S_{\mu_n}(z)$ be their Stieltjes transforms. Then
<ol type="i">
  <li>For $\mu_n \overset{w}\rightarrow \mu$, $\mathcal S_{\mu_n}(z) \rightarrow \mathcal S_\mu(z)$ for all $z \in \mathbb C \backslash \mathbb R$.</li>
  <li>For $\mathcal S_{\mu_n}(z) \rightarrow \mathcal S_\mu(z)$ for all $z \in \mathbb C \backslash \mathbb R$, $\mu_n \overset{w}\rightarrow \mu$.</li>
  <li>If $\mu_n$ are random, then if $\mathcal S_{\mu_n}(z) \overset{p}\rightarrow \mathcal S_\mu(z)$, we have $\mu_n \overset{w} \rightarrow \mu$ in probability.</li>
</ol>
</div>

We can now use these results on the Stieltjes transform to define the weak limit of the empirical spectral distribution of our $X^\top X$ where $X$ has iid Gaussian entries.

<div class="callout theorem"><span class="label">Theorem: The Marçenko-Pastur Law</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X \in \mathbb R^{n \times p}$ be a random matrix where $X_{ij} \overset{iid}\sim \mathcal N(0, \frac{1}{n})$. Let $\mu_n$ denote the empirical spectral distribution of $X^\top X$, i.e. $\mu_n \triangleq \frac{1}{p} \sum_{i = 1}^p \delta_{\lambda_i}$, where $\lambda_i$ are the eigenvalues of $X^\top X$.

Suppose for simplicity that $p < n$, i.e. the matrix $X^\top X$ is full rank. Then $\mu_n$ converges weakly in probability to the density
$$
f_{\gamma}(x) = \frac{1}{2\pi \gamma x} \sqrt{(\gamma_+ - x) (x - \gamma_-)} \cdot \mathbb I_{[\gamma_-, \gamma_+]},
$$
where, recall, $\frac{p}{n} \rightarrow \gamma$ is the aspect ratio, and $\lambda_{-} \triangleq (1 - \sqrt{\gamma})^2$, $\lambda_{+} \triangleq (1 + \sqrt{\gamma})^2$.
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
The proof is a bit long, and we'll provide heuristics on the way to guide our understanding and strategy.

For the matrix $X^\top X - z \cdot \mathbf I_p$, its eigenvalues are $\lambda_i - z$, and $(X^\top X - z\cdot \mathbf I_p)^{-1}$ has eigenvalues $\frac{1}{\lambda_i - z}$. We may compute the Stieltjes transform of $\mu_n$:
$$
\mathcal S_{\mu_n}(z) = \int \frac{1}{t - z} d\mu(t) = \frac{1}{p} \sum_{i = 1}^n \frac{1}{\lambda_i - z} = \frac{1}{p} \text{Tr}((X^\top X - z \cdot \mathbf I_p)^{-1}).
$$
We now focus on the matrix $(X^\top X - z \cdot \mathbf I_p)^{-1}$, and, for simplicity, consider only the first element. To facilitate this analysis, we write the matrix in a block matrix form.

$$&#9671;&#9672;&#9671;$$
An interlude: a useful linear algebra result. For a matrix
$$
M = \begin{bmatrix}
\ A & B \ \\
\ C & D \
\end{bmatrix}
$$
and defining $S \triangleq D - CA^{-1} B$, the inverse of $M$ takes the form
$$
M^{-1} = \begin{bmatrix}
\ A^{-1} + A^{-1} B S^{-1} C A^{-1} & -A^{-1} B S^{-1} \ \\
\ -S^{-1} CA^{-1} & S^{-1} \
\end{bmatrix}
$$
Now, utilizing the Woodbury identity on the top left block of $M$, we have that
$$
A^{-1} + A^{-1} B (D - C A^{-1} B)^{-1} C A^{-1} = (A + B D^{-1} C)^{-1}.
$$
That concludes our linear algebra interlude.
$$&#9671;&#9672;&#9671;$$

Returning to our matrix of interest and defining $x_1$ as the first column of $X$ and $X_{-1}$ to be the remaining $p - 1$ columns, we write it in the block form, isolating the first component:
$$
X^\top X - z \cdot \mathbf I_p = \begin{bmatrix}
\ x_1^\top x_1 - z & x_1^\top X_{-1} \ \\
\ X_{-1}^\top x_1 & X_{-1}^\top X_{-1} - z \cdot \mathbf I_{p - 1} \ \
\end{bmatrix}.
$$
Utilizing the result above about the block inversion, we have that
$$
(X^\top X - z \cdot \mathbf I_p)^{-1}_{11} = \bigg((x_1^\top x_1 - z) - x_1^\top X_{-1}(X_{-1}^\top X_{-1} - z \cdot \mathbf I_{p - 1})^{-1} X_{-1}^\top x_1\bigg)^{-1}.
$$
On this quantity, we can either apply concentration directly, or we can attempt to simplify it via linear algebra first. If we choose to do the latter (which is probably nicer), we can use the push-through identity, which gives us $-x_1^\top X_{-1} (X_{-1}^\top X_{-1} - z \cdot \mathbf I_{p - 1})^{-1}X_{-1}^\top x = x_1^\top x_1 + z \cdot x_1^\top (X_{-1} X_{-1}^\top - z \cdot \mathbf I_{n})^{-1} x_1$, and hence obtain
$$
(X^\top X - z \cdot \mathbf I_p)^{-1}_{11} = (-z - z x_1^\top (X_{-1} X_{-1}^\top - z \cdot \mathbf I_{p - 1})^{-1} x_1).
$$
Observe the term $x_1^\top (X_{-1}X_{-1}^\top - z \cdot \mathbf I)^{-1} x_1$. Because $x_1$ and $X_{-1}$ are independent, we have
$$
x_1^\top (X_{-1}X_{-1}^\top - z \cdot \mathbf I)^{-1} x_1 = \sum_{i = 1}^n (X_{-1} X_{-1}^\top - z \cdot \mathbf I)^{-1}_{ii} \cdot x_i^2 + \sum_{i \neq j} (X_{-1} X_{-1}^\top - z \cdot \mathbf I)^{-1}_{ij} \cdot x_i \cdot x_j
$$
Here, we can develop some heuristics for the resultant quantity, $\frac{1}{n} \text{Tr}(X_{-1} X_{-1}^\top - z \cdot \mathbf I)^{-1}$. We see that, in the limit, the quantity $(X_{-1} X_{-1}^\top - z \cdot \mathbf I) $ shouldn't be very different than $(XX^\top - z\cdot \mathbf I)$, and specifically, the trace is missing only one eigenvalue. Hence, it's reasonable to expect that $\text{Tr}(X_{-1} X_{-1}^\top - z \cdot \mathbf I)^{-1} \approx \text{Tr}(X X^\top - z \cdot \mathbf I)^{-1}$. If this is true, then we can expect that
$$
(X^\top X - z \cdot \mathbf I)_{11}^{-1} = (-z -z x_1^\top (X_{-1} X_{-1}^\top - z \cdot \mathbf I)^{-1} x_1)^{-1} \approx (-z - \frac{z}{n} \text{Tr}(X_{-1} X_{-1}^\top - z \cdot \mathbf I)^{-1}) \approx (-z - \frac{z}{n} \text{Tr}(X X^\top - z \cdot \mathbf I)^{-1}).
$$
What does this give us? Notice that when we apply our heuristic argument, the first element doesn't depend on any specific column! Hence, we should be able to say that 
$$
(X^\top X - z \cdot \mathbf I)_{11}^{-1} \approx (X^\top X - z \cdot \mathbf I)_{22}^{-1} \ldots \approx (X^\top X - z \cdot \mathbf I)_{pp}^{-1}
$$
and jump directly from a calculation on the first element of the diagonal to the entire trace. 

<details class="collapsible">
<summary>Proof of the Heuristic Argument.</summary>
<div class="collapsible__content">
Before we proceed further, we may take the time to prove our claim approximating $X_{-1} X_{-1}^\top - z \cdot \mathbf I$ with $X X^\top - z \cdot \mathbf I$. To that end, we'll apply the linear algebra result introduced (and proved) below.

Let $B$ be a Hermitian matrix, i.e. for $\mathsf{H}$ denoting the conjugate transpose, we have that $B = B^{\mathsf{H}}$. Then for $z \in \mathbb C \backslash \mathbb R$ and some conformable real vector $v$,
$$
\bigg| \text{Tr}(B + v v^\top - z \cdot \mathbf I)^{-1} - \text{Tr}(B - z \cdot \mathbf I)^{-1} \bigg| \leq \text{Im}(z)^{-1}.
$$
The quantity $v v^\top$ represents a rank-1 update to the matrix $B$. 

To prove this, we use the Woodbury identity on the first term of the left hand side and subsequently the trace trick on the denominator to get
$$
\bigg | \text{Tr}(B - z \cdot \mathbf I)^{-1} - \text{Tr} \left(\frac{(B - z \cdot \mathbf I)^{-1} v v^\top (B - z \cdot \mathbf I)^{-1}}{1 + v^\top (B - z \cdot \mathbf I)^{-1} v}\right) - \text{Tr}(B - z \cdot \mathbf I)^{-1}\bigg | = \left| \frac{v^\top (B - z \cdot \mathbf I)^{-2} v
}{1 + v^\top (B - z \cdot \mathbf I)^{-1} v}\right|.
$$
And since $B$ is Hermitian, it can be diagonalized by a unitary matrix $Q$ into $B = Q \Lambda Q^{\mathsf{H}}$. Writing the numerator of the above in terms of $Q$, we have
$$
v^\top (B - z \cdot \mathbf I)^{-2} v = \sum_{i = 1}^p \frac{|q_i^{\mathsf{H}} v|^2}{(\lambda_i - z)^2}.
$$

</div>
</details>

</div>
</details>






















<!-- 

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">

</div>
</details>




<div class="callout remark"><span class="label">Remark: Remark for Above Example</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Here is a remark.
</div>

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