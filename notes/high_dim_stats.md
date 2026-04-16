_Overview: Introduction of the section content._

## An Introduction to High-Dimensional Statistics

High-dimensional statistics can refer to a variety of statistical problems and frameworks. In general, it can refer to any problem where the dimension size of the problem, $p$, is large, especially with respect to the sample size $n$.

For example, one of the more popular settings, termed "proportional high dimensional asymptotics", pivots from the more "classical" asymptotics in the following way: while classical asymptotics holds the dimension of the problem $p$ fixed and takes $n \rightarrow \infty$, in PHA, as the sample size $n$ grows, $p$ also grows in some prescribed way <i>with</i> $n$: specifically, the proportion converges to a nonnegative real number, i.e. $\frac{n}{p} \rightarrow \delta \in (0, \infty)$. 

Consider the following simple motivating example: suppose we have a distribution that is multimodal in truth (but this is unknown to us). At the beginning, perhaps we choose a simple Gaussian model for the distribution. As we obtain more data, we would likely modify our model, as the data are highly non-gaussian. We can even make our model more and more complex as we better understand the distribution, which means there are more parameters to learn. In this case, the dimension of the parameters $p$ would grow with $n$.

Let's look at a few examples where "classical" tools may not be effective for analysis, and we'll need to look for something more tailored to high dimensional.

<div class="callout example"><span class="label">Example: Principle Component Analysis</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Suppose we have an observation
$$
X = \lambda u v_*^\top + Z,
$$
where $X \in \mathbb R^{n \times p}$, $\sqrt{\lambda}$ is the signal-to-noise ratio, $u_i \overset{iid}\sim \mathcal N \left(0, \frac{1}{n}\right)$, $Z_{ij} \overset{iid}\sim \mathcal N \left(0, \frac{1}{n}\right)$. The MLE, as we'll discover when we look at this problem in depth, is the first principal component, 
$$
\hat v_{MLE} = \arg\max_v v^\top X^\top X v.
$$
In the high-dimensional setting, as $n, p \rightarrow \infty$, suppose that $\frac{n}{p} \rightarrow \frac{1}{\alpha}$. Then, it turns out,
$$
\lim_{n \rightarrow \infty} | \langle \hat v, v \rangle | = \begin{cases} 0 & \text{if } \lambda \leq \sqrt{\alpha} \\ \sqrt{\frac{1 - \frac{\alpha}{\lambda^2}}{1 + \frac{\alpha}{\lambda}}} & \text{otherwise}\lambda\end{cases} \ .
$$
As we'll also see later, this problem has what's called a <i>phase transition</i>: at some values of an SNR the problem is possible or even easy, and very sharply, becomes impossible.
</div>

<div class="callout example"><span class="label">Example: Risk Estimation</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Take, for example, risk estimation, and specifically, parameter tuning. For example, if we have $k$-fold cross validation, biases persist even in high dimensions.
</div>

<div class="callout example"><span class="label">Example: Sparse Regression</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Another popular example that is of continued interest is the sparse linear regression model:
$$
y = X \beta^* + w,
$$
where $\beta^*$ is $k$-sparse, $X \in \mathbb R^{n \times p}$, and $w$ is some error vector. Commonly examined is the scenario where $n$, $p$, and $k$ all scale with one another in some prescribed way. Even for this seemingly simple problem, in the high-dimensional regime, the optimal $\hat\beta$ is unknown. 

Let's look at something specific. Suppose that $\frac{n}{p} = \delta$ is fixed. Let $\epsilon = \frac{\Vert \beta^* \Vert_0}{p}$, where $\Vert \beta \Vert_0 \triangleq |\{\beta_j : \beta_j > 0\}|$ is the zero norm. Suppose that $p > n$, i.e. that the solution is not unique (although given sparsity, this is also not definite). 

We could impose any number of (unreasonable) restrictions on the problem at the start just to play with it. Suppose we use a popular regularized estimator, the LASSO. Suppose that $w \equiv 0$. Note that for the LASSO, as $\lambda \rightarrow 0$, we converge to the min norm solution $\min_{\Vert \beta \Vert_1}$ for $y = X \beta$. With these crazy assumptions, we can ask: for what values $\delta$ and $\epsilon$ (i.e. the ratios between $n$, $p$, and $k$) can we recover the $k$-sparse solution exactly with the LASSO?
</div>

Now that we have a bit of motivation, let's begin our journey to developing some tools, intuition, and understanding of core topics in high-dimensional statistics.


<!-- 


<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">

</div>
</details>


<div class="callout remark"><span class="label">Remark: Remark for Above Example</span><br/>
<hr style="height:0.01px; visibility:hidden;" />

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