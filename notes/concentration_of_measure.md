_Overview: This section covers concentration of measure, which will be an underlying theme throughout the course._

## An Introduction to Concentration of Measure

Asymptotics is a useful place to turn to if finite-sample results are too difficult to obtain. There are many basic, essential tools for asymptotic analysis, such as the weak law of large numbers.

<div class="callout theorem"><span class="label">Theorem: Weak Law of Large Numbers</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\{X_n\}_{n = 1}^\infty$ be a sequence of integrable and independent and identically distributed random variables with mean $\mu$. Then the <i><strong>weak law of large numbers</i></strong> states that the sample mean, defined for each $n$ as $\bar X_n \triangleq \frac{1}{n} \sum_{i = 1}^n X_i$, converges in probability to $\mu$: for all $\epsilon > 0$,
$$
\mathbb P\left(|\bar X_n - \mu| > \epsilon \right) \rightarrow 0.
$$
</div>

One issue the WLLN fails to address is the rate at which the sample mean converges: it only states that it <i>will</i> converge. We might want to use asymptotics as an approximation of finite sample results, but how do we know how "good" that approximation is? Another useful tool, the central limit theorem, helps to address this.

<div class="callout theorem"><span class="label">Theorem: Central Limit Theorem</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\\{X_n\\}_{n = 1}^\infty$ be a sequence of square integrable and independent and identically distributed random variables with mean $\mu$ and variance $\sigma^2$. Then the <i><strong>central limit theorem</i></strong> states that for the sample mean, defined for each $n$ as $\bar X_n \triangleq \frac{1}{n} \sum_{i = 1}^n X_i$,
$$
\sqrt{n}(\bar X_n - \mu) \overset{d}\rightarrow \mathcal N(0, \sigma^2).
$$
</div>

Having the scaling and the limiting distribution seems to be very helpful for characterizing how good an asymptotic analysis can be. However, it's important to note that the central limit theorem itself is <i>also</i> an asymptotic analysis. How do we characterize how well <i>it</i> can do? We could keep going to higher-order approximations, using an Edgeworth expansion, for example. But instead, we'll pivot to <i>concentration of measure</i>, something that gives us finite sample results (in the form of strict upper bounds that hold at any $n$) at the cost of an exact characterization of error in asymptotics (in the form of telling us the shape of a distribution, for example).

What should these results look like? Let's provide a rough heuristic first to examine what we can look for. Let $\\{X_n\\}_{n = 1}^\infty$ be iid with mean $\mathbb E[X]$ and variance $0$.

The form we want is
$$
\mathbb P\left(\left| \frac{1}{n} \sum_{i = 1}^n X_i - \mathbb E[X]\right| > t\right) \leq [\text{bound}].
$$
What should this bound look like? What rate should it have? Let's base our heuristic on the central limit theorem and make the unreasonable assumption that, in some strange best-case scenario, the central limit has "worked" at some finite sample $n < \infty$. We would have
$$
\mathbb P\left(\sqrt{n}\left| \frac{1}{n} \sum_{i = 1}^n X_i - \mathbb E[X]\right| > \sqrt{n} t\right) \
\overset{(1)}= \ \mathbb P_{Z \sim \mathcal N(0, 1)} \bigg(\big|Z\big| > \sqrt{n} t\bigg) \
= \ 2 \cdot \int_{\sqrt{n}t}^\infty \frac{1}{\sqrt{2\pi}} e^{-\frac{z^2}{2}} dz \
\overset{(2)}\leq \ \frac{2}{\sqrt{2\pi}} \int_{\sqrt{n}t}^\infty \frac{z}{\sqrt{n}t} e^{-\frac{z^2}{2}} dz \ 
= \ \sqrt{\frac{2}{\pi}} \frac{1}{\sqrt{n}t} \exp \left\\{ -\frac{nt^2}{2}\right\\} \
\leq \ c_1 e^{-c_2 \cdot nt^2},
$$
where the equality $(1)$ is our unreasonable application of the CLT at finite $n$ and inequality $(2)$ is because $\frac{x}{\sqrt{n} t} > 1$ on our specifed domain of integration. Hence, in some very nice case scenarios, perhaps we can expect this exponential rate.

Let us now formally examine what sorts of random variables we can say that this holds for.


## Concentration Results

Here, we shall introduce several core results in concentration of measure, including those that are more generally applicable, and those that give us very good bounds.

The first result we state is Hoeffding's inequality.

<div class="callout theorem"><span class="label">Theorem: Hoeffding's Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\\{X_i\\}_{i = 1}^n$ be independent random variables where each is bounded: $a_i \leq X \leq b_i$. Define $X \triangleq \frac{1}{n} \sum_{i = 1}^n X_i$ and define $a \triangleq \frac{1}{n}\sum_{i = 1}^n a_i$ and $b \triangleq \frac{1}{n} \sum_{i = 1}^n b_i$. Then
$$
\mathbb P \left( \left|\frac{1}{n} \sum_{i = 1}^n - \mathbb E[X]\right| > t \right) \leq 2 \cdot \exp \left\\{ -\frac{2}{b - a} n t^2  \right\\}
$$
</div>

Let's compare this result to the simplest case: a finite-sample analysis of Gaussians.

<div class="callout example"><span class="label">Example: Gaussian Bound</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_1, X_2 \ldots X_n \sim \mathcal N(0, \sigma^2)$. Then
$$
\mathbb P \left(\left|\frac{1}{n} \sum_{i = 1}^n X_i\right| > t \right) = \mathbb P_{Z \sim \mathcal N(0, 1)} \left(|\sigma Z| > \sqrt{n} t\right) \leq 2 \int_{\sqrt{n}t}^\infty \frac{z}{\sqrt{n}t} \frac{1}{\sqrt{2\pi\sigma^2}} e^{-\frac{z^2}{2\sigma^2}} dz \ \ \propto \exp \left\\{-\frac{1}{2\sigma^2} nt^2\right\\}
$$
</div>

We notice that the variance, $\sigma^2$, appears in the bound, whereas in Hoeffding's inequality, only the <i>upper bound</i> of the variance, $\frac{1}{4}(b - a)^2$, shows up.

Hence, Hoeffding seems not to be the tightest bound we can get. Bernstein's inequality gives us this.

<div class="callout theorem"><span class="label">Theorem: Bernstein's Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_1, X_2 \ldots X_n$ be zero-mean iid random variables where $|X| < M$. Then
$$
\mathbb P \left(\left|\frac{1}{n} \sum_{i = 1}^n X_i\right| > t\right) \leq 2 \cdot \exp \left\\{-\frac{n^2t^2}{2 \left(n \text{Var}(X) + \frac{Mt}{3}\right)}\right\\}
$$
</div>

Note that the above two results are for bounded random variables. We can conduct a thought experiment: suppose I have $n$ iid random variables from a bounded distribution, for example if $X \in [-1, 1]$. Imagine removing one random variable: this can give us some idea of how the sample can vary. In the case where $X$ is bound within a very small range, removing 1 datapoint should be almost imperceptible. However, imagine removing a single datapoint from an unbounded distribution: although (if the distribution is sufficiently nice) the bulk should be within some "typical" range, there can be outliers that severely impact how the sample looks: that single extreme value could dominate the sum.

This motivates us to think about what should matter when giving concentration results: something that should matter a lot, given our thought experiment above, is the tail behavior of the distribution. How would $\max \\{|X_1|, |X_2| \ldots |X_n|\\}$ behave?

<div class="callout example"><span class="label">Example: Gaussian Tail Behavior</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_1, X_2 \ldots X_n \overset{iid}\sim \mathcal N(0, 1)$. Then for $X_{(n)} \triangleq \max \\{X_1, X_2 \ldots X_n\\}$, i.e. the $n$th order statistic,
$$
\frac{X_{(n)}}{\sqrt{2 \log n}} \overset{p}\rightarrow 1.
$$
</div>

We see from this example that the growth is order $\sqrt{\log n}$, which is <i>much</i> slower than linear. So the growth is not very crazy at all. What about for heavier tails, i.e. tails not as nice as Gaussian?

So, in our journey to extend concentration of measure results from bounded to unbounded random variables, we should turn our attention to studying tail behaviors. To this end, let us introduce a few key definitions first.

<div class="callout definition"><span class="label">Definition: Subgaussian Random Variable</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
A random variable $X$ is <i><strong>subgaussian</i></strong> if there exist constants $c, C$ such that 
$$
\mathbb P(|X| > t) \leq C \cdot e^{-ct^2}.
$$
</div>

<div class="callout definition"><span class="label">Definition: Subexponential Random Variable</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
A random variable $X$ is <i><strong>subexponential</i></strong> if there exist constants $c, C$ such that 
$$
\mathbb P(|X| > t) \leq C \cdot e^{-ct}.
$$
</div>

Using these definitions, we may establish concentrataion of measure results for these two types of (possibly) unbounded random variables.

<div class="callout proposition"><span class="label">Proposition: Bound for Subgaussian RVs</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_1, X_2 \ldots X_n$ be iid copies of a subgaussian random variable $X$. Then there exist positive constants $c_1, c_2$ such that
$$
\mathbb P\left(\left|\frac{1}{n} \sum_{i = 1}^n X_i - \mathbb E[X_i]\right| > t \right) \leq c_1 e^{-c_2 n t^2}.
$$
</div>

<div class="callout proposition"><span class="label">Proposition: Bound for Subexponential RVs</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_1, X_2 \ldots X_n$ be iid copies of a subexponential random variable $X$. Then there exist positive constants $c_1, c_2, c_3, c_4$, and $c$ such that
$$
\mathbb P\left(\left|\frac{1}{n} \sum_{i = 1}^n X_i - \mathbb E[X_i]\right| > t \right) \leq \begin{cases} c_1 e^{-c_2 n t^2} & t \leq c \\ c_3 e^{-c_4 nt} & t > c \end{cases} \ \ .
$$
</div>

Intuitively, we should interpret this result as: with heavier tails than a Gaussian, the possibility of the sum being governed by a single extreme datapoint away from the bulk is noticeably greater.

With these established results on simple random variables, we can move on to discussing concentration of more interesting objects.

## Concentration of Functions of Random Variables

Suppose we have a function $F : \mathbb R^n \rightarrow \mathbb R$ that takes a random vector as input. As before, we're interested in how far something can deviate away from its expectation: in this case, we want to study $|F(X_1, X_2 \ldots X_n) - \mathbb E[F(X_1, X_2 \ldots X_n)]|$. Should it hold that $F(X_1, X_2 \ldots X_n) \approx \mathbb E[F(X_1, X_2 \ldots X_n)]$? For what sorts of functions should this be true? 

For example, it obviously doesn't hold for $F(X_1, X_2 \ldots X_n) = X_1$, an extreme case that throws away everything but one datapoint... but we're interested in functions that have more notable applications. Another "extreme" most simple case is an example we've already thoroughly studied: $F(X_1, X_2 \ldots X_n) = \frac{1}{n} \sum_{i = 1}^n X_i$. What about eigenvalues or singular values of a random matrix?

A preliminary strategy worth examining is to apply tools we are familiar with to find an upper bound, i.e. find that
$$
\mathbb P (|F(X_1, X_2 \ldots X_n) - \mathbb E[F(X_1, X_2 \ldots X_n)]| > t) = o_p(1).
$$
The simplest tool we have is Markov's inequality.

<div class="callout theorem"><span class="label">Theorem: Markov's Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X$ be an integrable nonnegative random variable. Then
$$
\mathbb P(X > t) \leq \frac{\mathbb E[X]}{t}.
$$
This result is known as <i><strong>Markov's inequality</i></strong>.
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
We have that
$$
\mathbb P(X > t) \ \overset{(1)}= \ \mathbb E[\mathbb I_{\\{X > t\\}}] \ \overset{(2)}\leq \ \mathbb E\left[\mathbb I_{\\{X > t\\}}  \cdot \frac{X}{t}\right] \ \overset{(3)}\leq \ \frac{\mathbb E[X]}{t}.
$$
Equality $(1)$ is by definition of expectation over the probability measure, inequalities $(2)$ and $(3)$ are because $\frac{X}{t} > 1$ on the domain of integration and $\mathbb I_{\\{X > t\\}} \leq 1$, respectively, and then applying the nonnegativity of $X$ and the monotonicity of expectation.
</div>
</details>

<div class="callout corollary"><span class="label">Corollary: Chebyshev's Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X$ be a squarae integrable random variable. Then
$$
\mathbb P(|X - \mathbb E[X]| > t) \leq \frac{\text{Var}(X)}{t^2}.
$$
This result is known as <i><strong>Chebyshev's inequality</i></strong>.
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
The proof is a simple application of Markov's inequality.
$$
\mathbb P(|X - \mathbb E[X]| > t) = \mathbb P((X - \mathbb E[X])^2 > t^2) \leq \frac{\mathbb E[(X - \mathbb E[X])^2]}{t^2} = \frac{\text{Var}(X)}{t^2}.
$$
</div>
</details>

If we simply apply Chebyshev's inequality to our function, we obtain
$$
\mathbb P((F(X_1, X_2 \ldots X_n) - \mathbb E[F(X_1, X_2 \ldots X_n)])^2 > t^2) \leq \frac{\text{Var}(F(X_1, X_2 \ldots X_n))}{t^2}.
$$

We make an observation here: if we can control the (possibly the average) gradient of this function and ensure that it's not too steep, we shouldn't move too crazily with small perturbations of the input, and the variance could be bounded by something nice.

To this end, let's look at functions that have precisely this nice property.

<div class="callout definition"><span class="label">Definition: Lipschitz Continuity</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
A function $F : \mathbb R^n \rightarrow \mathbb R$ is <i><strong>Lipschitz continuous</i></strong> with constant $L$ if
$$
|F(x_1, x_2 \ldots x_n) - F(y_1, y_2 \ldots y_n)| \leq L \cdot \Vert x - y \Vert_2.
$$
</div>

And this definition directly plays into our gradient requirement via the following proposition.

<div class="callout proposition"><span class="label">Proposition: Mean Value Inequality on Lipschitz Functions</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let a function $F : \mathbb R^n \rightarrow \mathbb R$ be $L$-Lipschitz and differentiable. Then
$$
\Vert \nabla F \Vert_2 \leq L.
$$
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
We have by differentiability and Taylor's theorem that
$$
f(x + h) = f(x) + \nabla f(x)^\top h + o(\Vert h \Vert).
$$
Consider $\nabla f(x) \neq 0$ (the proof of the opposite case is obvious). Let $\epsilon > 0$ be arbitrary, and define
$$
h = \epsilon \cdot \frac{\nabla f(x)}{\Vert \nabla f(x) \Vert_2}.
$$
Then we have
$$
f(x + h) - f(x) = \epsilon \cdot \frac{\nabla f(x)^\top \nabla f(x)}{\Vert \nabla f(x) \Vert_2} + o(\Vert h \Vert_2) = \epsilon \cdot \Vert \nabla f(x) \Vert_2 + o(\epsilon).
$$
Hence we have
$$
\left|\Vert \nabla f(x) \Vert_2 + \frac{o(\epsilon)}{\epsilon}\right| = \left|\frac{f(x + h) - f(x)}{\epsilon}\right| \leq \frac{L \Vert x + h - x \Vert_2}{\epsilon} = L.
$$
As $\epsilon \rightarrow 0$, $\frac{o(\epsilon)}{\epsilon} \rightarrow 0$ by definition.
</div>
</details>

With these tools to describe the functions we analyze, we shall focus on an interesting phenomenon: self-averaging.

## Self-Averaging Functions

The notion of self-averaging comes from physics. The mathematical idea, for now, can be interpreted simply as: can functions that are not the sample mean still act like some sort of averaging? If so, for which types of functions does self-averaging occur?

We provide an example first. In fact, this example will be very useful for our future analysis in random matrix theory.

<div class="callout example"><span class="label">Example: Self-Averaging of the Quadratic Form</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Suppose we have a random vector (overloading notation a little bit) $\vec X_n \triangleq (X_i)_{i = 1}^n$, where elements $X_i \overset{iid}\sim N\left(0, \frac{1}{n} \right)$. The scaling $\frac{1}{n}$ will be explained in more detail later. Suppose further we have a matrix $A_n$ that is symmetric (it is deterministic for now, but we'll later look at random versions where $A_n$ is independent of $\vec X_n$) and that $|\lambda_{\max}(A)| \leq L$. Define the function
$$
F_n(X_n) \triangleq \vec X_n^\top A_n \vec X_n.
$$
We are interested in establishing an upper bound on the probability that $F_n$ deviates far from its mean:
$$
\mathbb P(|F_n(\vec X_n) - \mathbb E[F_n(\vec X_n)]| > t) \leq [\text{bound}].
$$
It turns out that we get this bound to be
$$
\mathbb P(|F_n(\vec X_n) - \mathbb E[F_n(\vec X_n)]| > t) \leq c \cdot \frac{L^2}{n t^2}.
$$
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
First we calculate the mean:
$$
\mathbb E[F_n(X_n)] = \mathbb E[X_n^\top A_n X_n] = \mathbb E[\text{Tr}(X_n^\top A_n X_n)] = \text{Tr}(A_n) \cdot \mathbb E[X_n^\top X_n] = \frac{1}{n} \mathbf I_n \cdot \text{Tr}(A_n).
$$
We can now look at the deviation from the mean:
$$
\mathbb P(|F_n(X_n) - \mathbb E[F_n(X_n)]| > t) = \mathbb P(|F_n(X_n) - \frac{1}{n} \text{Tr}(A_n)|) \leq \frac{1}{t^2} \mathbb E[(X_n^\top A_n X_n - \frac{1}{n} \text{Tr}(A_n))^2],
$$
where the last inequality is from applying Chebyshev. Now, we can do some algebraic manipulations on this quantity:
$$
\frac{1}{t^2} \mathbb E[(X_n^\top A_n X_n - \frac{1}{n} \text{Tr}(A_n))^2] = \frac{1}{t^2} \cdot \mathbb E\left[\left(\sum_{i = 1}^n x_i^2 A_{ii} + 2 \sum_{i < j} A_{ij}x_ix_j - \frac{1}{n}\sum_{i = 1}^n A_{ii}\right)^2\right] = \frac{1}{t^2} \mathbb E\left[\left(\sum_{i = 1}^n \left(x_i - \frac{1}{n}\right)A_{ii} + \sum_{i = 1}^n \sum_{j \neq i} A_{ij} x_i x_j\right)^2\right].
$$
Now, utilizing the fact that $2 \mathbb E[AB] \leq \mathbb E[A^2] + \mathbb E[B^2]$ (applying this identity on the cross term in the quadratic), we have that
$$
\mathbb P\left(\left|F_n(X_n) - \frac{1}{n} \text{Tr}(A_n)\right| > t\right) \leq \frac{2}{t^2} \mathbb E\left[\left(\sum_{i = 1}^n \sum_{j \neq i} A_{ij} x_i x_j\right)^2\right] + \frac{2}{t^2} \mathbb E\left[\left(\sum_{i = 1}^n \left(x_i^2 - \frac{1}{n}\right)A_{ii}\right)^2\right].
$$
Let us now evaluate this expression. We have, when expanding the second term,
$$
\mathbb E\left[\left(\sum_{i = 1}^n \left(x_i^2 - \frac{1}{n}\right)A_{ii}\right)^2\right] = \sum_{i = 1}^n \sum_{j = 1}^n A_{ii} A_{jj} \cdot \mathbb E\left[\left(x_{i}^2 - \frac{1}{n}\right)\left(x_j^2 - \frac{1}{n}\right)\right] \overset{(1)}= \sum_{i = 1}^n A_{ii}^2 \mathbb E\left[(x_i^2 - \frac{1}{n})^2\right] \overset{(2)}= \frac{2}{n^2} \sum_{i = 1}^n A_{ii}^2,
$$
where the first equality holds because independence between the elements of $x_i$ causes all cross terms to cancel. The second equality is because the expectation is precisely the variance of $x_i^2$, and we know the 4th moment of a zero-mean Gaussian.

Now, moving on to the first term in the expression, we again have by independence that "cross" terms disappear, and we're left with
$$
\mathbb E\left[\left(\sum_{i = 1}^n \sum_{j \neq i} A_{ij} x_i x_j\right)^2\right] = 2 \sum_{i = 1}^n \sum_{j \neq i} A_{ij}^2 \cdot \mathbb E\left[x_i^2 x_j^2\right] = \frac{2}{n^2} \sum_{i \neq j} A_{ij}^2.
$$
Combining the two terms we have that
$$
\mathbb P\left(\left|F_n(X_n) - \frac{1}{n} \text{Tr}(A_n)\right| > t\right) \leq \frac{2}{t^2} \left(\frac{2}{n^2} \sum_{i \neq j} A_{ij}^2 +  \frac{2}{n^2} \sum_{i = 1}^n A_{ii}^2\right) = \frac{4}{n^2t^2} \Vert A \Vert_F^2,
$$
where $\Vert \cdot \Vert_F$ denotes the Frobenius norm. Using that $\Vert A \Vert_F^2 = \sum_{i = 1}^n \lambda_i^2 \leq \sum_{i = 1}^n L^2 = n L^2$, we finally have that
$$
\mathbb P\left(\left|F_n(X_n) - \frac{1}{n} \text{Tr}(A_n)\right| > t\right) \leq \frac{4}{n^2 t^2} n L^2 = \frac{4L^2}{nt^2}.
$$
Hence, the absolute constant $c$ is computed to be 4.
</div>
</details>

With this example under our belt, we can broaden our scope of interest. Our more general goal is to show that this self-averaging phenomenon happens for generic Lipschitz functions (i.e. nonlinear), as long as the random variables in question have sufficiently light tails. For example, what if we move from the deterministic matrix in the previous example and look at the case when $A$ is a matrix of iid Gaussian random variables? Suppose we're interested in the largest (in magnitude) eigenvalue, i.e. for a sequence of random matrices $A_n$, we define $F_n \triangleq \lambda_{\max}(A_n)$. It turns out that this is nearly constant&mdash;and we may conclude that although this function is not a simple averaging, there is some form of self-averaging occuring.

For self-averaging functions, we can roughly consider two cases, or two types of results: we can get bounds of the form
$$
\mathbb P(|F_n(X_n) - \mathbb E[F_n(X_n)]| > t) \ \leq \ c_1 \cdot e^{-c_2 nt^2},
$$
or we can get bounds of the (less ideal) form
$$
\mathbb P(|F_n(X_n) - \mathbb E[F_n(X_n)]| > t) \ \leq \ \frac{\text{Var}(F_n(X_n))}{t^2} \ \simeq \ \frac{1}{\text{polylog}(n)}.
$$
The latter type, as we may surmise from the form, can be called Chebyshev-type inequalities, contrasted with the former bounds, which are exponential-type inequalities. This type is more general and requires fewer and less strong assumptions, at the cost of weaker results. Let us explore what tools we can use in conjunction with Chebyshev-type inequalities.

### Chebyshev-Type Inequalities

We see that a major task when looking at the Chebyshev-type inequalities will be to bound the variance, $\text{Var}(F_n)$. Let's look at two extreme toy examples to develop an intuition on how to proceed.

<div class="callout example"><span class="label">Example: Toy Examples for $\text{Var}(F)$</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Suppose for simplicity that $\text{supp}(F) \subseteq [0, 1]$.

First, in one extreme, take $F_n(X_n) \equiv c$ for a constant $c \in \mathbb R$. Then $\text{Var}(F)$ is 0, and the concentration is trivial.

On the opposite extreme, suppose that the function oscillates wildly on its domain $[0, 1]$ with tall spikes. Then it seems possible that there could be no concentration. Here, with even tiny movements in the inputs of $F_n$, there is a wild change in $F_n$. 

However, we don't expect the typical functions of interest to behave in this way: for our favorite example, the eigenvalues of a random matrix $A$ shouldn't spike in such an extreme way with slight modifications of the elements of $A$. Hence, for a function that doesn't change rapidly with small movements of the input, we should be able to upper bound. 
</div>

Two tools we shall use for bounding the variance $\text{Var}(F_n)$ will be the <i>Efron-Stein inequality</i> and the <i>Poincaré inequality</i>. Let's start with Efron-Stein.

<div class="callout theorem"><span class="label">Theorem: Efron-Stein Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $Z_n = F_n(X_n)$, where $X_n \triangleq (X_i)_{i = 1}^n$ are iid. Let $\tilde X_n \triangleq (\tilde X_i)_{i = 1}^n$ be independent copies of $X_i$, and denote $\tilde X_n^{(i)}$ denote the vector where $X_i$ has been exchanged with $\tilde X_i$, i.e. $\tilde X_n^{(i)} \triangleq (X_1 \ldots X_{i - 1}, \tilde X_i, X_{i + 1} \ldots X_n)$.
$$
\text{Var}(Z_n) \leq \frac{1}{2} \sum_{i = 1}^n \mathbb E \left[\left(F_n(X_n) - F_n(\tilde X_n^{(i)})\right)^2\right].
$$
</div>

Let's consider a toy example in applying the E-S inequality, U-statistics.

<div class="callout example"><span class="label">Example: E-S on U-Statistics</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_n = (X_i)_{i = 1}^n$ have iid components. Let $g : \mathbb R^2 \rightarrow \mathbb R$ be a bounded function, i.e. $\Vert G \Vert_\infty \leq C$ for some $C \in \mathbb R^+$.

Define the statistic
$$
F_n(X_n) \triangleq {n \choose 2}^{-1} \sum_{j < k} g(x_j, x_k).
$$
We have by Chebyshev that
$$
\mathbb P(|F_n - \mathbb E[F_n]| > t) \leq \frac{\text{Var}(F_n)}{t^2}.
$$
We now use E-S to bound the variance. For brevity, define the random variables
$$
\Delta(X_i) \triangleq F_n(X_n) - F_n(\tilde X_n^{(i)}).
$$
We now have
$$
\Delta(x_i) = {n \choose 2}^{-1} \left(\sum_{j \neq i} g(x_i, x_j) - \sum_{j \neq i}g(\tilde x_i, x_j)\right),
$$
because the terms that don't include $X_i$ are the same in the subtraction and are hence canceled out. This implies, by the triangle ineqaulity and because $\frac{1}{2}(n - 1)^2$ is an upper bound for $n \choose 2$, 
$$
|\Delta (x_i)| \leq \frac{2}{(n - 1)^2} \sum_{j \neq i} |g(x_i, x_j) - g(\tilde x_i, x_j)| \leq \frac{2}{(n - 1)^2} \cdot 4C \cdot (n - 1) = \frac{4C}{n - 1}.
$$
The last inequality is obtained by the boundedness of $g$.
</div>

Let us move on to the other tool we can use for bounding the variance: the Poincaré inequality. This work became foundational for subsequent concentration results after its publication. We'll provide a bit of intuition first: we hope that, for well-behaved differentiable functions,
$$
\text{Var}(F) \leq \mathbb E \left[\Vert \nabla F \Vert^2 \right].
$$
What can we say about random variables that satisfy this for every $F$?

<div class="callout definition"><span class="label">Definition: Poincaré Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
A random variable $X$ satisfies the <i><strong>Poincaré inequality</i></strong> with constant $c$ if, for every differentiable function $f$, we have that
$$
\text{Var}(f(X)) \leq c \cdot \mathbb E\left[f'(X)^2\right].
$$
</div>

We present a few examples of random variables that satisfy the Poincaré inequality.

<div class="callout example"><span class="label">Example: Random Variables Satisfying the Poincaré Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Here are random variables that satisfy the Poincaré inequality.
<ol type="1">
  <li>Uniform random variables $X \sim \text{Unif}([a, b])$ satisfy the Poincaré inequality with constant $\frac{1}{2}(b - a)$.</li>
  <li>Exponential random variables $X \sim \text{Exp}(\lambda)$ satisfy the Poincaré inequality. For $\lambda = 1$, $c = 4$.</li>
  <li>Gaussian random variables $X \sim \mathcal N(0, 1)$ satisfy the Poincaré inequality with constant $c = 1$.</li>
</ol>
That exponential random variables satisfy the Poincaré inequality illustrates that it applies not only to subgaussian random variables.
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
We present the proof of the first case, where $X \sim \text{Unif}([a, b])$: we show that, for all differentiable functions $f$, 
$$
\text{Var}(f(X)) \leq \frac{1}{2} (b - a) \cdot \mathbb E \left[(f'(X))^2\right].
$$
We shall use a standard ghost sample trick: let $Y$ be an independent copy of $X$, and we have that $\text{Var}(f(X)) = \frac{1}{2} \mathbb E_{X, Y} \left[(f(X) - f(Y))^2\right]$. Hence, our goal is to bound $\mathbb E_{X, Y} \left[(f(X) - f(Y))^2\right]$. Because $f$ is differentiable, we have that
$$
f(X) - f(Y) = \int_Y^X f'(t) dt.
$$
We then have that
$$
\bigg(f(X) - f(Y)\bigg)^2 \ 
= \ \bigg(\int_Y^X f'(t) dt\bigg)^2 \ 
\overset{(1)}\leq \ \int_Y^X \left(f'(t)\right)^2 dt \cdot \int_Y^X dt \ 
\overset{(2)}\leq \ \int_a^v \left(f'(t)\right)^2 dt \cdot \int_a^b dt \ 
= \ (b - a) \mathbb E\left[f'(W)^2\right],
$$
where $(1)$ is by Cauchy-Schwarz, $(2)$ is because the functions are nonnegative and $X, Y$ are bounded by $[a, b]$, and $W \sim \text{Unif}([a, b])$.
Hence, we conclude that
$$
\text{Var}(f(X)) \leq \frac{1}{2} \mathbb E\bigg[\bigg(f(X) - f(Y)\bigg)^2\bigg] \leq \frac{1}{2} (b - a) \cdot \mathbb E \left[(f'(X))^2\right].
$$
</div>
</details>

We might hope that the Poincaré inequality applies broadly to all subgaussian random variables... But this turns out not to be true. We provide an example below.

<div class="callout example"><span class="label">Example: Random Variables Satisfying the Poincaré Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Note that all bounded random variables are subgaussian. 

......

Hence, we can make an educated guess that the Poincaré inequality will not be helpful for bimodal or multimodal distributions.
</div>

One of the very useful properties of the Poincaré inequality is <i>tensorization</i>, meaning it extends very naturally into higher dimensions without the constant $c$ exploding.

<div class="callout theorem"><span class="label">Theorem: Tensorization of the Poincaré Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_1, X_2 \ldots X_n$ be independent random variables satisfying the Poincaré inequality with constants $c_i$. Then the random vector $X = (X_i)_{i = 1}^n$ satisfies the Poincaré inequality with constant $c = \max \\{c_1, c_2 \ldots c_n\\}$, i.e.
$$
\text{Var}(f(X)) \leq c \cdot \mathbb E\left[\Vert \nabla f(X) \Vert_2^2\right]
$$
This is called the <i><strong>tensorization property</i></strong> of the Poincaré inequality.
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
The proof of the tensorization property is via an application of Efron-Stein, which recall states that, for an independent copy of $\vec X_n$, denoted $\tilde X_n$, we have that 
$$
\text{Var}(f(X)) \leq \frac{1}{n} \sum \mathbb E \left[\left(f(\vec X_n) - f(\tilde X_n) \right)^2\right].
$$
Denote, as before, the vector $\vec X_n$ with sample $i$ replaced by its copy $\tilde X_i$ as $\tilde X_n^{(i)}$, and denote the one-dimensional sections of $f$ $g_i(X_i) \triangleq f(X_1, X_2 \ldots X_i \ldots X_n)$. We then have that, for each $i \in [n]$, 
$$
\mathbb E_{\vec X_n, \tilde X_i} \left[\left(f(\vec X_n) - f(\tilde X_n^{(i)})\right)^2\right] \
\overset{(1)}= \ \mathbb E_{\vec X_{-i}} \left[\mathbb E_{X_i, \tilde X_i} \left(g(X_i) - g(\tilde X_i)\right)^2\right] \
\overset{(2)}= \ 2 \cdot \mathbb E_{\vec X_{-i}} \left[\text{Var}(g_i(X_i))\right] \
\overset{(3)}\leq \ 2 c_i \cdot \mathbb E_{\vec X_{-1}} \left[\mathbb E_{X_i} \left[({g_i}^{'}(X_i))^2\right]\right] \
\overset{(4)} \leq \ 2c_i \cdot \mathbb E\left[\left(\frac{\partial}{\partial X_i} f(\vec X_n)\right)^2\right].
$$
$(1)$ is by the tower property and using our new notation to examine one-dimensional sections in order to apply Efron-Stein, $(2)$ is by using the ghost sample trick, $(3)$ is the application of Efron-Stein, and $(4)$ is by the tower property and because $f$ is differentiable.

Therefore, returning to our whole expression, we have
$$
\text{Var}(f(X)) \ 
\leq \ \frac{1}{2} \sum_{i = 1}^n 2 \cdot c_i \mathbb E_{\vec X} \left[\left(\frac{\partial}{\partial X_i} f(\vec X_n)\right)^2\right] \
\leq \max \\{c_1, c_2 \ldots c_n\\} \cdot \mathbb E_{\vec X_n} \left[\sum_{i = 1}^n \left(\frac{\partial}{\partial X_i} f(\vec X_n)\right)^2\right] \
= \ c \cdot \mathbb E_{\vec X_n} \bigg[\Vert \nabla f(X) \Vert_2^2\bigg].
$$
</div>
</details>

We return to our favorite example: the max singular value.

<div class="callout example"><span class="label">Example: Poincaré on a Gaussian Random Matrix</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $A$ be a random matrix with $A_{ij} \overset{iid}\sim \mathcal N(0, 1)$. We'd like to show something of the form
$$
\mathbb P(|\sigma_\max(A) - \mathbb E[\sigma_\max(A)] > t).
$$
This amounts to bounding the variance of a very complicated function of $A$, namely, $F(A) = \sup_{\Vert x \Vert_2 = 1} \Vert Ax \Vert_2$. We will use the Poincaré inequality, then use Lipschitzness of $F(A)$ to upper bound the gradient of the norm. In the end, we obtain a bound
$$
\mathbb P(|\sigma_\max - \mathbb E[\sigma_\max]| > t) \leq \frac{1}{t^2}.
$$
</div>

<details class="collapsible">
<summary>Proof.</summary>
<div class="collapsible__content">
We have immediately by Chebyshev's inequality that
$$
\mathbb P(|\sigma_\max - \mathbb E[\sigma_\max]| > t) \leq \frac{\text{Var}(\sigma_\max(A))}{t^2} \leq \frac{\mathbb E[\Vert \nabla \sigma_\max (A) \Vert_2]^2}{t^2},
$$
and the second inequality is our application of Poincaré, which for Gaussian random variables has constant $c = 1$.

Our objective, then, is to bound $\mathbb E\left[\Vert \nabla \sigma_\max(A) \Vert_2^2\right]$. If we're able to prove that $\sigma_\max(A)$ is an $L$-Lipschitz function, the norm of the gradient is upper bounded by $L$. It turns out that this is true.

<details class="collapsible">
<summary>Proof that $\sigma_\max(A)$ is Lipschitz.</summary>
<div class="collapsible__content">
Let $C$ and $D$ be matrices of the same dimension. Then 
$$
\sigma_\max(C + D) = \sup_{\Vert x \Vert_2 = 1} \Vert (C + D) x \Vert_2 \overset{(1)}\leq \sup_{\Vert x \Vert_2 = 1} \left(\Vert C x \Vert_2 + \Vert D x \Vert_2\right) \overset{(2)}\leq \sup_{\Vert x_1 \Vert_2 = 1} \Vert C x_1 \Vert_2 + \sup_{\Vert x_2 \Vert_2 = 1} \Vert D x_2 \Vert_2 = \sigma_\max(C) + \sigma_\max(D),
$$
where inequality $(1)$ is the triangle inequality and $(2)$ is by subadditivity of the supremum.

Now, define the matrices $C$ and $D$ to be $C \triangleq A - B$, and $D \triangleq B$. We then have that
$$
\sigma_\max(A) \leq \sigma_\max(B) + \sigma_\max(A - B) \implies |\sigma_\max(A) - \sigma_\max(B)| \leq \sigma_\max(A - B).
$$
Hence, the function $\sigma_\max(A)$ is 1-Lipschitz.
</div>
</details>
We have shown that the function is 1-Lipschitz, so we know that the norm of the gradient is also upper bounded by 1. Hence, we have that
$$
\mathbb P(|\sigma_\max - \mathbb E[\sigma_\max]| > t) \leq \frac{1}{t^2}.
$$
</div>
</details>

Two remarks about the above example.

<div class="callout remark"><span class="label">Remark: Scaling of the Gaussian Entries</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Note that the bound is $\frac{1}{t^2}$, which doesn't converge in $n$. We can change the scaling, e.g. to $A_{ij} \overset{iid}\sim \mathcal N\left(0, \frac{1}{n}\right)$, to obtain something that converges in $n$.
</div>

<div class="callout remark"><span class="label">Remark: Other Lipschitz Spectrum-Related Functions</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
It turns out that all functions of the type $\sigma(\boldsymbol\cdot)$ or $\lambda(\boldsymbol\cdot)$ are Lipschitz. The generalization is called Weyl's theorem.
</div>

This concludes our introduction to Chebyshev-type inequalities. Now, let us examine bounds that are stronger, but also require stronger assumptions about the random variables.

### Exponential-Type Inequalities

Exponential type inequalities have tighter bounds. Where the Chebyshev-type inequalities can give us $O((\text{polylog}(n))^{-1})$, exponential-type inequalities would give us bounds of the form $O(e^{-c n t^2})$.

We'll introduce two tools that form a sort of analog to the E-S and Poincaré inequalities. The first is McDiarmid's inequality, the exponential-type "version" of Efron-Stein.

<div class="callout theorem"><span class="label">Theorem: McDiarmid's Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_n = (X_i)_{i = 1}^n$ be a random vector of independent random variables. Let $\tilde X_n$ be an independent copy of $X_n$, and define $\tilde X_n^{i}$ to be the vector $X_n$ replacing $X_i$ with $\tilde X_i$.

Let $f$ be a function of $X_n$. Suppose that there exist constants $c_1, c_2 \ldots c_n$ such that
$$
|f(X_n) - f(\tilde X_n^{(i)})| \leq c_i.
$$
Then, for all $t > 0$, we have, defining $c = \sum_{i = 1}^n c_i$,
$$
\mathbb P (|f(X_n) - \mathbb E[f(X_n)]| > t) \leq e^{-\frac{2t^2}{c}}.
$$
This result is known as <i><strong>McDiarmid's inequality</i></strong>.
</div>

Let's compare this with Efron-Stein, which upper bounds the variance with $\sum_{i = 1}^n \mathbb E\left[\Delta(X_i)^2\right]$. Note that this result is a result over the expectation: if the function is nice enough <i>in expectation</i>, then we can get an upper bound. The idea with McDiarmid is that we can get a <i>better</i> upper bound if we are able to make a deterministic (or at least almost sure) statement about the function.

Now, we turn to the log-Sobolev inequality, the analog of Poincaré. We first define what it means for a distribution to satisfy the log-Sobolev inequality.

<div class="callout definition"><span class="label">Definition: Logarithmic-Sobolev Inequality</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
A probability measure $\mu$ on $\mathbb R^n$ satisfies the <i><strong>logarithmic Sobolev inequality</i></strong> with constant $c$ if, for all differentiable functions $f : \mathbb R^n \rightarrow \mathbb R$, we have that
$$
\mathbb E\left[f^2 \log f^2\right] - \mathbb E[f^2] \cdot \log (\mathbb E[f^2]) \leq 2c \cdot \mathbb E[\Vert \nabla f \Vert_2^2].
$$
</div>

To gain a bit of insight into the log-Sobolev inequality, let's do an interesting exercise. Define the function $\phi(z) = z \log z$ for $z > 0$. Observe that this exercise would functionally require $f^2 > 0$. Then we have that the log-Sobolev inequality takes the form
$$
\mathbb E[\phi(f^2)] - \phi(\mathbb E[f^2]) \leq 2c \cdot \mathbb E\left[\Vert \nabla f \Vert_2^2\right]. 
$$
This form looks somewhat reminiscent of Chebyshev. Note that the expression is bounded below by 0 by Jensen's inequality, as $\phi$ is convex.

There are a couple things to note here. First, note that this form captures variations of the function in a specific way: the expression on the left represents the global variation, while the one on the right represents the local variation. The aim to bound some type of variation in the function should be reminiscent of Poincaré. However, the stronger condition comes from the stronger emphasis on the values of $X$ that push $f$ toward $\infty$, and hence induces a stronger requirement on the tail behavior and is more stringent than Poincaré. 

What about the $\log$ term? First, using the $\log$ represents using the entropy rather than variance, and second, perhaps we can also interpret it as ensuring that exponential tails are excluded.

To show that log-Sobolev has a more stringent condition than Poincaré, we show that the former implies the latter.

<div class="callout proposition"><span class="label">Proposition: Log-Sobolev Implies Poincaré</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
For any square integrable function $f$,
$$
\text{Var}(f) \leq \int f^2 \log f^2 d\mu - \int f^2 d\mu \cdot \log \left(\int f^2 d\mu\right).
$$
Hence, log-Sobolev implies Poincaré, but not vice versa.
</div>

Let us now look at a few examples.

<div class="callout example"><span class="label">Example: Random Variables Satisfying Poincaré</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
A random variable $X \sim \text{Unif}([0, 1])$ satisfies the log-Sobolev inequality with constant $c = \frac{1}{2\pi^2}$.

A random variable $X \sim \text{Exp}(\lambda)$ does <i>not</i> satisfy log-Sobolev.

For a measure $\mu$ on $\mathbb R$ with a density $\mu(x) = \frac{1}{z} e^{-\nu(x)}$ where $\nu(x) - \frac{x^2}{2c}$ is convex, $\mu$ satisfies the log-Sobolev inequality with constant $c$. This essentially means that the tails are gaussian or lighter. For instance, $\mathcal N(0, 1)$ satisfies log-Sobolev with $c = 1$.
</div>

Just as with Poincaré, log-Sobolev enjoys a tensorization property.

<div class="callout theorem"><span class="label">Theorem: Tensorization Property of Log-Sobolev</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
For probability measures $\mu_1, \mu_2 \ldots \mu_n$ on $\mathbb R$, suppose that each $\mu_i$ satisfies the log-Sobolev inequality with constant $c_i$. Then $\mu_1 \otimes \mu_2 \ldots \otimes \mu_n$ satisfies log-Sobolev with constant $c = \max \\{c_1, c_2 \ldots c_n\\}$.
</div>

The proof of this is in Wainwright's book.

For example, take our previous proposition about random variables with densities of some exponential form.

<div class="callout example"><span class="label">Example: Tensorization on Exponential Densities</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
For a distribution $\mu$ on $\mathbb R^n$, if the density $\mu(x) = \frac{1}{2} e^{-\nu(x)}$ and $\nu(x) - \frac{\Vert x \Vert^2}{2c}$ is convex, then the distribution satisfies log-Sobolev with constant $c$. 
</div>

Now, we introduce a core theorem that makes the connection between the log-Sobolev inequality and subgaussian concentration.

<div class="callout theorem"><span class="label">Theorem: Herbst's Theorem</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\mu$ be a measure on $\mathbb R$ satisfying the log-Sobolev inequality with constant $c$. Then, for all $L$-Lipschitz functions $f$, we have that
$$
\mathbb P \left(\left|f(x) - \mathbb E[f(x)]\right| > t \right) \leq 2 \cdot \exp \left\\{-\frac{t^2}{2c L^2}\right\\}.
$$
This theorem is often referred to as <i><strong>Herbst's theorem</i></strong>.
</div>

We return to our favorite example, the largest singular value of a random matrix.

<div class="callout example"><span class="label">Example: Log-Sobolev on a Gaussian Random Matrix</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $A \in \mathbb R^{n \times p}$ where $A_{ij} \overset{iid}\sim \mathcal N(0, 1)$. Let $\sigma_\max(A)$ denote the largest (in magnitude) singular value of $A$. Then, because the entries are iid gaussian, each is log-Sobolev with constant 1, and the entire matrix hence satisfies log-Sobolev by tensorization.

We have shown before that $\sigma_\max(A)$ is 1-Lipschitz, i.e. $|\sigma_\max(A) - \sigma_\max(B)| \leq \Vert A - B \Vert_F$. Hence, we have that
$$
\mathbb P \left(\left|\sigma_\max(A) - \mathbb E[\sigma_\max(A)]\right| > t \right) \leq 2 \exp \left\\{-\frac{t^2}{2}\right\\}.
$$
It turns out that the analysis of the term $\mathbb E[\sigma_\max(A)]$ is nontrivial. We will eventually use something called Gordon's comparison inequality, Dudley's inequalities and discretization with union bounds, which give $\sqrt{n} - \sqrt{p} \leq \mathbb E[\sigma_\max(A)] \leq \sqrt{n} + \sqrt{p}$. For now, if we accept this, we have that
$$
\mathbb P \left(\left|\sigma_\max(A) - \mathbb E[\sigma_\max(A)]\right| > t \right) \leq 2 \exp \left\\{-\frac{t^2}{2}\right\\} \implies \mathbb P \left(\left|\frac{\sigma_\max(A)}{\sqrt{n}} - \frac{\mathbb E[\sigma_\max(A)]}{\sqrt{n}}\right| > t \right) \leq \mathbb P \left(\left|\frac{\sigma_\max(A)}{\sqrt{n}} - \left(1 + \sqrt{\frac{p}{n}}\right)\right| > \tilde t \right) \leq 2 \exp \left\\{-\frac{nt^2}{2}\right\\},
$$
where, in the last inequality, we have defined $\tilde t = \frac{t}{\sqrt{n}}$.

This gives us an exponential bound we were looking for.
</div>















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