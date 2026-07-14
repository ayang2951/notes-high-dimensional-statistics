_Overview: This section covers large deviation principle._

## Overview of Large Deviation Theory

Large deviation theory aims to characterize the speed of convergence to 0 for large deviations from the mean. While the Central Limit Theorem utilizes the "correct" scaling to get, via convergence in distribution, the asymptotics for certain deviations from the mean, large deviation is interested in scalings <i>larger</i> than the one used in the CLT.

Suppose we have iid random variables $X_1,X_2,\ldots,X_n$ where $\mathbb E[X_i]=0$ and $\text{Var}(X_i)=1$.
Let
$$
S_n\triangleq \sum_{i=1}^n X_i.
$$
The central limit theorem gives us that, for fixed $t$,
$$
\mathbb P(S_n>\sqrt n t) \ \rightarrow \ \mathbb P(Z>t), \qquad Z\sim\mathcal N(0,1).
$$
This is essentially a question of how far the sum deviates from its expectation (which is 0, here) at the "sweet spot" scaling $\sqrt{n}$ the CLT gives us. We can, however, ask: what about for other scalings?

For example, if we scale at $n^{-\frac{1}{4}}$, we would have
$$
\mathbb P\left(S_n > n^{-\frac{1}{4}} \cdot t\right) \ = \ \mathbb P\left(\frac{1}{\sqrt{n}}S_n > \frac{t}{n^{\frac{1}{4}}}\right) \rightarrow \frac{1}{2}.
$$
This is a <i>small</i> deviation, i.e. with a scaling $o(\sqrt{n})$, and the result is easy via the CLT.

What about a <i>large</i> deviation, with a scaling $\Omega(\sqrt{n})$, such as $n$? We also know via the CLT that 
$$
\mathbb P(S_n>nt) = \mathbb P\left(\frac{S_n}{\sqrt n}>\sqrt n t\right)  \rightarrow 0.
$$
This is a large deviation, decaying to 0 instead of converging to something positive and fixed.

For more complicated functions than the sample sum, we might want to understand this decay: knowing that it converges to 0 is not enough; we want to know the rate at which is does. Large deviation theory was developed precisely to address this: providing more accurate asymptotics for scalings larger than the CLT.

<div class="callout remark"><span class="label">Remark: Concentration vs. Large Deviations</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Concentration of measure and large deviation theory both study tail probabilities, but they have slightly different goals. Concentration usually gives finite-sample upper bounds, often of the form
$$
\mathbb P(\text{bad event})\leq e^{-cn}.
$$
Large deviation theory is an asymptotic result, and instead tries to identify the exact exponential rate: when we take the log and divide by $n$, we should have results that look like
$$
\frac{1}{n}\log \mathbb P(\text{bad event})\rightarrow -\text{cost}.
$$
So it's important to keep in mind that large deviation results are asymptotic .
</div>

Let's now look at our favorite example: iid Gaussian random variables. This is a good first place to look, since the Gaussian tails are well-studied.


<div class="callout example"><span class="label">Example: Gaussian Large Deviation</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_1,X_2,\ldots,X_n\overset{iid}\sim \mathcal N(0,1)$. Then $S_n=\sum_{i=1}^n X_i\sim \mathcal N(0,n)$. Therefore, we have for the probability that the sum deviates from its mean for a large devitation is
$$
\mathbb P(S_n>nt) = \mathbb P(\sqrt n Z>nt) = \mathbb P(Z>\sqrt n t), \qquad Z\sim\mathcal N(0,1).
$$
Gaussian tail bounds give us that
$$
\mathbb P(Z>x)\sim \frac{1}{x\sqrt{2\pi}}e^{-x^2/2},
$$
and taking $x=\sqrt n t$, we get
$$
\mathbb P(S_n>nt) \sim \frac{1}{\sqrt n t\sqrt{2\pi}}e^{-nt^2/2}.
$$
Observe that there are two terms in the above: the exponential term in $n$, as well as the square root term in $n$. Usually, in large deviation theory, we do not care about the latter, the polynomial term. To remove these terms that don't matter to us, a common theme will be taking the log of the probability so that only the originally-exponential term remains. We then divide by $n$, the rate.

Doing so on the Gaussian form we obtained above,
$$
\frac{1}{n}\log\mathbb P(S_n>nt) \sim \frac{1}{n}\log\left(\frac{1}{\sqrt n t\sqrt{2\pi}}e^{-nt^2/2}\right) = -\frac{t^2}{2}-\frac{1}{2n}\log n-\frac{1}{n}\log(t\sqrt{2\pi}) \rightarrow -\frac{t^2}{2}.
$$
We see that only the originally exponential term (asymptotically) remains after taking the log.
</div>

Gaussian sums are some of the easiest to work with. Further work involves proving results for more difficult distributions and, eventually, more difficult functions.

To that end, let's take a meta-look at large deviation theory: large deviation theory and convergence in distribution results both deal with asymptotic shapes of distributions. Recall that, with convergence in distribution, we did not need to prove every result from the definition: establishing useful theorems such as the continuous mapping theorem, the delta method, and the method of moments gave us much more flexibility in analyzing convergence in distribution of different random variables and functions of random variables.

Similarly, we can prove more complicated and flexible large deviation results, analogous to those above, that will help us widen the scope of our results most efficiently.

### Cramér's Theorem

A first basic building block of large deviation, before we get to more complicated functions, is Cramér's theorem: a result that gives us the large deviation rate of convergence for sums of iid random variables.

<div class="callout theorem"><span class="label">Theorem: Cramer's Theorem</span><br/>
<hr style="height:0.01px; visibility:hidden;" />

Let $ X_1,\ldots,X_n $
be iid random variables satisfying
$$
\phi(t)\triangleq \mathbb E[e^{tX_1}]<\infty,
$$
and assume for simplicity that $\mathbb E[X_1]=0$.
Let
$$
S_n=\sum_{i=1}^n X_i.
$$
Then we have that (for large sufficiently large $a$)
$$
\lim_{n\rightarrow\infty}\frac{1}{n}\log\mathbb P(S_n\geq na)=-I(a),
$$
where
$$
I(a)=\sup_{t\in\mathbb R}\left\\{ta-\log\phi(t)\right\\}
$$
is called the rate function.
</div>

<details class="collapsible">
<summary>Proof of the upper bound.</summary>
<div class="collapsible__content">
Consider first $a > 0$, and let $\theta>0$ also. Then
$$
\mathbb P(S_n\geq na) = \mathbb P(e^{\theta S_n}\geq e^{n\theta a}) \overset{(1)}\leq e^{-n\theta a}\mathbb E[e^{\theta S_n}] = e^{-n\theta a}\mathbb E\left[e^{\theta\sum_{i=1}^n X_i}\right] = e^{-n\theta a}\prod_{i=1}^n\mathbb E[e^{\theta X_i}] = e^{-n\theta a}\phi(\theta)^n = \exp\left\{-n(\theta a-\log\phi(\theta))\right\}.
$$
Inequality $(1)$ is Markov's inequality applied to the nonnegative random variable $e^{\theta S_n}$, and by independence of the $X_i$ and algebraic manipuation, we end with the exponential form.

Since this holds for every $\theta>0$, we can optimize: we want the upper bound to be as small as posisble (to be tight as possible). Optimizing over $\theta$, we have
$$
\mathbb P(S_n\geq na) \leq \exp\left\\{-n\sup_{\theta>0}(\theta a-\log\phi(\theta))\right\\}.
$$
It remains to show that the sup is equivalent for $\theta \in \mathbb R^+$ and $\theta \in \mathbb R$.

<details class="collapsible">
<summary>Proof that the sup is equivalent over $\mathbb R^+$ and $\mathbb R$.</summary>
<div class="collapsible__content">
We have by Jensen's inequality, as $f(x) = \exp(x)$ is a convex function,
$$
     \Lambda(\theta) \triangleq \log \left(\mathbb E[e^{\theta X}] \right) \geq \log \left(\exp \{\theta \cdot \mathbb E[X]\}\right) = 0,   
$$
as $\mathbb E[X] = 0$. Hence, we have that $\Lambda$ is nonnegative. Here is where we need to use that $a > 0$: we have that
$$
\theta\cdot a - \Lambda(\theta) \leq \theta \cdot a \leq 0
$$ 
if $\theta \leq 0$, whereas if $\theta > 0$, we have
$$
\theta \cdot a - \Lambda(\theta) \leq \theta \cdot a \geq 0.
$$
Hence, the supremum over $\theta > 0$ is at least as large as the supremum over $\theta \leq 0$.
</div>
</details>

That the supremums are equivalent is because we are studying the upper tail. The lower tail, if we are interested in that, is the reverse: take $a < 0$ large and $\theta < 0$.
</div>
</details>

<!-- <div class="callout remark"><span class="label">Remark: The Lower Bound Heuristic</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
The upper bound comes from Markov's inequality and is relatively easy. The lower bound is more subtle.

The heuristic is that we want to tilt the distribution so that the rare event becomes typical. If there is a $\theta^\*$ such that
$$
\Lambda'(\theta^\*)=a,
$$
then we define a tilted distribution by
$$
\frac{d\mathbb P_{\theta^\*}}{d\mathbb P}(x)=\exp\{\theta^\*x-\Lambda(\theta^*)\}.
$$
Under this tilted distribution, the mean is $a$. Thus, $\bar X_n\approx a$ is no longer rare. Changing back to the original distribution produces the cost
$$
\exp\{-n(\theta^\*a-\Lambda(\theta^\*))\}=e^{-nI(a)}.
$$
</div> -->


### Rates and Rate Functions

We saw from the above that we have a particular scaling (faster than the CLT scaling) and our results are given in terms of a function $I(a)$, the rate function. 

In general, our results will look like
$$
\lim_{n \rightarrow \infty} \frac{1}{n} \log \mathbb P(Z_n \geq n \cdot a) = -I(a),
$$
where the scaling $\frac{1}{n}$ and $n$ in the probability are called the <i><strong>rate</strong></i> or speed parameter, and the function $I(a)$ is called the <i><strong>rate function</strong></i>.



<div class="callout remark"><span class="label">Remark: Log MGF Notation</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
It is common to define the log moment generating function
$$
\Lambda(t)\triangleq \log\phi(t)=\log\mathbb E[e^{tX_1}].
$$
Then the rate function is
$$
I(a)=\sup_{t\in\mathbb R}\{ta-\Lambda(t)\}.
$$
This is the Legendre transform, or convex conjugate, of $\Lambda$.
</div>


With our large deviation looking like
$$
\mathbb P(S_n \geq n a) = \mathbb P(\frac{1}{n}S_n \geq a) \sim \frac{e^{-nI(a)}}{\sqrt{n}},
$$
we can even think about the density: ignoring the square root factor and taking the derivative, we have that the density is (when ignoring the polynomial terms)
$$
f_{Z_n}(a) \sim I'(a) \cdot \exp \left\\{-nI(a)\right\\},
$$
so we get the same scaling.

Even more generally, if we do not think about the cumulative distribution function <i>or</i> the density but rather the law itself: for an arbitrary event $A$, we have
$$
\mathbb P(Z_n \in A)= \int_A f_{Z_n}(x)dx \sim \int_A e^{-nI(x)}dx.
$$
Here, we can invoke Laplace's method, and we get that this integral should be dominated by the point in $A$ where $I(x)$ is smallest:
$$
\int_A e^{-nI(x)}dx \approx e^{-n\inf_{x\in A}I(x)}.
$$
Therefore,
$$
\frac{1}{n}\log\mathbb P(\bar X_n\in A)\rightarrow -\inf_{x\in A}I(x).
$$
This is the large deviation principle in heuristic form: the probability of an event is dominated by the least costly way for the event to happen, where the density is the highest.

### The Large Deviation Principle

We now make the previous heuristic more formal. This first statement will still be a rough one, meant to evoke a specific image and intuition, but ignoring technical details.

A sequence of probability measures $\\{\mathbb P_n\\}$ satisfies a large deviation principle with rate function $I$ and rate $n$ if, roughly, for every event $A$, we have that
$$
\lim_{n\rightarrow\infty}\frac{1}{n}\log \mathbb P_n(A)=-\inf_{x\in A}I(x).
$$
There's still a technical issue here: we need to be careful with discontinuity points and boundaries. Convergence in distribution does not imply that $\mathbb P_n(A) \rightarrow \mathbb P(A)$ because of discontinuity points. Hence, we need to discuss upper and lower bounds with open and closed sets of A.

<div class="callout definition"><span class="label">Definition: Large Deviation Principle</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\{\mathbb P_n\}$ be a sequence of probability measures on a space $\mathcal X$. We say that $\{\mathbb P_n\}$ satisfies the <i><strong>large deviation principle</strong></i> with speed $n$ and rate function $I:\mathcal X\rightarrow[0,\infty]$ if:

$$
-\inf_{x\in A^\circ}I(x)\leq \liminf_{n\rightarrow\infty}\frac{1}{n}\log \mathbb P_n(A)\leq \limsup_{n\rightarrow\infty}\frac{1}{n}\log \mathbb P_n(A)\leq -\inf_{x\in \bar A}I(x),
$$
where $A^\circ$ is the interior of the set and $\bar A$ is the closure.
<!-- <ol type="1">
  <li>$I$ is lower semicontinuous.</li>
  <li>For every closed set $F\subseteq\mathcal X$,
  $$
  \limsup_{n\rightarrow\infty}\frac{1}{n}\log P_n(F)\leq -\inf_{x\in F}I(x).
  $$</li>
  <li>For every open set $G\subseteq\mathcal X$,
  $$
  \liminf_{n\rightarrow\infty}\frac{1}{n}\log P_n(G)\geq -\inf_{x\in G}I(x).
  $$</li>
</ol> -->
</div>

Note that if
$$
\inf_{x\in A^\circ}I(x)=\inf_{x\in \overline A}I(x),
$$
then the limit exists and
$$
\lim_{n\rightarrow\infty}\frac{1}{n}\log P_n(A)=-\inf_{x\in A}I(x),
$$
so our "heuristic definition" will naturally hold.

<!-- <div class="callout remark"><span class="label">Remark: Why Open and Closed Sets Appear</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
The lower bound is stated for open sets because, to prove that an event has at least a certain probability, we need the event to contain a neighborhood of a point. The upper bound is stated for closed sets because, to prove that an event has at most a certain probability, we need to control all possible limiting points. The issue is boundary behavior. If the boundary of $A$ has a different exponential cost from the interior, then a single clean formula for $A$ can fail.
</div> -->

With our formal large deviation principle statement, we can rewrite Cramér's theorem in the following way. 

<div class="callout theorem"><span class="label">Theorem: Cramér's Theorem, Large Deviation Principle</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $X_1,X_2,\ldots$ be iid random variables with $\mathbb E[X_1]=0$
and moment generating function $\phi(t)=\mathbb E[e^{tX_1}]$, where $\phi(t)<\infty$. Then the distribution of
$$
\bar X_n=\frac{1}{n}\sum_{i=1}^n X_i
$$
satisfies an LDP with speed $n$ and rate function
$$
I(x)=\sup_{t\in\mathbb R}\{tx-\log\phi(t)\}.
$$
</div>

In particular, for a nice set $A$,
$$
\mathbb P(\bar X_n\in A)\approx \exp\left\\{-n\inf_{x\in A}I(x)\right\\},
$$
from our previous heuristic calculations.

<!-- ## Interpreting the Rate Function

The rate function
$$
I(x)=\sup_{t\in\mathbb R}\{tx-\log\phi(t)\}
$$
is the cost of forcing the empirical mean to be near $x$.

A few basic facts are useful.

<div class="callout proposition"><span class="label">Proposition: Basic Properties of the Cramér Rate Function</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let
$$
I(x)=\sup_{t\in\mathbb R}\{tx-\Lambda(t)\},\qquad \Lambda(t)=\log\mathbb E[e^{tX}].
$$
Then:
<ol type="1">
  <li>$I(x)\geq 0$.</li>
  <li>If $\mu=\mathbb E[X]$, then $I(\mu)=0$.</li>
  <li>$I$ is convex, because it is a supremum of affine functions of $x$.</li>
  <li>Larger values of $I(x)$ correspond to rarer events.</li>
</ol>
</div>

To check that $I(\mu)=0$, note that
$$
I(\mu) \overset{(1)}= \sup_{t\in\mathbb R}\{t\mu-\Lambda(t)\} \overset{(2)}\geq 0.
$$
Equality $(1)$ is the definition of $I$. Inequality $(2)$ is because taking $t=0$ gives $t\mu-\Lambda(t)=0-\log 1=0$.

On the other hand, Jensen's inequality gives
$$
\Lambda(t)=\log\mathbb E[e^{tX}] \overset{(1)}\geq \log e^{t\mathbb E[X]} \overset{(2)}= t\mu.
$$
Inequality $(1)$ is Jensen's inequality applied to the convex function $e^x$, equivalently $\mathbb E[e^{tX}]\geq e^{t\mathbb E[X]}$. Equality $(2)$ is because $\mu=\mathbb E[X]$.

Thus $t\mu-\Lambda(t)\leq 0$ for every $t$, and so $I(\mu)=0$.

For example, if $X\sim\mathcal N(0,1)$, then
$$
\phi(t)=e^{t^2/2},\qquad \Lambda(t)=\frac{t^2}{2}.
$$
Thus,
$$
I(x)=\sup_{t\in\mathbb R}\left\{tx-\frac{t^2}{2}\right\}.
$$
The maximizing value satisfies
$$
x-t=0,
$$
so $t=x$. Hence,
$$
I(x)=x^2-\frac{x^2}{2}=\frac{x^2}{2}.
$$
This agrees with the direct Gaussian calculation. -->

## Tools for Building Large Deviation Principles

Cramér's theorem is one way to prove an LDP, but it is specialized to sums of iid random variables. Here, as promised, we introduce tools for building LDPs in more complicated settings.

As stated before, our goal will be to construct tools analogous to those commonly used in convergence in distribution results:
<ol type="1">
  <li>Gärtner-Ellis theorem: an analogue of using MGFs.</li>
  <li>Varadhan's theorem and Bryc's theorem: analogues of Portmanteau.</li>
  <li>Contraction principle: an analogue of the continuous mapping theorem.</li>
</ol>

### Gärtner-Ellis Theorem

First, we shall introduce the Gartner-Ellis theorem, which is the analog of the method of moments for convergence in distribution: under some regularity conditions (such as the moments uniquely characterizing the distribution),
$$
\mathbb E[e^{tX_n}]\rightarrow \mathbb E[e^{tX}] \quad \implies \quad X_n\overset{d}\rightarrow X.
$$
Suppose that the sequence of random variables $Z_1, Z_2 \ldots Z_n$ satisfies the large deviation principle with rate $n$ and rate function $I(x)$. Then we have that the density heuristic
$$
f_n(x) \sim e^{-nI(x)}
$$
gives us that, for the MGF $\phi(t)$ of $Z_n$,
$$
\mathbb E[e^{tZ_n}] = \int e^{tx}f_n(x)dx \sim \int e^{tx}e^{-nI(x)}dx.
$$
If we use the Laplace method on this directly, we see that the term involving $g(x) \triangleq e^{tx}$ disappears: we would have that
$$
\int g(x) e^{-nI(x)}dx \rightarrow e^{tx_0} \cdot e^{-nI(x_0)},
$$
where $x_0$ is the minimizer of the rate function $I(x)$. Hence, this scaling is "off", in that it's not useful for our analysis when using the moment generating function: after mapping by $\frac{1}{n}\log$, the $tx$ term disappears.

Because the MGF is precisely what we want to analyze, we naturally need that the term doesn't disppear. We deduce that the scaling of the exponential term needs to have $n$ as well if we don't want it to disappear: from that, we would obtain via Laplace's method (assuming that the density is nice enough, such as that it has a unique minimizer and that the density around the minimum is well-approximated by a quadratic):
$$
\mathbb E[e^{ntZ_n}] = \int e^{ntx}f_n(x)dx \sim \int e^{ntx}e^{-nI(x)}dx = \int e^{n(tx-I(x))}dx.
$$
Applying Laplace's method on the scaled version of the MGF is much more fruitful, as we see from the following:
$$
\frac{1}{n}\log\phi_n(nt) \approx \sup_x\{tx-I(x)\}.
$$
We remark that the right hand side here is another Legendre transform, this time of the rate itself. 

With the understanding of what scaling the MGF needs in order to be "visible" in the end result, we introduce the quintessential MGF result on large deviation; the analog of the method of moments, called the Gärtner-Ellis theorem.

<div class="callout theorem"><span class="label">Theorem: Gärtner-Ellis Theorem</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\{\mathbb P_n\}$ be a sequence of probability measures, and let $\phi_n(t)\triangleq \mathbb E[e^{ntZ_n}]$. Additionally, define the scaled log of the MGF as
$$
\Lambda_n(t)\triangleq \frac{1}{n}\log\phi_n(t)=\frac{1}{n}\log\mathbb E[e^{ntZ_n}].
$$
Suppose that $\Lambda(t)=\lim_{n\rightarrow\infty}\Lambda_n(t)$ exists and is differentiable. 

Then $\\{\mathbb P_n\\}$ satisfies an LDP with rate $n$ and rate function
$$
I(z)=\sup_{t\in\mathbb R}\{tz-\Lambda(t)\}.
$$
This result is known as the <i><strong>Gärtner-Ellis</strong></i> theorem.
</div>

<div class="callout remark"><span class="label">Remark: Be Careful About the Conditions</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Differentiability is not an easily-removed assumption, and is in fact quite important. 

<!-- Intuitively, nondifferentiability means that there may be more than one competing exponential mechanism. -->
</div>

Now, we introduce an example of where Gärtner-Ellis fails to apply, in consideration of the remark above.

<div class="callout example"><span class="label">Example: Gartner-Ellis on Gaussian Mixture</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Consider the Gaussian mixture distribution
$$
X_n\sim \frac{1}{2}\mathcal N\left(1,\frac{1}{n}\right)+\frac{1}{2}\mathcal N\left(-1,\frac{1}{n}\right).
$$
We shall show that
<ol type="a">
  <li>$X_n$ satisfies the LDP from the definition, and $I(x)$ is not convex.</li>
  <li>Gärtner-Ellis gives us a different function than $I(x)$.</li>
</ol>

First, from the density, the probability of seeing $X_n$ near $x$ is dominated by whichever Gaussian component is closer to $x$. The density is proportional to
$$
e^{-\frac{n}{2}(x-1)^2}+e^{-\frac{n}{2}(x+1)^2}.
$$
Therefore, the true rate function is
$$
I(x)=\min\left\\{\frac{1}{2}(x-1)^2,\frac{1}{2}(x+1)^2\right\\}.
$$
This has two wells, one at $x=1$ and one at $x=-1$, and is not convex.

Now compute the scaled log MGF. If $Y\sim \mathcal N(\mu,1/n)$, then
$$
\mathbb E[e^{ntY}]\overset{(1)}= \exp\left\\{nt\mu+\frac{1}{2}n^2t^2\cdot\frac{1}{n}\right\\}\overset{(2)}=e^{n(t\mu+t^2/2)}.
$$
Equality $(1)$ is the MGF formula for a Gaussian. Equality $(2)$ is simplifying the variance term.

Thus,
$$
\mathbb E[e^{ntX_n}] \overset{(1)}= \frac{1}{2}e^{n(t+t^2/2)}+\frac{1}{2}e^{n(-t+t^2/2)} \overset{(2)}= e^{nt^2/2}\cdot \frac{e^{nt}+e^{-nt}}{2}.
$$
Equality $(1)$ is applying the Gaussian MGF formula to the two mixture components. Equality $(2)$ factors out $e^{nt^2/2}$.

Therefore,
$$
\Lambda(t) \overset{(1)}= \lim_{n\rightarrow\infty}\frac{1}{n}\log\mathbb E[e^{ntX_n}] \overset{(2)}= \lim_{n\rightarrow\infty}\left[\frac{t^2}{2}+\frac{1}{n}\log\left(\frac{e^{nt}+e^{-nt}}{2}\right)\right] \overset{(3)}= \frac{t^2}{2}+|t|.
$$
Equality $(1)$ is the definition of the limiting scaled log MGF. Equality $(2)$ uses the expression above. Equality $(3)$ is because $\frac{1}{n}\log(e^{nt}+e^{-nt})\rightarrow |t|$, while $\frac{1}{n}\log 2\rightarrow 0$.

This function is not differentiable at $t=0$.

If we take the Legendre transform anyway, we get
$$
\Lambda^\*(x)=\sup\_{t\in\mathbb R}\{tx-\Lambda(t)\}.
$$
This produces the convexified rate
$$
\Lambda^\*(x)=\begin{cases}\frac{1}{2}(x-1)^2, & x>1,\\ 0, & -1\leq x\leq 1,\\ \frac{1}{2}(x+1)^2, & x<-1.\end{cases}
$$
This differs from the true rate function.
</div>

<div class="callout remark"><span class="label">Filling in the Gap: Why This Does Not Contradict Part (a)</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
<span style="color:#b00020;">The handwritten notes ask why part (b) contradicts part (a). The resolution is that it does not contradict part (a), because the differentiability hypothesis in Gärtner-Ellis fails. The limiting scaled log MGF is
$$
\Lambda(t)=\frac{t^2}{2}+|t|,
$$
which is not differentiable at $t=0$. Therefore the simplified Gärtner-Ellis theorem is not applicable.</span>
</div>

This example is useful because it shows that true rate functions need not be convex, while Legendre transforms always produce convex functions. Moment-generating-function methods naturally see the convex envelope of the behavior.

### Varadhan's Theorem and Bryc's Theorem

Suppose we're interested in using convergence in distribution as a model for our LDP results. We have that, by the Portmanteau theorem,
$$
X_n\overset{d}\rightarrow X \quad \iff \quad \mathbb E[g(X_n)]\rightarrow \mathbb E[g(X)]
$$
for every bounded continuous function $g\in C_b$.

It turns out that the Vardhan-Bryc theorem is the analog of Portmanteau for LDP: Varadhan's theorem is one direction (forward), and Bryc's theorem is the other (backward).

Suppose $X_n\sim \mathbb P_n$ and $\\{\mathbb P_n\\}$ satisfies an LDP with speed $n$ and rate function $I(x)$. If we look at
$$
\mathbb E[g(X_n)]=\int g(x)d\mathbb P_n(x),
$$
then using the heuristic $d\mathbb P_n(x)\approx e^{-nI(x)}dx$ gives
$$
\mathbb E[g(X_n)] = \int g(x)d\mathbb P_n(x) \approx \int g(x)e^{-nI(x)}dx.
$$
By Laplace's method, this is dominated by the minimizer of $I$, and the function $g$ does not really show up at exponential scale:
$$
\frac{1}{n}\log\mathbb E[g(X_n)]\rightarrow -\inf_x I(x).
$$
So a function $g(X_n)$ is again not very interesting without the correct scaling.

If the function is scaled exponentially, however, we would get that
$$
\mathbb E\left[e^{nF(X_n)}\right] = \int e^{nF(x)}d\mathbb P_n(x) \approx \int e^{nF(x)}e^{-nI(x)}dx = \int e^{n(F(x)-I(x))}dx.
$$
By Laplace's method, we should find what the part that dominates this integral: this should be where $F(x) - I(x)$ is maximized, i.e. $\sup_x F(x) - I(x)$, and we would obtain
$$
\int e^{n(F(x)-I(x))}dx \approx \exp\left\\{n\sup_x(F(x)-I(x))\right\\}.
$$

<div class="callout theorem"><span class="label">Theorem: Varadhan's Theorem</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\mathbb P_n$ satisfy an LDP with rate $n$ and rate function $I(x)$. Let $F$ be a bounded continuous function. Then
$$
\lim_{n\rightarrow\infty}\frac{1}{n}\log\int e^{nF(x)}d\mathbb P_n(x)=\sup_x\{F(x)-I(x)\}.
$$
</div>

## Bryc's Theorem

Bryc's theorem is the reverse direction of Varadhan's theorem, but it is far less useful and used in practice.

We give a preliminary definition, then state the theorem.

<div class="callout definition"><span class="label">Definition: Exponential Tightness</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
A sequence of probability measures $\{\mathbb P_n\}$ is <i><strong>exponentially tight</strong></i> with speed $n$ if, for every $L<\infty$, there exists a compact set $K_L$ such that
$$
\limsup_{n\rightarrow\infty}\frac{1}{n}\log P_n(K_L^c)\leq -L.
$$
</div>

<div class="callout theorem"><span class="label">Theorem: Bryc's Theorem</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
For a bounded continuous function $F$, define
$$
\Lambda_n(F)\triangleq \frac{1}{n}\log\int e^{nF(x)}d\mathbb P_n(x).
$$
If $\{\mathbb P_n\}$ is exponentially tight and $\Lambda(F)=\lim_{n\rightarrow\infty}\Lambda_n(F)$ exists for all $F\in C_b$ (the class of continuous and bounded functions), then $\mathbb P_n$ satisfies an LDP with speed $n$ and rate function
$$
I(x)=\sup_{F\in C_b}\{F(x)-\Lambda(F)\}.
$$
</div>

Bryc's theorem is clean and conceptually useful, but it is not seen as much in practice because the rate function is hard to compute: optimizing the quantity
$$
I(x)=\sup_{F\in C_b}\{F(x)-\Lambda(F)\}.
$$
over all bounded continuous functions is usually not possible. Analogously, for the convergence in distribution results, we don't usually <i>literally</i> show that all continuous and bounded functions converge, but rather a useful subset, such as the moment generating function or characteristic function.

<div class="callout remark"><span class="label">Remark: Varadhan vs. Bryc</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Varadhan's theorem says:
$$
\text{LDP}\quad\Longrightarrow\quad\text{asymptotics of exponential integrals}.
$$
Bryc's theorem says, under exponential tightness:
$$
\text{asymptotics of exponential integrals}\quad\Longrightarrow\quad\text{LDP}.
$$
</div>

### Contraction Mapping Theorem

The next tool is the large-deviation version of the continuous mapping theorem, called the contraction mapping theorem.

Recall the continuous mapping theorem: if $X_n\overset{d}\rightarrow X$ and $g$ is continuous, then
$$
g(X_n)\overset{d}\rightarrow g(X).
$$
We now show the large deviation analog of the continuous mapping theorem. In essence, with a sequence $X_n$ that satisfies a large deviation principle with some rate function and a continuous map $T$, we want to derive the rate function of the sequence $Y_n = T(X_n)$.

Suppose $X_n\sim \mathbb P_n$ and that we have a continuous mapping
$$
T:\mathcal X\rightarrow\mathcal Y.
$$
For the random variable
$$
Y_n \triangleq T(X_n),
$$
let $\mathbb Q_n$ denote the distribution of $Y_n$.

Assume that $\\{\mathbb P_n\\}$ satisfies an LDP with speed $n$ and rate function $I(x)$. We ask the following questions:
<ol type="1">
  <li>Does $\mathbb Q_n$ satisfy an LDP?</li>
  <li>If yes, what are the rate and the rate function of $\mathbb Q_n$?</li>
</ol>

To build intuition, assume for a moment (ha) that the density of $X_n$ has that $X_n \sim e^{-nI(x)}$.

Further suppose for now that $\mathbb Q_n$ satisfies LDP with rate function $J(y)$. We want to understand the probability that $Y_n$ is near $y$.

If we have that
$$
\mathbb Q_n([y, y + \Delta y]) \approx e^{-nJ(y)} \Delta y,
$$
and call the event $A_y \triangleq T^{-1}(y)$, then write everything in terms of $x$, we have
$$
\int_{A_y} e^{-n(x)}dx \cdot \Delta x \approx e^{-\inf_{x \in A_y} I(x)},
$$
which implies that the rate function is $J(y) = \inf_{x \in A_y} I(x)$.

<div class="callout remark"><span class="label">More Rigorous Argument Using the Pushforward Measure</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
To do this more rigorously, we can avoid using densities and instead use the pushforward measure: we know that $\mathbb Q_n(B) = \mathbb P_n\big(T^{-1}(B)\big)$. Since we have the LDP for $\mathbb P_n$, we can apply it to the preimage $T^{-1}(B)$. Note that continuity of $T$ ensures that for open or closed $B$, $T^{-1}(B)$ remains open and closed, respectively.
</div>

Hence, we give the formal theorem for the contraction mapping theorem.

<div class="callout theorem"><span class="label">Theorem: Contraction Mapping Theorem</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Let $\{\mathbb P_n\}$ be a sequence of probability measures satisfying an LDP with rate $n$ and rate function $I(x)$. Let
$$
T:\mathcal X\rightarrow\mathcal Y
$$
be continuous, and let $\mathbb Q_n$ denote the pushforward measures on $\mathcal Y$, $\mathbb  Q_n=P_n\circ T^{-1}$.

Then $\mathbb Q_n$ satisfies an LDP with speed $n$ and rate function
$$
J(y)=\inf_{x:T(x)=y}I(x).
$$
This is known as the <i><strong>contraction mapping theorem</strong></i>.
</div>

<!-- <div class="callout example"><span class="label">Example: Squaring a Large-Deviation Sequence</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Suppose $X_n$ satisfies an LDP on $\mathbb R$ with rate function $I$, and define
$$
Y_n=X_n^2.
$$
Here
$$
T(x)=x^2.
$$
The contraction principle says that $Y_n$ satisfies an LDP with rate function
$$
J(y)=\inf_{x:x^2=y}I(x).
$$
Therefore,
$$
J(y)=\begin{cases}\min\{I(\sqrt y),I(-\sqrt y)\}, & y\geq 0,\\ \infty, & y<0.\end{cases}
$$
</div> -->

<!-- ## Summary

The main theme in both Laplace's method and large deviation theory is that exponential scaling turns analysis into optimization.

Laplace's method says that
$$
\int e^{nh(x)}dx
$$
is dominated by the point where $h$ is largest, while
$$
\int e^{-nh(x)}dx
$$
is dominated by the point where $h$ is smallest.

Large deviation theory says that
$$
\mathbb P(Z_n\in A)\approx \exp\left\{-n\inf_{x\in A}I(x)\right\}.
$$
So rare events happen through the least costly mechanism.

The tools fit together as follows:
<ol type="1">
  <li><strong>Laplace's method:</strong> Exponential integrals are dominated by maxima or minima.</li>
  <li><strong>Cramér's theorem:</strong> Empirical means of iid random variables satisfy an LDP with rate function given by the Legendre transform of the log MGF.</li>
  <li><strong>Gärtner-Ellis theorem:</strong> More general sequences can satisfy LDPs if their scaled log MGFs converge nicely.</li>
  <li><strong>Varadhan's theorem:</strong> An LDP implies asymptotics for exponential integrals.</li>
  <li><strong>Bryc's theorem:</strong> Under exponential tightness, enough exponential integral asymptotics imply an LDP.</li>
  <li><strong>Contraction principle:</strong> Continuous maps preserve LDPs, with the new rate function obtained by minimizing over preimages.</li>
</ol>

The recurring slogan is:
$$
\text{at exponential scale, the dominant term wins.}
$$ -->