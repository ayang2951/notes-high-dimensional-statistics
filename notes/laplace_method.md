_Overview: This section covers the laplace method, an approximation for hard-to-compute integrals that will be important in the discussion on large deviation principle and, subsequently, the replica method._

## Laplace's Method

Laplace's method is a clever way to approximate complex integrals. The central theme is to isolate the parts of the integral that dominate the overall behavior, then do a polynomial approximation in that region.

There are many different variants of Laplace's method, most of which we do not have time to cover. We'll introduce a few variants that are relatively simple but quite useful.

Suppose we have an integral of the form
$$
\int_{-\infty}^\infty e^{-n h(x)} dx
$$
for a function $h : \mathbb R \rightarrow \mathbb R^+$ and $n$ that is large. If $h$ is complicated, this integral may be hard to compute. Take note of something important, though: since the integrand is exponential and the exponent is further scaled by some large $n$, it should be possible to approximate the integral very accurately by understanding where $h$ is minimized, as the integral is only non-negligible when $h$ is very small.

The key intuition that powers the Laplace method is thus: where the function $h(x)$ takes values even slightly larger than its minimum, $e^{-n h(x)}$ is <i>exponentially</i> smaller. Hence, for large $n$, most of the contribution comes from a very small neighborhood around the minimizer of $h$.

For now, suppose that $h$ has a unique minimizer at the point $a$ and that $h(a) = 0$. Then $e^{-n h(a)} = 1$, while $e^{-n h(x)}$ is tiny when $h(x)$ is bounded away from $0$.

Let $\delta > 0$. We will consider a $\delta$ neighborhood around $a$, which should be the only area that ultimately matters for the integral. Let's look at one side first: we have
$$
\int_{a + \delta}^\infty e^{-n h(x)} dx = \int_{a + \delta}^\infty e^{-(n - 1) h(x)} \cdot e^{-h(x)} dx \overset{(1)}\leq \int_{a + \delta}^\infty e^{-(n - 1)h(a + \delta)} e^{-h(x)}dx = e^{-n h(a + \delta)} \cdot e^{h(a + \delta)} \cdot \int_{a + \delta}^\infty e^{-h(x)}dx \overset{(2)}= O(e^{-cn}).
$$
Inequality $(1)$ holds under the assumption that $h(x) \geq h(a+\delta) > 0$ for all $x \in [a + \delta, \infty)$, and hence $e^{-(n-1)h(x)} \leq e^{-(n-1)h(a+\delta)}$. $(2)$ is because the remaining integral is constant in $n$, while the term outside the integral decays exponentially in $n$.

More generally, if $h$ is not monotone on the right side, we can replace $h(a+\delta)$ with the right-tail infimum
$$
c_{\delta,+} \triangleq \inf_{x\geq a+\delta}h(x).
$$
As long as $c_{\delta,+}>0$, we get
$$
\int_{a+\delta}^{\infty}e^{-(n-1)h(x)}e^{-h(x)}dx \overset{(1)}\leq e^{-(n-1)c_{\delta,+}}\int_{a+\delta}^{\infty}e^{-h(x)}dx \overset{(2)}= O(e^{-c_{\delta,+}n}).
$$
Here, $(1)$ is because $h(x)\geq c_{\delta,+}$ for all $x\geq a+\delta$, and equality $(2)$ gives us the same order of bound as when $h$ is monotone.

The same holds for the left side $(-\infty,a-\delta]$ by symmetry. If
$$
c_{\delta,-} \triangleq \inf_{x\leq a-\delta}h(x)>0,
$$
then
$$
\int_{-\infty}^{a-\delta}e^{-nh(x)}dx = O(e^{-c_{\delta,-}n}).
$$
From this, we see that exponentially small terms should not matter for $n$ large when $h(x)$ is not near 0, i.e. when $x$ is near $a$. On the other hand, the contribution near $a$ is <i>not</i> exponentially small, but rather <i>polynomially</i> small, and hence matters.

To show that the neighborhood around $a$ is of an order more important than the rest of the integral, we'll prove that
$$
\int_{a-\delta}^{a+\delta}e^{-nh(x)}dx = \Omega\left(\frac{1}{\sqrt n}\right).
$$
To facilitate using a Taylor expansion, let us assume that $h$ is three times differentiable.

Taylor's theorem gives, because $h(a) = h'(a) = 0$,
$$
h(x) \ = \ h(a)+h'(a)(x-a)+\frac{1}{2}h''(a)(x-a)^2+O((x-a)^3) \ = \ \frac{1}{2}h''(a)(x-a)^2+O((x-a)^3).
$$
So, near $a$, the integral behaves like
$$
\int e^{-nh(x)}dx \approx \int \exp\left\\{-\frac{n}{2}h''(a)(x-a)^2\right\\}dx.
$$
This is a Gaussian integral, and we see that the variance is $\frac{1}{nh''(a)}$, which should give us a hint as the to the correct order.

<!-- This already explains why the answer should be of order $\frac{1}{\sqrt n}$. For example, if $h(x)\leq C(x-a)^2$ on a small neighborhood of $a$, then
$$
\int_{a-\frac{1}{\sqrt n}}^{a+\frac{1}{\sqrt n}}e^{-nh(x)}dx \overset{(1)}\geq \int_{a-\frac{1}{\sqrt n}}^{a+\frac{1}{\sqrt n}}e^{-Cn(x-a)^2}dx \overset{(2)}= \frac{1}{\sqrt n}\int_{-1}^{1}e^{-Cu^2}du \overset{(3)}= \Omega\left(\frac{1}{\sqrt n}\right).
$$
Inequality $(1)$ is because $h(x)\leq C(x-a)^2$ near $a$, so $-nh(x)\geq -Cn(x-a)^2$. Equality $(2)$ is the change of variables $u=\sqrt n(x-a)$. Equality $(3)$ is because $\int_{-1}^{1}e^{-Cu^2}du$ is a positive constant independent of $n$. -->

We also notice something interesting: we can probably scale the domain of integration <i>with</i> the sample size $n$: if we shrink the domain at the correct speed, the part of the integral we lose in the approximation can disappear faster than it matters. 

For example, let's choose
$$
\delta_n = n^{-2/5}.
$$
Note that this, while it obviously converges to 0, scales less slowly than $\frac{1}{\sqrt{n}}$, i.e. $\delta_n = \omega\left(\frac{1}{\sqrt{n}}\right)$.

On the interval $|x-a|\leq \delta_n$, Taylor's theorem gives
$$
\left|h(x) - h(a) - h'(a)(x - a) - \frac{1}{2}h''(a)(x - a)^2 \right| = O\left(h^{(3)}(a) \cdot \left(n^{-\frac{2}{5}}\right)^3\right) = O\left(h^{(3)}(a) \cdot n^{-\frac{6}{5}}\right)
$$
Multiplying by $n$, we get
$$
nh(x) = \frac{n}{2}h''(a)(x-a)^2+O(n^{-1/5}).
$$
The scaling we chose, then, is a reasonable one, since, after multiplying by $n$, we still have that the error is $o(1)$.
Therefore, on this shrinking interval,
$$
e^{-nh(x)} \approx \exp\left\\{-\frac{n}{2}h''(a)(x-a)^2\right\\}.
$$
So the main part of the integral should satisfy
$$
\int_{a-\delta_n}^{a+\delta_n}e^{-nh(x)}dx \approx \int_{a-\delta_n}^{a+\delta_n}\exp\left\\{-\frac{n}{2}h''(a)(x-a)^2\right\\}dx.
$$
We have a choice here: either we can try to use Gaussian integrals, using the inverse variance and Gaussian tail bounds, to compute it exactly, or we could do some more clever approximation: use the change of variables $u=\sqrt{nh''(a)}(x-a)$. Then
$$
\int_{a-\delta_n}^{a+\delta_n}\exp\left\\{-\frac{n}{2}h''(a)(x-a)^2\right\\}dx \ = \ \frac{1}{\sqrt{nh''(a)}}\int_{-\sqrt{nh''(a)}\delta_n}^{\sqrt{nh''(a)}\delta_n}e^{-u^2/2}du \ \overset{(1)}\sim \ \frac{1}{\sqrt{nh''(a)}}\int_{-\infty}^{\infty}e^{-u^2/2}du = \sqrt{\frac{2\pi}{nh''(a)}}.
$$
Asymptotic equivalence $(1)$ is because $\sqrt{n}\delta_n=n^{1/10}\rightarrow\infty$, so the limits of integration go to $\pm\infty$. Here is where we see that the specific scaling we chose was a good choice.

So, from this illustrative example, we see that the overall idea of the Laplace method is
$$
\int_{-\infty}^{\infty}e^{-nh(x)}dx \approx \int_{a-\delta_n}^{a+\delta_n}e^{-nh(x)}dx \approx \int_{a-\delta_n}^{a+\delta_n}e^{-\frac{n}{2}h''(a)(x-a)^2}dx \approx \int_{-\infty}^{\infty}e^{-\frac{n}{2}h''(a)(x-a)^2}dx = \sqrt{\frac{2\pi}{nh''(a)}}.
$$
The first approximation is because the nonlocal tails are exponentially small. The second approximation is because Taylor expansion gives a quadratic near $a$. The third approximation is because the shrinking interval is still much wider than the Gaussian scale $n^{-1/2}$.

<div class="callout theorem"><span class="label">Theorem: Laplace's Method</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
Suppose $h:\mathbb R\rightarrow\mathbb R^+$ is three times continuously differentiable near $0$, and suppose
$$
h(0)=0,\qquad h'(0)=0,\qquad h''(0)\geq c >0.
$$
Assume also that $0$ is the unique minimizer and that the tails are sufficiently well-behaved so that the contribution away from $0$ is exponentially small. Then
$$
\lim_{n\rightarrow\infty}\frac{\int_{-\infty}^{\infty}e^{-nh(x)}dx}{\sqrt{\frac{2\pi}{nh''(0)}}}=1.
$$
Equivalently,
$$
\int_{-\infty}^{\infty}e^{-nh(x)}dx \sim \sqrt{\frac{2\pi}{nh''(0)}}.
$$
</div>



<details class="collapsible">
<summary>Proof sketch.</summary>
<div class="collapsible__content">
Because $h(0)=h'(0)=0$, we have by Taylor's theorem that
$$
\lim_{x\rightarrow 0}\frac{h(x)-\frac{1}{2}h''(0)x^2}{x^2}=0.
$$
Thus, for every $\epsilon>0$, there exists $\delta>0$ such that, for all $|x|<\delta$,
$$
\left|h(x)-\frac{1}{2}h''(0)x^2\right|\leq \epsilon x^2.
$$
Equivalently,
$$
\left(\frac{1}{2}h''(0)-\epsilon\right)x^2 \leq h(x) \leq \left(\frac{1}{2}h''(0)+\epsilon\right)x^2.
$$
Therefore, on the interval $[-\delta,\delta]$, we have
$$
\int_{-\delta}^{\delta}e^{-n(\frac{1}{2}h''(0)+\epsilon)x^2}dx \leq \int_{-\delta}^{\delta}e^{-nh(x)}dx \leq \int_{-\delta}^{\delta}e^{-n(\frac{1}{2}h''(0)-\epsilon)x^2}dx.
$$
For any fixed $\alpha>0$,
$$
\int_{-\delta}^{\delta}e^{-n\alpha x^2}dx \overset{(1)}= \frac{1}{\sqrt n}\int_{-\delta\sqrt n}^{\delta\sqrt n}e^{-\alpha u^2}du \overset{(2)}\sim \frac{1}{\sqrt n}\int_{-\infty}^{\infty}e^{-\alpha u^2}du \overset{(3)}= \sqrt{\frac{\pi}{n\alpha}}.
$$
Equality $(1)$ is the change of variables $u=\sqrt n x$. Asymptotic equivalence $(2)$ is because $\delta\sqrt n\rightarrow\infty$. Equality $(3)$ is the Gaussian integral.

Applying this with $\alpha=\frac{1}{2}h''(0)+\epsilon$ and $\alpha=\frac{1}{2}h''(0)-\epsilon$, and then letting $\epsilon\downarrow 0$, gives
$$
\int_{-\delta}^{\delta}e^{-nh(x)}dx \sim \sqrt{\frac{2\pi}{nh''(0)}}.
$$
Finally, the contribution outside $[-\delta,\delta]$ is exponentially small by the tail argument above, and hence it is negligible compared to $n^{-1/2}$.
</div>
</details>

<div class="callout remark"><span class="label">Remark: Recentered Laplace's Method</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
If the unique minimizer is $a$ instead of $0$, and if $h(a)=0$, then translating the variable gives
$$
\int_{-\infty}^{\infty}e^{-nh(x)}dx \sim \sqrt{\frac{2\pi}{nh''(a)}}.
$$
</div>

In the next few sections, we will look at some variants of Laplace's method for simple variants of the form above.

### Laplace's Method with a Prefactor

A useful variant of Laplace's method is when the integral has a multiplicative factor: for a thrice differentiable function $h$ where $h(0)=h'(0)=0$, $h''(0)>0$, and a continuous function $g$ where $g(0)\neq 0$,
$$
\int_{-\infty}^{\infty}g(x)e^{-nh(x)}dx.
$$
Since the non-negligible part of the integral is concentrated near $0$ and $g$ is continuous, we should be able to replace $g(x)$ by $g(0)$. Hence,
$$
\int_{-\infty}^{\infty}g(x)e^{-nh(x)}dx \sim g(0)\sqrt{\frac{2\pi}{nh''(0)}}.
$$
To show this, we split the domain of integration again. For the right tail,
$$
\int_{\delta}^{\infty}g(x)e^{-nh(x)}dx \ \leq \ e^{-(n-1)h(\delta)}\int_{\delta}^{\infty}|g(x)|e^{-h(x)}dx = O(e^{-cn}).
$$
In the region of the domain where the integral matters, we have that $g(x)\approx g(0)$ and $h(x)\approx \frac{1}{2}h''(0)x^2$. Therefore,
$$
\int_{-\delta}^{\delta}g(x)e^{-nh(x)}dx \approx \int_{-\delta}^{\delta}g(0)e^{-\frac{n}{2}h''(0)x^2}dx \approx g(0)\sqrt{\frac{2\pi}{nh''(0)}}.
$$
We see that, if $g$ doesn't scale similarly to the exponential term, the effect of $g$ is only visible through the constant.

### Laplace's Method on Functions Bounded Away from 0

What happens if the function $h(x) : \mathbb R \rightarrow \mathbb R^+$ does not have a minimizer at exactly 0? This would mean that, as $n$ grows, the integral should be completely exponential.

In this case, we simply factor out the value of $h$ at the minimum:
$$
\int_{-\infty}^{\infty}e^{-nh(x)}dx = \int_{-\infty}^{\infty}e^{-n[h(a)+(h(x)-h(a))]}dx = e^{-nh(a)}\int_{-\infty}^{\infty}e^{-n(h(x)-h(a))}dx \sim e^{-nh(a)}\sqrt{\frac{2\pi}{nh''(a)}}.
$$
We can combine this result with the one above including a prefactor to obtain
$$
\int_{-\infty}^{\infty}g(x)e^{-nh(x)}dx \sim g(a)e^{-nh(a)}\sqrt{\frac{2\pi}{nh''(a)}}.
$$
We see that the polynomial term with the (most) nontrivial part of the integral also has an exponentially decaying factor in the case that the minimum $h(a) > 0$.

### Laplace's Method on Minimizers with a Nonzero Derivative

Suppose we want 
$$
\int_{0}^\infty e^{-nh(x)}dx 
$$
for a function $h$ that is minimized at $h(0)$, but that $h'(0) > 0$. We would have to modify the previous Taylor's theorem approach: the first order is nonzero.

We have
$$
\int_0^\infty e^{-n h(x)}dx = \int_0^\infty e^{-n[h(0) + h'(0)x]}dx = \int_{0}^\infty = \int_0^\infty e^{-nh'(0)x} = \Theta(\frac{1}{n}).
$$
In conclusion, the Laplace method is a useful way to approximate complex integrals, and will play a critical role in our next topic of discussion: large deviations.

<!-- The previous version assumes that the minimizer is an interior point. Suppose instead that we want
$$
\int_0^\infty e^{-nh(x)}dx
$$
and suppose $h$ is minimized at $0$, but
$$
h'(0)>0.
$$
This is different because the first nonzero term in the Taylor expansion is now linear:
$$
h(x)=h(0)+h'(0)x+o(x).
$$
Therefore,
$$
\int_0^\infty e^{-nh(x)}dx \overset{(1)}\approx \int_0^\infty e^{-n(h(0)+h'(0)x)}dx \overset{(2)}= e^{-nh(0)}\int_0^\infty e^{-nh'(0)x}dx \overset{(3)}= \frac{e^{-nh(0)}}{nh'(0)}.
$$
Approximation $(1)$ is the first-order Taylor approximation near the boundary point $0$. Equality $(2)$ is because $e^{-nh(0)}$ is constant in $x$. Equality $(3)$ is the elementary integral $\int_0^\infty e^{-cx}dx=\frac{1}{c}$ with $c=nh'(0)$.

Thus,
$$
\int_0^\infty e^{-nh(x)}dx \sim \frac{e^{-nh(0)}}{nh'(0)}.
$$

<div class="callout remark"><span class="label">Remark: Boundary vs. Interior Scale</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
For an interior minimum, the first nonzero Taylor term is quadratic:
$$
h(x)-h(a)\approx \frac{1}{2}h''(a)(x-a)^2.
$$
The non-negligible region is therefore $|x-a|=O(n^{-1/2})$. For a boundary minimum with $h'(0)>0$, the first nonzero Taylor term is linear:
$$
h(x)-h(0)\approx h'(0)x.
$$
The non-negligible region is therefore $x=O(n^{-1})$. So the boundary version has order $1/n$, while the interior version has order $1/\sqrt n$.
</div> -->

<!-- ## Maximum Version of Laplace's Method

Sometimes the exponent has the opposite sign:
$$
\int_{-\infty}^{\infty}e^{nh(x)}dx.
$$
In this case, the integral is dominated by the maximizer of $h$, not the minimizer.

Suppose $t$ is the unique maximizer of $h$. Then
$$
h'(t)=0,\qquad h''(t)<0.
$$
Taylor expanding around $t$ gives
$$
h(x) \overset{(1)}= h(t)+\frac{1}{2}h''(t)(x-t)^2+O((x-t)^3) \overset{(2)}= h(t)-\frac{1}{2}|h''(t)|(x-t)^2+O((x-t)^3).
$$
Equality $(1)$ is Taylor expansion around the maximizer $t$, using $h'(t)=0$. Equality $(2)$ is because $h''(t)<0$, so $h''(t)=-|h''(t)|$.

Therefore,
$$
\int_{-\infty}^{\infty}e^{nh(x)}dx \approx e^{nh(t)}\int_{-\infty}^{\infty}e^{-\frac{n}{2}|h''(t)|(x-t)^2}dx \overset{(1)}= e^{nh(t)}\sqrt{\frac{2\pi}{n|h''(t)|}}.
$$
Equality $(1)$ is again the Gaussian integral.

Thus,
$$
\int_{-\infty}^{\infty}e^{nh(x)}dx \sim e^{nh(t)}\sqrt{\frac{2\pi}{n|h''(t)|}}.
$$
More generally,
$$
\int_{-\infty}^{\infty}g(x)e^{nh(x)}dx \sim g(t)e^{nh(t)}\sqrt{\frac{2\pi}{n|h''(t)|}}.
$$

## Logarithmic Laplace Principle

Large deviation theory usually cares about exponential rates, not exact polynomial prefactors. Therefore, it is common to take $\frac{1}{n}\log$ of the integral.

From Laplace's method,
$$
\int e^{-nh(x)}dx \sim e^{-nh(a)}\sqrt{\frac{2\pi}{nh''(a)}}.
$$
Taking logarithms and dividing by $n$,
$$
\frac{1}{n}\log\int e^{-nh(x)}dx \overset{(1)}= -h(a)+\frac{1}{2n}\log\left(\frac{2\pi}{nh''(a)}\right)+o(1) \overset{(2)}\rightarrow -h(a).
$$
Equality $(1)$ is obtained by applying $\frac{1}{n}\log$ to the Laplace asymptotic. Convergence $(2)$ is because $\frac{1}{n}\log n\rightarrow 0$.

Therefore,
$$
\frac{1}{n}\log\int e^{-nh(x)}dx \rightarrow -\inf_x h(x).
$$
For the maximum version,
$$
\frac{1}{n}\log\int e^{nh(x)}dx \rightarrow \sup_x h(x).
$$

<div class="callout remark"><span class="label">Remark: Laplace's Method Turns Integration into Optimization</span><br/>
<hr style="height:0.01px; visibility:hidden;" />
At exponential scale, integrals behave like maxima or minima. The integral averages over all possible $x$, but the exponential scaling forces the average to be dominated by the best $x$. This is the bridge to large deviation theory: rare-event probabilities are also governed by optimization problems. The probability of an event is dominated by the cheapest way for the event to occur.
</div> -->