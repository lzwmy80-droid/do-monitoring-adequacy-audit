# Supplementary Theory Proofs

## S1. Scope and notation

This note proves Theorem 1 in the main text for the exact-inlier replacement
model. Let \(M\ge2\), \(k\in\mathbb Z_{\ge0}\), \(L\ge0\), and \(H>0\),
and define

\[
t_0<\cdots<t_{M-1},\qquad
N=M-1,\qquad
\delta_i=t_{i+1}-t_i,\qquad
h=\max_{0\le i<N}\delta_i,
\]

and let \(T=t_{M-1}-t_0\). Write
\([M]=\{0,\ldots,M-1\}\). For a report vector
\(y=(y_0,\ldots,y_{M-1})\in\mathbb R^M\), define the unrestricted-state
feasible class

\[
\mathcal F_{L,k}(y)=
\bigcup_{\substack{E\subseteq[M]\\ |E|\le k}}
\left\{
x:[t_0,t_{M-1}]\to\mathbb R:
\operatorname{Lip}(x)\le L,\ 
x(t_i)=y_i\ \text{for every }i\notin E
\right\}.
\tag{S1}
\]

The set \(E\) is a replacement witness: reports indexed by \(E\) impose no
constraint on that candidate path. Calling an index “retained” or “missing”
below refers only to this mathematical witness; it is not a claim that a
particular sensor report is known to be correct or erroneous.

Let

\[
D_H(x)=\int_{t_0}^{t_{M-1}}(H-x(t))_+\,dt,
\qquad
\mathcal R_D(y)=
\{D_H(x):x\in\mathcal F_{L,k}(y)\}.
\tag{S2}
\]

For a nonempty functional image, we use

\[
\operatorname{diam}\mathcal R_D
:=
\sup_{u,v\in\mathcal R_D}|u-v|
=
\sup\mathcal R_D-\inf\mathcal R_D.
\]

We suppress the argument \(y\) below. All retained reports in (S1) are
matched exactly. Theorem 1, the endpoint algorithm, its verification, and
the empirical analysis in the main text all use this exact-inlier model.
Nothing below asserts an endpoint algorithm or an empirical result for a
nonzero inlier tolerance.

## S2. Common retained anchors and the pairwise geometry bound

Take any \(x,z\in\mathcal F_{L,k}(y)\), and choose replacement witnesses
\(E_x,E_z\) satisfying \(|E_x|,|E_z|\le k\). Their **common retained-anchor
set** is

\[
C=[M]\setminus(E_x\cup E_z).
\tag{S3}
\]

Thus \(i\in C\) means that both candidate explanations retain report \(i\),
and consequently

\[
x(t_i)=z(t_i)=y_i.
\tag{S4}
\]

When \(C\ne\varnothing\), list its elements as
\(c_0<\cdots<c_s\), and define

\[
a=t_{c_0}-t_0,\qquad
b=t_{M-1}-t_{c_s},\qquad
g_j=t_{c_j}-t_{c_{j-1}}\quad(1\le j\le s).
\tag{S5}
\]

Here \(a\) and \(b\) are the two one-sided boundary gaps, and the \(g_j\)
are the gaps between adjacent common anchors.

**Lemma S1 (pairwise common-anchor bound).** If \(C\ne\varnothing\), then

\[
|D_H(x)-D_H(z)|
\le
L\left(
a^2+b^2+\frac12\sum_{j=1}^{s}g_j^2
\right).
\tag{S6}
\]

**Proof.** The hinge map \(u\mapsto(H-u)_+\) is 1-Lipschitz, so

\[
|D_H(x)-D_H(z)|
\le
\int_{t_0}^{t_{M-1}}|x(t)-z(t)|\,dt.
\tag{S7}
\]

Set \(q=x-z\). Because both paths are \(L\)-Lipschitz, \(q\) is
\(2L\)-Lipschitz. Equation (S4) gives \(q(t_i)=0\) at every common anchor.
Between two adjacent common anchors at \(u<v\), with \(g=v-u\),

\[
|q(t)|
\le
2L\min\{t-u,v-t\}.
\tag{S8}
\]

Integration over that interval gives

\[
\int_u^v |q(t)|\,dt\le \frac{Lg^2}{2}.
\tag{S9}
\]

On a left boundary segment of length \(a\), only its right endpoint is a
common anchor, so

\[
|q(t)|\le2L(t_{c_0}-t)
\quad\Longrightarrow\quad
\int_{t_0}^{t_{c_0}}|q(t)|\,dt\le La^2.
\tag{S10}
\]

The right boundary contributes at most \(Lb^2\) by the same argument.
Summing the disjoint segments proves (S6). \(\square\)

## S3. Missing runs

The **missing-anchor set** for the pair \((x,z)\) is

\[
U=E_x\cup E_z=[M]\setminus C,
\qquad
m=|U|.
\tag{S11}
\]

Because each witness has at most \(k\) indices,

\[
m=|E_x\cup E_z|
\le |E_x|+|E_z|
\le 2k.
\tag{S12}
\]

A **missing run** is a maximal consecutive block of indices in \(U\). A run
is a boundary run if it precedes the first common anchor or follows the last
common anchor; otherwise it is an internal run between two common anchors.
The length \(r\) of a run is its number of missing anchors. A boundary run
of length \(r\) joins \(r\) sampling intervals into one one-sided boundary
gap. An internal run of length \(r\) joins \(r+1\) sampling intervals into
one gap between common anchors.

For comparison, if every sampling point were a common anchor, the geometric
quantity in Lemma S1 would be

\[
G_0=\frac12\sum_{i=0}^{N-1}\delta_i^2.
\tag{S13}
\]

Consider first a boundary run of length \(r\), and denote its \(r\) interval
lengths by \(d_1,\ldots,d_r\). Its contribution after the anchors are
missing is \((\sum_p d_p)^2\), whereas those intervals contributed
\(\frac12\sum_p d_p^2\) to \(G_0\). Since \(0<d_p\le h\), the incremental
cost is at most

\[
\left(\sum_{p=1}^{r}d_p\right)^2
-\frac12\sum_{p=1}^{r}d_p^2
\le
\left(r^2-\frac r2\right)h^2
\equiv b(r).
\tag{S14}
\]

The expression on the left is increasing in each \(d_p>0\), so its maximum
over \(d_p\le h\) occurs at \(d_1=\cdots=d_r=h\).

For an internal run of length \(r\), let
\(d_0,\ldots,d_r\) be the \(r+1\) interval lengths between its two bounding
common anchors. The incremental cost is

\[
\frac12\left(\sum_{p=0}^{r}d_p\right)^2
-\frac12\sum_{p=0}^{r}d_p^2
=\sum_{0\le p<q\le r}d_pd_q
\le
\frac{r(r+1)}2h^2
\equiv i(r).
\tag{S15}
\]

Set \(b(0)=i(0)=0\).

**Lemma S2 (run merging).** Suppose \(k\ge1\), \(C\ne\varnothing\), and the
total number of missing anchors is \(m\le2k\). The total incremental
geometric cost of all missing runs relative to \(G_0\) is at most

\[
b(2k)=(4k^2-k)h^2.
\tag{S16}
\]

**Proof.** Let \(R\) be the total number of missing anchors in the at most
two boundary runs, and let \(S\) be the total number in all internal runs.
Thus \(m=R+S\). For nonnegative integers \(r,s\),

\[
b(r+s)-b(r)-b(s)=2rs\,h^2\ge0,
\tag{S17}
\]

and

\[
i(r+s)-i(r)-i(s)=rs\,h^2\ge0.
\tag{S18}
\]

Therefore, all boundary-run costs sum to at most \(b(R)\), and all
internal-run costs sum to at most \(i(S)\).

If \(m\ge2\) and \(R=0\), then

\[
b(m)-i(m)=\frac{m(m-2)}2h^2\ge0.
\tag{S19}
\]

If \(R\ge1\) and \(S\ge1\), then

\[
b(R+S)-b(R)-i(S)
=
S\left(2R+\frac S2-1\right)h^2
\ge0.
\tag{S20}
\]

If \(S=0\), (S17) already gives a total cost no larger than
\(b(R)=b(m)\). Hence every configuration with \(m\ge2\) has total
incremental cost at most \(b(m)\): concentrating all missing anchors into
one boundary run is the worst of these configurations under the
interval-wise bound \(h\).

The case \(m=1\) needs separate treatment. A single boundary missing anchor
costs at most \(b(1)=h^2/2\), whereas a single internal missing anchor costs
at most \(i(1)=h^2\). Since \(k\ge1\),

\[
i(1)=h^2\le b(2k).
\tag{S21}
\]

Finally, (S12) gives \(m\le2k\), and \(b(r)\) is increasing for integer
\(r\ge1\). Thus configurations with \(m\ge2\) satisfy
\(b(m)\le b(2k)\), the case \(m=1\) is covered by (S21), and \(m=0\) has
zero incremental cost. This proves (S16). \(\square\)

## S4. Proof of Theorem 1

Assume that \(\mathcal F_{L,k}(y)\) is nonempty and \(M>2k\). For any two
candidate paths \(x,z\), (S12) implies

\[
|C|=M-|E_x\cup E_z|
\ge M-2k>0.
\tag{S22}
\]

Thus the common-anchor set is nonempty and Lemma S1 always applies.

If \(k=0\), then \(E_x=E_z=\varnothing\). Every sampling point is a common
anchor, so \(a=b=0\), the adjacent common-anchor gaps are precisely the
\(\delta_i\), and Lemma S1 gives

\[
|D_H(x)-D_H(z)|
\le
\frac L2\sum_{i=0}^{M-2}\delta_i^2.
\tag{S23}
\]

Taking the supremum over all pairs \(x,z\in\mathcal F_{L,0}(y)\) proves

\[
\operatorname{diam}\mathcal R_D
\le
\frac L2\sum_{i=0}^{M-2}\delta_i^2.
\tag{S24}
\]

Now let \(k\ge1\). Equations (S13)--(S16) show that the geometric quantity
in Lemma S1 is bounded by

\[
a^2+b^2+\frac12\sum_{j=1}^{s}g_j^2
\le
\frac12\sum_{i=0}^{M-2}\delta_i^2
+(4k^2-k)h^2.
\tag{S25}
\]

Consequently, for every feasible pair,

\[
|D_H(x)-D_H(z)|
\le
L\left[
\frac12\sum_{i=0}^{M-2}\delta_i^2
+(4k^2-k)h^2
\right].
\tag{S26}
\]

The diameter of a nonempty real-valued functional image equals the
supremum of its pairwise absolute differences. Taking that supremum in
(S26) proves

\[
\operatorname{diam}\mathcal R_D
\le
L\left[
\frac12\sum_{i=0}^{M-2}\delta_i^2
+(4k^2-k)h^2
\right],
\tag{S27}
\]

which is the \(k\ge1\) branch of Theorem 1. This completes the proof.
\(\square\)

For completeness,

\[
\sum_{i=0}^{M-2}\delta_i^2
\le
h\sum_{i=0}^{M-2}\delta_i
=hT.
\tag{S28}
\]

On a regular grid with \(N=M-1\) intervals and \(h=T/N\), (S27) becomes

\[
\operatorname{diam}\mathcal R_D
\le
\frac{LT^2}{2N}
+L(4k^2-k)\frac{T^2}{N^2}.
\tag{S29}
\]

## S5. Equality construction on a regular grid

This section establishes sharpness only for the report-uniform regular-grid
design bound in the model without an external state bound. It does not say
that every fixed report vector attains (S29).

Let \(t_i=ih\), \(0\le i\le N\), and first suppose \(k\ge1\) and
\(N\ge2k\), equivalently \(M=N+1>2k\). Choose a constant \(c\) low enough
that

\[
c+2kLh<H.
\tag{S30}
\]

This choice is always possible in the unrestricted-state model. Define the
reports by

\[
y_i=
\begin{cases}
c+L(2k-i)h, & 0\le i<k,\\
c-L(2k-i)h, & k\le i<2k,\\
c, & 2k\le i\le N.
\end{cases}
\tag{S31}
\]

Define \(x\) and \(z\) on the initial segment \(0\le t\le2kh\) by

\[
x(t)=c+L(2kh-t),
\qquad
z(t)=c-L(2kh-t).
\tag{S32}
\]

For each remaining grid interval \([t_i,t_{i+1}]\), \(i\ge2k\), let

\[
\tau_i(t)=\min\{t-t_i,t_{i+1}-t\},
\tag{S33}
\]

and set

\[
x(t)=c+L\tau_i(t),
\qquad
z(t)=c-L\tau_i(t).
\tag{S34}
\]

Both paths are \(L\)-Lipschitz. Path \(x\) matches every report after using

\[
E_x=\{k,\ldots,2k-1\},
\tag{S35}
\]

and path \(z\) matches every report after using

\[
E_z=\{0,\ldots,k-1\}.
\tag{S36}
\]

Each witness contains exactly \(k\) indices. Condition (S30) places both
paths below \(H\) everywhere. Therefore,

\[
D_H(z)-D_H(x)
=\int_0^{Nh}[x(t)-z(t)]\,dt.
\tag{S37}
\]

On the initial boundary segment,

\[
\int_0^{2kh}[x(t)-z(t)]\,dt
=
\int_0^{2kh}2L(2kh-t)\,dt
=4Lk^2h^2.
\tag{S38}
\]

Each of the \(N-2k\) remaining intervals contributes

\[
\int_{t_i}^{t_{i+1}}2L\tau_i(t)\,dt
=\frac{Lh^2}{2}.
\tag{S39}
\]

It follows that

\[
D_H(z)-D_H(x)
=
Lh^2\left[
4k^2+\frac{N-2k}{2}
\right]
=
Lh^2\left(\frac N2+4k^2-k\right).
\tag{S40}
\]

The value in (S40) equals the right-hand side of (S29). The theorem supplies
the matching upper bound, so the constructed report vector has diameter
equal to the regular-grid bound.

The \(k=0\) coefficient is also attainable. Take \(y_i=c\) at every grid
point and choose \(c+Lh/2<H\). On every interval, let \(x\) be the upper
\(L\)-tent \(c+L\tau_i\) and \(z\) the lower \(L\)-tent
\(c-L\tau_i\). Both paths retain every report, and

\[
D_H(z)-D_H(x)
=N\frac{Lh^2}{2}
=\frac{LT^2}{2N}.
\tag{S41}
\]

The ability to choose \(c\) without an externally imposed lower state bound
is part of these constructions. If a state floor is imposed, a construction
may cease to be feasible. No general floor-constrained sharpness claim is
made here.

## S6. Why \(M>2k\) is needed for a report-uniform unrestricted-state bound

The following examples concern a nondegenerate monitoring horizon \(T>0\)
and the unrestricted-state model.

### S6.1. The case \(k<M\le2k\)

Partition \([M]\) into two nonempty sets \(I_-\) and \(I_+\) with

\[
|I_+|=k,
\qquad
|I_-|=M-k\le k.
\tag{S42}
\]

For any \(A>0\), define

\[
y_i=
\begin{cases}
H-A, & i\in I_-,\\
H+A, & i\in I_+.
\end{cases}
\tag{S43}
\]

The constant path \(x_-(t)=H-A\) is feasible with witness \(E_-=I_+\),
and the constant path \(x_+(t)=H+A\) is feasible with witness \(E_+=I_-\).
Both are \(L\)-Lipschitz for every prescribed \(L\ge0\). Their deficit
values are

\[
D_H(x_-)=AT,
\qquad
D_H(x_+)=0.
\tag{S44}
\]

Hence

\[
\operatorname{diam}\mathcal R_D\ge AT,
\tag{S45}
\]

which can be made arbitrarily large by changing the report amplitude \(A\)
while holding the design, \(L\), \(H\), and \(k\) fixed. Thus no finite
bound independent of report magnitude exists in this regime.

Because \(k<M\), every feasible explanation retains at least one report and
every nonempty fixed-record feasible class has finite endpoints. Indeed, a
retained report \(i\) gives
\(x(t)\ge y_i-LT\ge\min_j y_j-LT\) throughout the monitoring horizon, which
uniformly bounds \(D_H(x)\) for that fixed finite report vector. If \(L=0\)
in the construction above, every feasible path is constant and must equal
one of the two retained report levels. For that particular fixed record,

\[
\mathcal R_D=\{0,AT\}.
\tag{S46}
\]

Thus its endpoints are finite even though report-uniform control fails, and
its functional image is disconnected. The endpoint hull should not be
confused with the functional image itself.

### S6.2. The case \(k\ge M\)

If \(k\ge M\), the witness \(E=[M]\) removes every report constraint. For
each \(A>0\), the constant path

\[
x_A(t)=H-A
\tag{S47}
\]

is feasible and satisfies

\[
D_H(x_A)=AT.
\tag{S48}
\]

Therefore, even for one fixed report vector,

\[
\sup\mathcal R_D=+\infty.
\tag{S49}
\]

This is stronger than the preceding failure: for \(k<M\le2k\), each
feasible explanation retains at least one anchor and every nonempty
fixed-record feasible class has finite endpoints, but those endpoints cannot
be bounded uniformly over report magnitudes; for \(k\ge M\), the no-anchor
branch is unbounded for the fixed record itself.

## S7. Optional state-floor cap

Suppose an external state floor \(x(t)\ge B\) is added and the resulting
feasible class is nonempty. Pointwise,

\[
0\le(H-x(t))_+\le(H-B)_+.
\tag{S50}
\]

Consequently,

\[
0\le D_H(x)\le(H-B)_+T,
\qquad
\operatorname{diam}\mathcal R_D\le(H-B)_+T.
\tag{S51}
\]

In particular, if \(B<H\), the cap is \((H-B)T\); if \(B\ge H\), every
deficit is zero. In the stability argument, the main text uses the floor
only through this independent cap. When \(M>2k\), the cap and Theorem 1 may
be combined by taking the smaller upper bound. When \(M\le2k\), Theorem 1
is unavailable, but (S51) still holds. The unrestricted-state sharpness
construction and counterexamples above are not asserted to remain valid
after a floor is imposed.

## S8. Strict occupation without a threshold margin

The cumulative-deficit result does not transfer to strict threshold
occupation without an additional margin assumption. Define

\[
O_H(x)=
\int_{t_0}^{t_{M-1}}\mathbf 1\{x(t)<H\}\,dt.
\tag{S52}
\]

**Proposition S1 (no-margin obstruction).** For any finite grid
\(t_0<\cdots<t_{M-1}\) and any \(L>0\), take \(k=0\) and
\(y_i=H\) for every \(i\). In the exact-inlier unrestricted-state class,

\[
\inf_{x\in\mathcal F_{L,0}(y)}O_H(x)=0,
\qquad
\sup_{x\in\mathcal F_{L,0}(y)}O_H(x)=T.
\tag{S53}
\]

**Proof.** Let \(G=\{t_0,\ldots,t_{M-1}\}\) and let

\[
d(t,G)=\min_{0\le i<M}|t-t_i|.
\tag{S54}
\]

Distance to a fixed set is 1-Lipschitz. Hence the two paths

\[
u(t)=H+L\,d(t,G),
\qquad
\ell(t)=H-L\,d(t,G)
\tag{S55}
\]

are \(L\)-Lipschitz and satisfy
\(u(t_i)=\ell(t_i)=y_i=H\) at every report time. Because \(u(t)\ge H\),
its strict occupation is zero. On every open sampling interval,
\(d(t,G)>0\), so \(\ell(t)<H\); equality occurs only at the finitely many
sampling times, which have zero Lebesgue measure. Thus
\(O_H(\ell)=T\). Since every occupation value lies in \([0,T]\), both
endpoints in (S53) are sharp. \(\square\)

The proposition is an exact-inlier, \(k=0\), no-margin counterexample. It
shows that merely driving the maximum grid interval \(h\) to zero does not
give a report-uniform contraction guarantee for strict occupation. It does
not claim that occupation is unstable for every report vector, and it does
not rule out data-dependent guarantees under an explicit threshold margin.

## S9. Closed-form evaluation of local costs

Each boundary envelope is one affine segment, and each pair envelope is the
selected maximum or minimum of two affine segments. To evaluate a local cost,
partition its time interval at the cone intersection, at every intersection
with the state floor for a lower envelope, and at every intersection with
\(H\). On each resulting open subinterval \((a,b)\), the floor-adjusted
envelope has one affine expression \(r(t)=\alpha t+\beta\) and does not cross
\(H\). Its strict-occupation contribution is

\[
J_O(a,b)=
(b-a)\mathbf 1\{r((a+b)/2)<H\},
\tag{S56}
\]

and its cumulative-deficit contribution is

\[
J_D(a,b)=
\begin{cases}
\dfrac{b-a}{2}\{2H-r(a)-r(b)\},
& r((a+b)/2)<H,\\
0,&r((a+b)/2)\ge H.
\end{cases}
\tag{S57}
\]

Summing these contributions gives \(A_{e,F}\), \(C_{e,F}\), or \(R_{e,F}\)
as defined in the main text. Isolated threshold contacts have zero Lebesgue
measure. If an affine piece is identically \(H\), both its strict-occupation
and deficit contributions are zero. Coincident cones when \(L=0\) are
represented by one constant piece. These rules are algebraically identical
to the segment construction used by the reference implementation.
