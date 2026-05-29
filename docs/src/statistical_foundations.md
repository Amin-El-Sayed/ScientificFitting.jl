# JuFitter Statistical Foundations

Stand: 2026-05-07

Dieses Dokument beschreibt die mathematische Grundlage, auf der JuFitter
modular erweitert werden soll. Least Squares ist dabei ein Spezialfall einer
Likelihood-basierten Sicht, nicht die allgemeine Basis aller Fits.

## 1. Architekturprinzip

Ein Fit soll konzeptionell in getrennte Schichten zerfallen:

$$
\mathrm{Data} + \mathrm{Model} + \mathrm{UncertaintyModel}
+ \mathrm{CostFunction} + \mathrm{Optimizer}
\longrightarrow
\mathrm{FitResult}.
$$

Der Optimizer bestimmt nicht die Statistik. Er minimiert nur die vorher sauber
definierte Kostenfunktion.

## 2. Notation

Wir verwenden:

- Datenvektor $d \in \mathbb{R}^n$
- Modellvorhersage $m(\theta) \in \mathbb{R}^n$
- Fitparameter $\theta \in \mathbb{R}^p$
- Residuum $r(\theta)=d-m(\theta)$
- Kovarianzmatrix $V(\theta) \in \mathbb{R}^{n\times n}$
- negative Log-Likelihood in der Konvention $\mathrm{NLL}=-2\log L$

Die Konvention $-2\log L$ ist praktisch, weil Likelihood-Ratio-Differenzen
unter regulaeren Bedingungen asymptotisch Chi-Quadrat-verteilte Groessen liefern
koennen.

## 3. Gauss-Likelihood

Wenn

$$
d \sim \mathcal{N}(m(\theta), V(\theta)),
$$

ist die volle negative Log-Likelihood

$$
\mathrm{NLL}(\theta)
=
n\log(2\pi)
+
\log\det V(\theta)
+
r(\theta)^T V(\theta)^{-1}r(\theta).
$$

Der quadratische Term ist

$$
\chi^2(\theta)=r(\theta)^T V(\theta)^{-1}r(\theta).
$$

Falls $V$ parameterunabhaengig ist, sind $n\log(2\pi)$ und $\log\det V$ fuer die
Lage des Minimums konstant. Dann ist die Minimierung von $\chi^2$ aequivalent
zur Maximierung der Gauss-Likelihood. Falls $V$ von $\theta$ abhaengt, muss
$\log\det V(\theta)$ in der Kostenfunktion bleiben.

## 4. Diagonale und korrelierte Fehler

Fuer unkorrelierte bekannte Fehler gilt

$$
V = \operatorname{diag}(\sigma_1^2,\ldots,\sigma_n^2)
$$

und damit

$$
\chi^2(\theta)
=
\sum_i
\left(
\frac{d_i-m_i(\theta)}{\sigma_i}
\right)^2.
$$

Fuer korrelierte Fehler bleibt die Matrixform

$$
\chi^2(\theta)=r(\theta)^T V^{-1}r(\theta)
$$

massgeblich. Numerisch soll $V^{-1}$ nicht explizit gebildet werden. Stattdessen
wird ein Faktor $F$ mit

$$
V=FF^T
$$

verwendet. Dann loest man

$$
Fz=r
$$

und berechnet

$$
\chi^2=z^Tz.
$$

## 5. Fehlerquellen

Langfristig sollte JuFitter Fehlerquellen als Komponenten modellieren:

$$
V_{\mathrm{total}}(\theta)=\sum_k V_k(\theta).
$$

Wichtige Komponenten sind:

- absolute unkorrelierte Fehler
- absolute korrelierte Fehler
- datenrelative Fehler
- modellrelative Fehler
- x-Fehler
- externe Parameterconstraints

Modellrelative Fehler machen $V(\theta)$ parameterabhaengig und erfordern daher
die volle Gauss-NLL.

## 6. x-Unsicherheiten

Fuer ein 1D-Modell $y=f(x,\theta)$ und kleine x-Unsicherheiten kann die lokale
effektive Varianz erster Ordnung verwendet werden:

$$
V_{\mathrm{eff}}(\theta)
=
V_y + J_x(\theta)V_xJ_x(\theta)^T.
$$

Dabei ist

$$
(J_x)_{ij}=\frac{\partial f(x_i,\theta)}{\partial x_j}.
$$

Fuer punktweise Modelle mit unkorrelierten x-Fehlern folgt

$$
\sigma_{\mathrm{eff},i}^2(\theta)
=
\sigma_{y,i}^2
+
\left(
\frac{\partial f(x_i,\theta)}{\partial x}
\right)^2
\sigma_{x,i}^2.
$$

Diese Approximation ist sinnvoll fuer kleine x-Fehler und glatte Modelle. Fuer
grosse x-Fehler oder starke Nichtlinearitaet ist ein echtes
Errors-in-Variables-Modell bzw. Orthogonal Distance Regression sauberer.

## 7. Parameter-Priors und fixed parameters

Ein unkorrelierter Gaussian prior fuer Parameter $\theta_j$ hat die Form

$$
C_{\mathrm{prior},j}(\theta)
=
\left(
\frac{\theta_j-\mu_j}{\tau_j}
\right)^2.
$$

In der vollen NLL kommt die Normierung dazu:

$$
\mathrm{NLL}_{\mathrm{prior},j}(\theta)
=
\log(2\pi\tau_j^2)
+
\left(
\frac{\theta_j-\mu_j}{\tau_j}
\right)^2.
$$

Asymmetrische Priors verwenden lokal

$$
\tau_j =
\begin{cases}
\tau_j^- & \theta_j < \mu_j,\\
\tau_j^+ & \theta_j \ge \mu_j.
\end{cases}
$$

Fixed parameters sind etwas anderes als enge Priors. Ein fixed parameter wird
aus dem optimierten Parametervektor entfernt. Eine angegebene Unsicherheit eines
fixed parameter beschreibt eine externe reportbare Unsicherheit, sie zieht den
Fit aber nicht zum Wert hin, weil der Wert bereits festgehalten wird.

## 8. Parameterkovarianz

Fuer Least-Squares-Fits mit lokal linearem Modell gilt mit gewichtetem Jacobian
$J_w$

$$
\operatorname{Cov}(\hat\theta)
\approx
(J_w^TJ_w)^{-1}.
$$

Falls die Messunsicherheiten unbekannt waren und aus den Residuen geschaetzt
werden, wird oft skaliert mit

$$
s^2=\frac{\chi^2}{\mathrm{ndf}}.
$$

Dann gilt

$$
\operatorname{Cov}(\hat\theta)
\approx
s^2(J_w^TJ_w)^{-1}.
$$

Bei bekannten experimentellen Unsicherheiten soll diese Skalierung nicht
automatisch erzwungen werden. JuFitter nutzt deshalb die Policy
`scale_covariance=:auto | :never | :always`.

Fuer allgemeine NLL-Fits ist die lokale Kruemmung massgeblich:

$$
\operatorname{Cov}(\hat\theta)
\approx
2H^{-1},
$$

wobei

$$
H = \nabla^2 \mathrm{NLL}(\hat\theta)
$$

in der Konvention $\mathrm{NLL}=-2\log L$ ist.

## 9. Profile und Contours

Lokale symmetrische Fehler koennen bei nichtlinearen Modellen irrefuehrend
sein. Ein Profil fuer Parameter $\theta_i$ fixiert diesen Parameter auf einen
Rasterwert $a$ und minimiert alle anderen freien Parameter neu:

$$
C_{\mathrm{prof}}(a)
=
\min_{\theta_{-i}} C(\theta_i=a,\theta_{-i}).
$$

Die relevante Groesse ist

$$
\Delta C(a)=C_{\mathrm{prof}}(a)-C_{\min}.
$$

Fuer zwei Parameter entsteht entsprechend

$$
C_{\mathrm{prof}}(a,b)
=
\min_{\theta_{-(i,j)}} C(\theta_i=a,\theta_j=b,\theta_{-(i,j)}).
$$

Bei $C=-2\log L$ sind typische Wilks-Level $\Delta C=1.0$ fuer ein
Ein-Parameter-Intervall und $\Delta C=2.30$ fuer eine gemeinsame 1-Sigma-Region
in zwei Parametern.

## 10. Poisson-Modelle

Fuer Zaehlstatistik mit beobachteten Counts $n_i$ und Erwartungswerten
$\mu_i(\theta)$ ist

$$
P(n_i\mid\mu_i)=e^{-\mu_i}\frac{\mu_i^{n_i}}{n_i!}.
$$

Bis auf datenabhaengige Konstanten ist

$$
-2\log L(\theta)
=
2\sum_i\left(\mu_i(\theta)-n_i\log\mu_i(\theta)\right).
$$

Eine Goodness-of-Fit-Statistik ist die Poisson-Deviance

$$
D(\theta)=
2\sum_i
\left[
\mu_i(\theta)-n_i+n_i\log\frac{n_i}{\mu_i(\theta)}
\right],
$$

wobei der Term $n_i\log(n_i/\mu_i)$ fuer $n_i=0$ als $0$ definiert wird.

## 11. Histogram- und Unbinned-Fits

Histogram-Fits koennen ueber erwartete Bin-Inhalte $\mu_i(\theta)$ mit
Poisson-Likelihood beschrieben werden. Fuer eine Dichte $f(x\mid\theta)$ und
Bin $[a_i,b_i]$ gilt

$$
\mu_i(\theta)=N\int_{a_i}^{b_i}f(x\mid\theta)\,dx.
$$

Unbinned Fits verwenden direkt

$$
\log L(\theta)=\sum_i \log f(x_i\mid\theta)
$$

oder als extended likelihood zusaetzlich den Poisson-Term fuer die Gesamtzahl.

## 12. NDF, p-Werte, AIC und BIC

Fuer einfache Gauss-Fits verwendet JuFitter

$$
\mathrm{ndf}=n_{\mathrm{obs}}-n_{\mathrm{free}}.
$$

Priors koennen als zusaetzliche Beobachtungsterme zaehlen; fixed parameters
zaehlen nicht zu den freien Parametern. Allgemeine Constraints, aktive Bounds
und nicht-Gauss-Likelihoods machen die Interpretation von `chi2/ndf` und
p-Werten schwieriger und muessen im Report transparent gemacht werden.

AIC und BIC werden aus der NLL berechnet:

$$
\mathrm{AIC}=\mathrm{NLL}_{\min}+2k,
$$

$$
\mathrm{BIC}=\mathrm{NLL}_{\min}+k\log n.
$$

Sie sind nur sinnvoll vergleichbar, wenn Modelle dieselbe Datenbasis und eine
kompatible Likelihood-Konvention nutzen.

## 13. Implementierungsregeln

- Keine explizite Matrixinverse in Kostenfunktionen.
- Analytische Ableitungen bevorzugen, danach Automatic Differentiation,
  Finite Differences nur als Fallback.
- Parameterkovarianz muss zur tatsaechlich minimierten Kostenfunktion passen.
- Dense Kovarianzen sind mathematisch unterstuetzt, aber fuer grosse Daten nur
  begrenzt praktikabel.
- Reports muessen Kostenfunktion, NDF-Definition, Kovarianzskalierung und
  relevante Warnungen sichtbar machen.

## 14. Ziel-Modulstruktur

Langfristig sollte der Kern auf abstrakte Schichten hinauslaufen:

- `AbstractFitData`: XY, indexed, histogram, unbinned
- `AbstractModel`: vectorized model, density model, count model, ODE/PDE model
- `AbstractUncertainty`: absolute, relative, covariance, x-error, parameter constraint
- `AbstractCostFunction`: chi-square, Gaussian NLL, Poisson NLL, unbinned NLL, custom
- `FitProblem`: data, model, uncertainty, parameters, options
- `FitResult`: estimates, covariance, statistics, diagnostics, report

## 15. Referenzen

- kafe2 Mathematical Foundations: <https://etpwww.etp.kit.edu/~quast/kafe2/htmldoc/parts/mathematical_foundations.html>
- kafe2 API Documentation: <https://etpwww.etp.kit.edu/~quast/kafe2/htmldoc/parts/api_documentation/index.html>
- Particle Data Group, Probability and Statistics reviews: <https://pdg.lbl.gov/>
- G. Cowan, *Statistical Data Analysis*, Oxford University Press, 1998.
- G. Cowan, K. Cranmer, E. Gross, O. Vitells, *Asymptotic formulae for likelihood-based tests of new physics*, Eur. Phys. J. C 71, 1554, 2011.
- F. James and M. Roos, *MINUIT*, Computer Physics Communications 10, 343-367, 1975.
- P. T. Boggs et al., ODRPACK / Orthogonal Distance Regression, NIST.
