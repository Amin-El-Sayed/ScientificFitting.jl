# JuFitter Backend Design

Stand: 2026-05-07

Dieses Dokument beschreibt den neu eingefuehrten Backend-Kern. Formeln sind in
LaTeX gesetzt, damit sie in Markdown-Renderern sauber lesbar sind.

## Ziel

Der Solver soll nicht mehr implizit definieren, welche Statistik minimiert
wird. JuFitter trennt daher:

$$
\text{FitProblem} + \text{CostFunction} + \text{Optimizer}
\longrightarrow \text{FitResult}.
$$

Aktuell sind zwei Kostenfunktionen implementiert:

$$
\texttt{:chi2}
$$

und

$$
\texttt{:gaussian\_nll}.
$$

`cost=:auto` waehlt derzeit:

$$
\texttt{:gaussian\_nll}
$$

falls die effektive Kovarianz parameterabhaengig ist, insbesondere bei
x-Unsicherheiten. Sonst wird
$$
\texttt{:chi2}
$$

verwendet.

## Chi-Quadrat-Kostenfunktion

Fuer Residuen

$$
r(\theta)=y-m(x,\theta)
$$

und Kovarianzmatrix

$$
V
$$

ist

$$
\chi^2(\theta)=r(\theta)^T V^{-1} r(\theta).
$$

Bei diagonalen Unsicherheiten gilt:

$$
\chi^2(\theta)=
\sum_i
\left(
\frac{y_i-m(x_i,\theta)}{\sigma_i}
\right)^2.
$$

Gaussian parameter priors werden als additive Strafterme behandelt:

$$
\chi^2_{\mathrm{prior}}(\theta)
=
\sum_j
\left(
\frac{\theta_j-\mu_j}{\tau_j}
\right)^2.
$$

Die minimierte Chi-Quadrat-Kostenfunktion ist:

$$
C_{\chi^2}(\theta)
=
\chi^2_{\mathrm{data}}(\theta)
+
\chi^2_{\mathrm{prior}}(\theta).
$$

## Volle Gauß-NLL

Bei multivariat normalverteilten Daten gilt:

$$
y \sim \mathcal{N}(m(x,\theta), V(\theta)).
$$

JuFitter verwendet die Konvention

$$
\mathrm{NLL}(\theta)=-2\log L(\theta).
$$

Damit ist:

$$
\mathrm{NLL}_{\mathrm{data}}(\theta)
=
n\log(2\pi)
+
\log\det V(\theta)
+
r(\theta)^T V(\theta)^{-1} r(\theta).
$$

Der Term

$$
\log\det V(\theta)
$$

ist konstant, wenn die Kovarianz parameterunabhaengig ist. Er ist aber
notwendig, wenn die Kovarianz von den Parametern abhaengt.

Gaussian priors werden in der NLL inklusive Normierung addiert:

$$
\mathrm{NLL}_{\mathrm{prior}}(\theta)
=
\sum_j
\left[
\log(2\pi\tau_j^2)
+
\left(
\frac{\theta_j-\mu_j}{\tau_j}
\right)^2
\right].
$$

Die volle Kostenfunktion ist:

$$
C_{\mathrm{NLL}}(\theta)
=
\mathrm{NLL}_{\mathrm{data}}(\theta)
+
\mathrm{NLL}_{\mathrm{prior}}(\theta).
$$

## Effektive Varianz bei x-Unsicherheiten

Fuer kleine x-Unsicherheiten verwendet JuFitter weiterhin die lokale
effektive-Varianz-Approximation:

$$
V_{\mathrm{eff}}(\theta)
=
V_y
+
J_x(\theta) V_x J_x(\theta)^T.
$$

Punktweise diagonal wird daraus:

$$
\sigma_{\mathrm{eff},i}^2(\theta)
=
\sigma_{y,i}^2
+
\left(
\frac{\partial m(x_i,\theta)}{\partial x}
\right)^2
\sigma_{x,i}^2.
$$

Weil diese Kovarianz von

$$
\theta
$$

abhaengen kann, waehlt `cost=:auto` hier die volle Gauß-NLL.

## FitResult-Semantik

`FitResult.stats` unterscheidet jetzt:

$$
\texttt{cost}
$$

die gewaehlte Kostenfunktion,

$$
\texttt{cost\_min}
$$

den tatsaechlich minimierten Wert,

$$
\texttt{nll\_min}
$$

die volle Gauß-NLL am Minimum, soweit fuer das aktuelle XY-Modell definiert,
und

$$
\chi^2
$$

als Goodness-of-Fit-Groesse fuer die vorhandenen Gauß-Residualterme.

AIC und BIC werden aus der NLL berechnet:

$$
\mathrm{AIC}=\mathrm{NLL}_{\min}+2k
$$

und

$$
\mathrm{BIC}=\mathrm{NLL}_{\min}+k\log n.
$$

## Solver-Auswahl

`LsqFit` wird nur fuer unbeschraenkte statische Chi-Quadrat-Fits genutzt.
Dabei werden analytische Jacobians, falls angegeben, jetzt an den Solver
weitergereicht.

Allgemeinere Kostenfunktionen, Bounds, Constraints, Priors und
parameterabhaengige Kovarianzen laufen ueber `Optimization.jl`.

## Implementierter Backend-Schritt

Dieser Backend-Schritt ist das Fundament, aber noch nicht das Ende:

- `scale_covariance` ist jetzt eine explizite Policy
  `:auto | :never | :always`.
- Fixed parameters haben eine eigene Parameterraum-Abbildung. Der Optimizer
  sieht nur freie Parameter.
- Asymmetrische Unsicherheiten fixer Parameter und asymmetrische Gaussian
  parameter priors werden lokal unterstuetzt.
- Profile und 2D-Contours re-minimieren die gewaehlte Kostenfunktion mit
  fixierten Parametern.

## Noch offen

Die naechsten mathematisch wichtigen Erweiterungen sind:

- Poisson-NLL, Histogram-Fits, Unbinned- und Extended-Unbinned-Fits,
  IndexedFit, CustomFit, MultiFits mit Parameter-Mapping und korrelierte
  Gaussian parameter constraints sind jetzt implementiert.
- Benannte/deaktivierbare Fehlerkomponenten, Diagnoseplots und robuste inverse
  Fallbacks sind jetzt implementiert.
- Fuer spaetere Versionen bleiben Komfort-APIs fuer komplexe Fehlerbudgets und
  weitere numerische Spezialfall-Benchmarks sinnvoll.
