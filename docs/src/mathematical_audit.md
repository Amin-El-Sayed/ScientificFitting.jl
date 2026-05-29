# JuFitter Mathematical Audit vs kafe2

Stand: 2026-05-07

## Kurzfazit

JuFitter hat inzwischen einen deutlich saubereren mathematischen Kern als die
erste Version. Der aktuelle Backend-Stand ist fuer 1D-XY-Fits mit Gauss-Fehlern,
vollen Kovarianzmatrizen, lokaler x-Fehler-Projektion, Bounds, Constraints,
Parameter-Priors, fixed parameters, Profilen und 2D-Contours als Fundament
brauchbar.

JuFitter ist weiterhin nicht vollstaendig identisch mit kafe2, hat aber fuer
eine v1 nun die wichtigsten wissenschaftlichen Grundfunktionen: strukturierte
Fehlerkomponenten, Poisson-, Histogram-, Unbinned-, Extended-Unbinned-, Indexed-,
Custom-, MultiFit- und korrelierte Parameterconstraint-Grundlagen,
Profile/Contours, Diagnoseplots und numerische Diagnosewarnungen.

Quellen fuer den Funktionsvergleich:

- kafe2 Documentation: <https://etpwww.etp.kit.edu/~quast/kafe2/htmldoc/index.html>
- kafe2 Mathematical Foundations: <https://etpwww.etp.kit.edu/~quast/kafe2/htmldoc/parts/mathematical_foundations.html>
- kafe2 User Guide, Profiles und Contours: <https://etpwww.etp.kit.edu/~quast/kafe2/htmldoc/parts/user_guide.html>

## Aktueller JuFitter-Status

Vorhanden:

- 1D-XY-Fits mit skalarem Modelloutput
- Diagonale und volle Kovarianzen fuer y
- Diagonale und volle Kovarianzen fuer x ueber lokale effektive Varianz
- Klassisches gewichtetes Least Squares mit

$$
\chi^2(\theta)=r(\theta)^T V^{-1}r(\theta)
$$

- Volle Gauss-NLL mit Logdeterminante

$$
-2\log L(\theta)
=
n\log(2\pi)+\log\det V(\theta)+r(\theta)^T V(\theta)^{-1}r(\theta)
$$

- Automatische Kostenfunktionswahl per `cost=:auto`
- Parameter-Priors mit symmetrischer oder asymmetrischer Unsicherheit
- Fixed parameters mit symmetrischer oder asymmetrischer reportbarer
  Unsicherheit
- Bounds und allgemeine Gleichungs-/Ungleichungsconstraints
- Optimierung freier Parameter in reduziertem Parameterraum
- Analytische Jacobians im LsqFit-Pfad
- Parameterkovarianz aus lokaler Kruemmung/Jacobian-Approximation
- Korrelationsmatrix, `chi2`, `chi2/ndf`, p-Wert, AIC und BIC
- `scale_covariance=:auto | :never | :always`
- Profilscans und 2D-Contours durch Re-Minimierung
- Plotting fuer Fit, Profile und Contours
- Programmatisch extrahierbarer `fit_report`
- `FitDiagnostics` mit Warnungen, Konditionszahlen und aktiven Bounds
- Poisson-, Histogram-, Unbinned- und Extended-Unbinned-Likelihoods
- IndexedFit-artige Datenmodelle
- CustomFit-API fuer frei definierte Kostenfunktionen
- einfache MultiFits mit Parameter-Mapping und benannten Parametern
- strukturierte Fehlerkomponenten fuer absolute, relative, modellrelative und
  Kovarianz-Fehlerquellen
- Pull-, Residual- und Ratio-Diagnoseplots
- numerische Diagnosewarnungen und robustere inverse/Hessian-Fallbacks

Noch offen fuer spaetere Versionen:

- Komfort-API fuer benannte Fit-Komponenten auf kafe2-Niveau
- automatische Report-Warnungen pro Fehlerkomponente
- erweiterte Visualisierungen wie kombinierte Fit-/Pull-Layouts

## Mathematische Bewertung

### Kostenfunktionen

Die Trennung zwischen Kostenfunktion und Optimizer ist jetzt korrekt angelegt.
Fuer statische Gauss-Kovarianzen ist die Minimierung von `:chi2` mathematisch
aequivalent zur Maximierung der Gauss-Likelihood, da konstante
Normalisierungsterme die Lage des Minimums nicht aendern.

Sobald die effektive Kovarianz parameterabhaengig ist, zum Beispiel durch
x-Unsicherheiten,

$$
V_{\mathrm{eff}}(\theta)=V_y + J_x(\theta)V_xJ_x(\theta)^T,
$$

muss der Logdeterminanten-Term in der Likelihood beruecksichtigt werden.
`cost=:auto` schaltet dafuer auf `:gaussian_nll`. Das beseitigt die wichtigste
mathematische Luecke der ersten Implementierung.

### Parameterunsicherheiten

Die Standardfehler aus `param_stderr` bleiben lokale, parabolische Fehler. Das
ist fuer gut konditionierte, lokal lineare Probleme korrekt, kann aber bei stark
nichtlinearen Problemen zu optimistisch oder symmetrisch verzerrt sein.

Profile und 2D-Contours sind deshalb wichtig. JuFitter re-minimiert dabei die
Kostenfunktion, waehrend ein oder zwei Parameter fixiert werden. Fuer eine
NLL-Konvention mit

$$
C(\theta)=-2\log L(\theta)
$$

entsprechen typische Konturlevel zum Beispiel

$$
\Delta C=1.0
$$

fuer ein ungefaehres 1-Sigma-Intervall in einem Parameter und

$$
\Delta C=2.30
$$

fuer die gemeinsame 1-Sigma-Region zweier Parameter.

### Fixed Parameters

Fixed parameters werden jetzt nicht mehr als kuenstlicher Datenpunkt oder Prior
missbraucht. Sie werden aus dem Optimierungsvektor entfernt und erst beim
Auswerten des Modells wieder in den vollen Parametervektor eingebettet. Dadurch
bleiben Optimierung, Kovarianzmatrix und NDF-Logik konsistent.

Eine angegebene Unsicherheit eines fixed parameter ist keine Constraint im Fit,
sondern eine reportbare lokale Zusatzunsicherheit fuer diesen Parameter. Das ist
die richtige Semantik, wenn ein Parameter wirklich festgehalten werden soll,
aber sein externer Fehler im Ergebnis nicht verloren gehen darf.

### NDF und Goodness of Fit

JuFitter verwendet aktuell

$$
\mathrm{ndf}=n_{\mathrm{obs}}-n_{\mathrm{free}},
$$

wobei unkorrelierte Parameter-Priors als zusaetzliche Beobachtungsterme zaehlen
und fixed parameters nicht zu den freien Parametern gehoeren. Allgemeine
Constraints und aktive Bounds sind statistisch subtiler; dort sind p-Werte und
`chi2/ndf` nur mit Vorsicht zu interpretieren.

### AIC und BIC

AIC und BIC werden aus der vollen NLL berechnet:

$$
\mathrm{AIC}=\mathrm{NLL}_{\min}+2k,
$$

$$
\mathrm{BIC}=\mathrm{NLL}_{\min}+k\log n.
$$

Dabei ist

$$
k=n_{\mathrm{free}}
$$

und

$$
n=n_{\mathrm{obs}}.
$$

AIC/BIC sind nur sinnvoll vergleichbar, wenn die verglichenen Modelle zur selben
Datenbasis und kompatiblen Likelihood-Konvention gehoeren.

## Backend-Bewertung

### LsqFit-Pfad

Geeignet fuer:

- unbeschraenkte Least-Squares-Probleme
- statische y-Unsicherheiten
- schnelle Standardfits

Nicht geeignet fuer:

- NLL mit parameterabhaengiger Kovarianz
- Bounds, Constraints oder Parameter-Priors
- Histogram-, Poisson- oder Unbinned-Likelihood-Fits im `LsqFit`-Pfad

### Optimization.jl-Pfad

Geeignet fuer:

- skalare Kostenfunktionen
- Bounds
- Constraints
- Parameter-Priors
- x-Unsicherheiten mit voller Gauss-NLL
- Poisson-, Histogram-, Unbinned-, Multi- und spaeter Custom-Likelihoods

Die AD-Konfiguration fuer constrained optimization nutzt eine explizite
Second-Order-Konfiguration. Damit ist die fruehere Warnung zur Hessian/Gradient
Kombination beseitigt.

## Funktionsvergleich gegen kafe2

| Bereich | kafe2 | JuFitter aktuell | Status |
| --- | --- | --- | --- |
| XY Gauss-Fit | ja | ja | vorhanden |
| Volle y-Kovarianz | ja | ja | vorhanden |
| x-Unsicherheiten | ja | lokale effektive Varianz + NLL | vorhanden als Approximation |
| Bounds | ja | ja | vorhanden |
| Fixed parameters | ja | ja | vorhanden |
| Gaussian constraints | ja | unkorreliert, asymmetrisch | teilweise |
| Korrelierte Parameterconstraints | ja | ja | vorhanden |
| Profile-Likelihood | ja | ja | vorhanden |
| 2D-Contours | ja | ja | vorhanden |
| Poisson-NLL | ja | ja | vorhanden |
| Histogram-Fit | ja | ja | vorhanden |
| Unbinned-Fit | ja | ja | vorhanden, auch extended |
| IndexedFit | ja | ja | vorhanden |
| CustomFit | ja | ja | vorhanden |
| MultiFit | ja | ja | Parameter-Mapping vorhanden |
| Fehlerkomponenten | ja | named absolute/relative/model-relative/covariance | vorhanden |
| Pull/Ratio/Residual-Plots | ja | ja | vorhanden |
| Text-/Datenreport | ja | ja | vorhanden, ausbaubar |

## Naechste robuste Ausbauschritte

1. `FitData`-Abstraktion einfuehren, damit XY-, Histogram-, Unbinned- und
   spaeter ODE/PDE-Fits dieselbe Cost/Optimizer-Schicht nutzen.
2. Fehlerquellen als eigene Objekte modellieren: absolut/relativ,
   datenrelativ/modellrelativ, korreliert/unkorreliert, aktivierbar.
3. Komfort-API und Dokumentation fuer benannte Workflows ausbauen.
4. Kombinierte Fit-/Pull-Layouts ergaenzen.
5. Numerische Spezialfaelle mit weiteren Referenzproblemen benchmarken.
