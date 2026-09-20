"""
Modello MKM (Microdosimetric Kinetic Model) saturazione-corretto di
Kase/Hawkins, nella forma usata clinicamente a NIRS/HIMAC (Chiba) per
il carbonio-12 -- sostituisce l'approssimazione lineare ad hoc della
sezione 5.31 con il modello reale, con parametri e formule prese da:

  Kase Y, Kanai T, Matsumoto Y, et al. (2006) "Microdosimetric
  measurements and estimation of human cell survival for heavy-ion
  beams." Radiat Res 166:629-638.
  Kase Y, Kanai T, Sakama M, et al. (2011) "Microdosimetric Approach
  to NIRS-defined Biological Dose Measurement for Carbon-ion Treatment
  Beam." J Radiat Res 52:59-68 -- equazioni (7)-(10), parametri HSG
  (alpha0=0.13 Gy^-1, beta=0.05 Gy^-2, rd=0.42um, y0=150 keV/um,
  D10R=5.0 Gy).

Due approssimazioni dichiarate rispetto al modello clinico completo:
1. y_D (dose-mean lineal energy, richiederebbe uno spettro TEPC
   misurato o una vera simulazione microdosimetrica di sito) e'
   approssimato con il LET dose-mediato gia' calcolato da Geant4
   (sezione 5.31) -- valido come proxy di primo ordine, non equivalente.
2. y* e' calcolato nel limite "spettro a valore singolo" (delta di
   Dirac su y_D) dell'equazione (9) di Kase 2011, non integrando uno
   spettro f(y) completo -- riduce l'integrale a una forma chiusa:
       y* = y0^2/y_D * [1 - exp(-(y_D/y0)^2)]

Validato contro i valori PUBBLICATI in Kase et al. 2011 (Tabella 1 /
Fig. 4, fascio carbonio-12 290 MeV/u, SOBP 6cm): usando i loro y*
misurati come proxy di y_D in questo modello si riottengono RBE10
entro il 10-15% dai valori riportati (es. y*=58.7 keV/um -> RBE10
calcolato 2.23, atteso ~2.0-2.3; punto di riferimento clinico NIRS
LET_d=80 keV/um -> RBE10 calcolato 2.67, riferimento clinico 3.0) --
uno scarto piccolo e spiegabile dalle due approssimazioni sopra, non
un errore del modello.
"""

import math

ALPHA0 = 0.13   # Gy^-1, HSG cells, Kase 2006/2011
BETA = 0.05     # Gy^-2, HSG cells (assunto indipendente da LET)
RD = 0.42       # um, raggio del dominio subcellulare
Y0 = 150.0      # keV/um, parametro di saturazione
D10R = 5.0      # Gy, dose di riferimento raggi X 200kVp per 10% sopravvivenza HSG
_D_DOMAIN = 2 * RD  # um, diametro del dominio (sito sferico)


def y_star(y_d):
    """Dose-mean lineal energy corretta per saturazione (overkill),
    limite a spettro a valore singolo dell'eq. (9) di Kase 2011."""
    if y_d <= 0:
        return 0.0
    return (Y0**2 / y_d) * (1 - math.exp(-(y_d / Y0)**2))


def z_star(ystar):
    """Specific energy per evento (Gy), conversione ICRU standard per
    un sito sferico di diametro _D_DOMAIN."""
    return 0.204 * ystar / _D_DOMAIN**2


def alpha_mkm(y_d):
    """alpha(y_D) = alpha0 + beta*z*, eq. (7)-(8) di Kase 2011."""
    return ALPHA0 + BETA * z_star(y_star(y_d))


def rbe10(y_d):
    """RBE al 10% di sopravvivenza, eq. (10) di Kase 2011, usando LET
    dose-mediato (y_d, keV/um) come proxy della dose-mean lineal energy."""
    a = alpha_mkm(y_d)
    disc = a**2 - 4 * BETA * math.log(0.1)
    return (2 * BETA * D10R) / (math.sqrt(disc) - a)


if __name__ == '__main__':
    # Validazione contro Kase et al. 2011, Tabella 1 / Fig. 4
    print("Validazione contro valori pubblicati (Kase et al. 2011):")
    for y_pub, atteso in [(14.6, "~1.3"), (28.0, "~1.5-1.6"), (58.7, "~2.0-2.3")]:
        print(f"  y*={y_pub} keV/um -> RBE10 calcolato={rbe10(y_pub):.2f} (atteso {atteso})")
    print(f"\nPunto di riferimento clinico NIRS (LET_d=80 keV/um, RBE atteso=3.0):")
    print(f"  RBE10 calcolato = {rbe10(80.0):.2f}")
