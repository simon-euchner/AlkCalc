/* -------------------------------------------------------------------------- *
 * AlkCalc: Rechner für Alkalimetallatome und Erdalkalimetallionen            *
 *                                                                            *
 * Autor dieser Datei: Simon Euchner                                          *
 * -------------------------------------------------------------------------- */

(Englische Version: 'READMEeng.txt')


Kontakt.

    Bei etwaigen Zweifeln, Fragen und/oder Anregungen zögern Sie bitte nicht,
    den Autor zu kontaktieren.

        Simon Euchner
        Elektronische Post: <euchner.se@gmail.com>


Einführung.

        Gebundene Einelektronenzustände von Alkalimetallatomen und
    Erdalkalimetallionen sind von erheblicher Bedeutung in der
    Quanteninformationsverarbeitung, der Quantensimulation, der Rydbergphysik,
    der Quantenoptik und mehr. Doch abgesehen vom Sonderfall des
    Wasserstoffatoms und wasserstoffartiger Ionen lässt sich die
    zeitunabhängige Schrödingergleichung nicht in geschlossener Form lösen. Aus
    diesem Grund müssen Eigenenergien und Eigenzustände numerisch ermittelt
    werden. Die Bibliothek AlkCalc ist genau zu diesem Zweck geschaffen. Sie
    erlaubt es dem Anwender, den relevanten Teil des gebundenen Energiespektrums
    samt der zugehörigen Eigenzustände zu berechnen.

        AlkCalc selbst ist in zwei Teile aufgeteilt. Der erste Teil berechnet
    die Eigenenergien und Eigenzustände ein einziges Mal und legt sie auf der
    Festplatte ab. Der zweite Teil besteht aus den Bibliotheksfunktionen, welche
    nützliche Werkzeuge zum Einlesen und Verarbeiten der vorab berechneten
    atomphysikalischen Daten (Eigenenergien und Eigenzustände) darstellen. Mit
    diesen Funktionen lassen sich zum Beispiel radiale Übergangsmatrixelemente
    und Oszillatorstärken berechnen.

        AlkCalc ist sehr "zum Anfassen" gehalten, in dem Sinne, dass es einfach
    genug ist, um unmittelbar mit dem Quellkode selbst zu arbeiten. Da es zudem
    in der maschinennahen Prorammiersprache C geschrieben ist, und sich streng
    an den C99-Standard hält, ist es nahezu plattformunabhängig. Da höhere
    Sprachen üblicherweise leicht zu handhabende Wege bieten, um C-Code
    auszuführen (z. B. ctypes oder Cython für Python), ist es unkompliziert,
    eine Schnittstelle zur höheren Sprache der eigenen Wahl zu schreiben
    (etwa Python, Julia, MATLAB usw. --- was auch immer heutzutage Mode ist).
    AlkCalc sollte daher nicht allein als fertiges, einsatzbereites
    Softwarepaket verstanden werden, sondern weiter gefasst als Grundlage für
    die Entwicklung anderer, in höheren Sprachen geschriebener, Softwarepakete
    für die numerische Behandlung radialsymmetrischer Potentiale. Eine
    Schnittstelle für Python wird mit AlkCalc bereits mitgeliefert (siehe
    AlkCalc/pyalkcalc).

        Zur Modellierung der Atome und Ionen verwendet AlkCalc sogenannte
    parametrische Modellpotentiale, insbesondere jene, die 1994 von Marinescu u.
    a. in Verw. [Mar1994] sowie 1996 von Aymar u. a. in Verw. [Aym1996]
    eingeführt wurden. Diese parametrischen Modellpotentiale werden ausführlich
    in theory/theory.pdf beschrieben. Anders als übliche Bibliotheken zur
    Behandlung von Rydbergatomen, wie ARC [Sib2017] oder PairInteraction
    [Web2017], stützt sich AlkCalc nicht auf die Quantendefekttheorie.
    Stattdessen berechnet AlkCalc alle Eigenenergien und Eigenzustände
    selbstkonsistent aus dem parametrischen Modellpotential. Dies bringt den
    Vorteil mit sich, dass alle radialen Eigenzustände einem gemeinsamen
    Ursprung entstammen, und daher alle abgeleiteten Größen (Oszillatorstärken,
    Matrixelemente, Lebenszeiten, usw.) auf einer konsistenten Basis von wahren
    Eigenzuständen des parameterischen Modellpotentials basieren. Ein weiterer
    Vorteil ist, dass auch energetisch niedrig gelegene Zustände konsistent
    berechnet werden können, ohne den Teil der Wellenfunktion die den effektiven
    Kern durchdringt an eine Coulomb-Whittaker Wellenfunktion im äußeren Bereich
    anpassen zu müssen, wie es üblich ist in der Quantendefekttheorie.
    Schließlich verleiht der Verzicht auf die Quantendefekttheorie AlkCalc nicht
    nur Konsistenz, sondern auch die Fähigkeit, im Grunde jedes Problem mit
    radialer Symmetrie zu behandeln. Unter anderem ist diese Anpassungsfähigkeit
    der Grund dafür, dass es mit AlkCalc so einfach ist sowohl Rydbergzustände
    in neutralen Atomen, als auch in Ionen, von Haus aus zu behandeln --- man
    beachte, dass ARC und PairInteraction für neutrale Rydbergatome ausgelegt
    sind. Um die Schnittstelle sauber und einfach zu halten, sind alle Parameter
    der parametrischen Modellpotentiale in einer einzigen Textdatei gesammelt.
    Tatsächlich werden sämtliche zu AlkCalc gehörigen Daten in reinen
    Textdateien abgelegt, was den Vorteil hat, dass keine weiteren
    Softwarevoraussetzungen zum Lesen binärer Datenformate und dergleichen
    entstehen.

        Letztendlich ist AlkCalc für Forscher gedacht, denen die volle Kontrolle
    über die Daten hinter ihrer Forschung wichtig ist. Weiter gefasst ist
    AlkCalc für jeden gedacht, dem durchsichtige Software wichtig ist: nichts
    ist "versteckt", alles liegt offen, es werden keine unnötigen binären
    Datenformate verwendet (nur reine Textdateien), und die Software bleibt
    einfach und leichtgewichtig, ohne dabei die Vollständigkeit zu verlieren ---
    vollständig ist AlkCalc in dem Sinne, dass die Eigenenergien und
    Eigenzustände berechnet werden können, sprich, das Problem also gelöst ist.
    Zudem ist AlkCalc völlig quelloffen, und Anwender sind ausdrücklich dazu
    ermutigt, den Quellkode einzusehen und mit diesem zu arbeiten. Außerdem
    kommt AlkCalc ohne Abhängigkeiten von Drittsoftware aus, abgesehen von einem
    C-Compiler, einem FORTRAN-Compiler und einer C-Bibliothek. Der Grund dafür
    ist, dass alle sonstige benötigte Software direkt in AlkCalc selbst
    eingearbeitet ist. Die eingebaute Drittsoftware ist dabei jahrzehntelang
    erprobter Code von Netlib [Don1987].

        Neben dieser README gibt die Datei theory/theory.pdf eine gründliche
    Einführung in die Physik der Hamiltonoperatoren einzelner Atome und Ionen,
    welche AlkCalc diagonalisiert. Die Datei enthält außerdem eine vollständige
    Beschreibung der numerischen (B-Spline-)Methode, mit der AlkCalc das radiale
    Eigenwertproblem auf ein verallgemeinertes Matrixeigenwertproblem reduziert.
    Schließlich enthält theory/theory.pdf ein vollständiges Nachschlagewerk zu
    den Bibliotheksfunktionen von AlkCalc, das alle Angaben umfasst, die ein
    Anwender benötigt.


Aufbau dieses README.

        Diese README gliedert sich in die folgenden Abschnitte. Der Abschnitt
    "Softwareanforderungen" führt die benötigte Drittsoftware auf. Der Abschnitt
    "Installation" gibt Anweisungen zur richtigen Installation von AlkCalc. Der
    Abschnitt "Erzeugung der Daten" erklärt, wie die Eigenenergien und radialen
    Eigenzustände mit AlkCalc berechnet und auf der Festplatte abgelegt werden.
    Schließlich behandelt der Abschnitt "Wichtige weitere Hinweise" technische
    Gesichtspunkte, die vor der Verwendung von AlkCalc zu beachten sind.


Softwareanforderungen.

    -----------------------------------------------------------------------
    Software            Beispiel (erprobt)                        Verweis
    -----------------------------------------------------------------------
    -----------------------------------------------------------------------
    C-Bibliothek        glibc 2.44+r24+g16be1518495f-1            [GLC]

    C-Compiler          gcc 15.2.1+r22+gc4e96a094636-1            [GCC]

    FORTRAN-Compiler    gcc-fortran 15.2.1+r22+gc4e96a094636-1    [GCC]
    -----------------------------------------------------------------------

    - AlkCalc hält sich an den C99-Standard, Verw. [C99], mit der zusätzlichen
      Forderung, dass die Ganzzahltypen fester Breite int8_t, int32_t und
      int64_t vorhanden sind. Der C99-Standard behandelt diese Typen als
      optional (siehe Abschn. 7.18.1.1 in Verw. [C99]). Die meisten neueren
      C-Bibliotheken (z. B. die C-Bibliothek von GNU, Verw. [GLC]) legen sie
      jedoch fest, und der C99-Standard versichert, dass diese Typen die
      Zweierkomplementdarstellung ohne Füllbits verwenden (siehe
      Abschn. 7.18.1.1 in Verw. [C99]). AlkCalc nutzt diese Typen und ihre
      Eigenschaften, um Eindeutigkeit zu gewährleisten und Überlaufprüfungen bei
      der Ganzzahlarithmetik zu vereinfachen. Sind die Ganzzahltypen fester
      Breite auf Ihrem System nicht festgelegt, meldet der Compiler einen
      Fehler.


Installation.

    Der Installationsvorgang ist in acht Schritte aufgeteilt.

    1. Wechseln Sie in den Ordner GAUSSQ und führen Sie die Makefile aus. Dies
       erstellt das Programm GAUSSQ [Gol1969] (Quellkode von Netlib [Don1987]),
       welches zur Berechnung von Gaußquadraturen dient. Die Makefile ist für
       GNU make geschrieben. Steht GNU make auf Ihrem System nicht zur
       Verfügung, passen Sie die Makefile entsprechend an, oder führen Sie die
       darin enthaltenen Schritte von Hand aus.

    2. Wechseln Sie in den Ordner BSPLINES und führen Sie die Makefile aus (für
       die Makefile gelten dieselben Bedingungen wie in Schritt 1). Dies
       erstellt eine Version von de Boors Algorithmus [dBo2001], geschrieben von
       D. E. Amos [Amo1993].

    3. Wechseln Sie in den Ordner EIGLAPACK und führen Sie die Makefile aus
       (auch hier gelten dieselben Bedingungen wie in Schritt 1). Dies erstellt
       die beiden Eigenwertlöser DSBGVX und DSBEVX, welche Teil von LAPACK
       [And1999] sind (Quellkode von Netlib [Don1987]).

    4. Wechseln Sie in den Ordner AlkCalc/interface und setzen Sie die Variable
       PATH_TO_STATES, welche den Ort festlegt, an dem die Daten der radialen
       Eigenzustände abgelegt werden. Üblicherweise ist es unbedenklich, die
       radialen Eigenzustände unmittelbar im Datenverzeichnis von AlkCalc
       abzulegen.

--- Anmerkung (*)

    5. Wechseln Sie in den Ordner MVMBLAS und führen Sie die Makefile aus
       (dieselben Bedingungen wie in Schritt 1 gelten). Dies erstellt die
       Routine DSBMV, eine BLAS-Routine der Stufe 2. Diese Routine ist Teil von
       LAPACK [And1999,Don1987] und erlaubt es, Produkte von Matrizen und
       Vektoren effektiv numerisch zu berechnen.

    6. Wechseln Sie in den Ordner AlkCalc/interface und setzen Sie die
       relevanten Pfade in der Datei alkcalk.h. Die Variable PATH_TO_ALKCALC
       muss der vollständige Pfad zum Ort des Ordners AlkCalc sein. Die Variable
       PATH_TO_STATES legt fest, von wo auf der Festplatte die radialen
       Eigenzustände gelesen werden. Die radialen Eigenzustände können einen
       beachtlichen Speicherplatz einnehmen (~5 GB), was auf älterer Hardware zu
       einer Einschränkung werden kann. Die freie Wahl des Speicherorts der
       radialen Eigenzustände schafft hier Abhilfe; sie könnten zum Beispiel auf
       einer weiteren Festplatte abgelegt werden. Auf neuerer Hardware ist der
       Speicherplatz meist ausreichend (es stehen typischerweise >> 5 GB, eher
       ~1 TB, zur Verfügung).

    7. Führen Sie die Makefile im Ordner AlkCalc mit dem Argument "lib" aus.
       Dies erstellt die Bibliotheksfunktionen von AlkCalc, welche zum Arbeiten
       mit den von AlkCalc berechneten Eigenenergien und radialen Eigenzuständen
       dienen. Ein vollständiges Nachschlagewerk zu den Bibliotheksfunktionen
       befindet sich in theory/theory.pdf.

    8. Bevor die Bibliotheksfunktionen verwendet werden können, müssen die
       zugehörigen Daten (das heißt die Eigenenergien und die radialen
       Eigenzustände) berechnet werden. Dieser Vorgang wird im nächsten
       Abschnitt beschrieben.

(*) Dies kennzeichnet den Punkt, an dem die Eigenenergien und radialen
    Eigenzustände bereits berechnet werden können. Die folgenden Schritte sind
    nur nötig, wenn auch die Bibliotheksfunktionen von AlkCalc gebraucht werden
    (eine vollständige Beschreibung der Bibliotheksfunktionen findet sich in
    theory/theory.pdf).


Erzeugung der Daten.

        Dieser Abschnitt beschreibt, wie AlkCalc verwendet wird, um den
    vollständigen Hamiltonoperator eines einzelnen Atoms oder Ions zu
    diagonalisieren, der in theory/theory.pdf beschrieben ist. Dieser Teil von
    AlkCalc dient der Erzeugung der Daten (Eigenenergien und radiale
    Eigenzustände), auf welche die Bibliotheksfunktionen (siehe
    Abschn. Installation) angewiesen sind. Im Folgenden liegt das Augenmerk auf
    einer allgemeinen Atom- oder Ionenart X. Die vier nachstehenden Schritte
    beschreiben, wie die Daten für X mit Hilfe von AlkCalc erzeugt werden.

    1. Wechseln Sie in den Ordner AlkCalc/interface und öffnen Sie die Datei
       species.dat. Stellen Sie sicher, dass diese Datei die nötigen Daten für
       die Art X enthält. Beim Hinzufügen einer neuen Art achten Sie darauf,
       dass der Aufbau des Eintrags mit den vorhandenen Einträgen übereinstimmt.
       Beim Hinzufügen neuer Arten beachten Sie bitte die folgenden Regeln und
       Annahmen:

          (1) Achten Sie darauf, dass die Einträge nach aufsteigender
              Bahndrehimpulsquantenzahl anzuordnen sind, z. B. müssen die Daten
              für P-Zustände vor denen für F-Zustände stehen.

          (2) Die Daten für manche Bahndrehimpulsquantenzahlen können
              übersprungen werden, z. B. ist es zulässig, nur Daten für P-, D-
              und G-Zustände anzugeben. Dies ist nützlich, wenn etwa nur
              P-Zustände oder bestimmte Kreiszustände (sog. Circular states)
              berechnet werden sollen.

          (3) In species.dat muss für mindestens eine Bahndrehimpulsquantenzahl
              ein Eintrag vorhanden sein.

          (4) Angenommen, in settings.c wird eine Bahndrehimpulsquantenzahl l
              verlangt, die größer ist als die größte in species.dat angegebene,
              l0. Für l > l0 verwendet AlkCalc dann einfach l, jedoch mit den zu
              l0 gehörigen Daten aus species.dat. Dieses Verhalten ist in
              Übereinstimmung mit dem in den Verwweisen [Mar1994,Aym1996].

          (5) In species.dat wird für jedes l eine niedrigste Hauptquantenzahl,
              nl, angegeben. Dabei gibt es zwei Möglichkeiten: entweder folgt nl
              dem wasserstoffartigen Gesetz (das heißt nl = l + 1), oder nl ist
              anormal, in dem Sinne, dass nl > l + 1 gilt. Das richtige nl lässt
              sich unmittelbar aus der Elektronenanordnung der Atom- oder
              Ionenart X ablesen. Rubidium etwa besitzt die Elektronenanordnung
              [Kr]5s1, also für S-Zustände (l = 0) gilt nl = n0 = 5 > 0 + 1.
              Dies ist der anormale Fall und muss in species.dat ausdrücklich
              angegeben werden. Wird für eine verlangte
              Bahndrehimpulsquantenzahl kein Eintrag in species.dat gefunden,
              so nimmt AlkCalc an, dass das wasserstoffartige Gesetz gilt (d. h.
              nl = l + 1).

    2. Öffnen Sie die Datei settings.c und legen Sie die Parameter fest. Die
       Kennung (z. B. RB für Rubidium) der Art X wird in species.dat festgelegt.
       Ein besonders wichtiger Parameter, offset, muss von Hand und mit etwas
       Sorgfalt gewählt werden. Er bestimmt, welche Eigenenergie der niedrigsten
       Hauptquantenzahl nl zugeordnet wird. Dies ist nötig, weil zum Einen bei
       anormalen nl die niedrigste berechnete Eigenenergie nicht immer die
       richtige ist: Das Potential kann Eigenenergien unterhalb der wahren
       Grundzustandsenergie beherbergen. Zum Anderen kann es, falls das
       verallgemeinerte Eigenwertproblem schlecht konditioniert ist, dazu
       kommen, dass stark negative, unphysikalische Eigenwerte auftreten. Um
       diese unerwünschten Eigenwerte auszuschließen, gehen Sie wie folgt vor:
       Beginnen Sie mit offset = -nl + 1. Dies ordnet nl den niedrigsten vom
       Potential beherbergten Eigenwert zu, ob physikalisch oder nicht (prüfen
       Sie data/energies-X-...). Ist die zu nl gehörige Eigenenergie nicht die
       richtige, erhöhen Sie offset um eins (offset += 1) und berechnen Sie die
       Eigenenergien neu (siehe Punkt 3). Fahren Sie damit fort, bis alle
       unphysikalischen Eigenwerte abgeschnitten sind. In der Praxis liegen die
       unphysikalischen Eigenwerte meist weit von den physikalischen entfernt
       und sind daher leicht auszumachen.

    3. Führen Sie die Makefile im Ordner AlkCalc mit dem Argument "solve" aus.
       Dies erzeugt die Eigenenergien und radialen Eigenzustände. Die radialen
       Eigenzustände werden an dem vom Anwender angegebenen Ort abgelegt (siehe
       Abschn. Installation), und die Eigenenergien werden im Datenverzeichnis
       von AlkCalc abgelegt. Beim Ausführen der Makefile mit dem Argument
       "solve" gibt AlkCalc außerdem die Konditionszahl der Massenmatrix des
       verallgemeinerten Eigenwertproblems aus (siehe theory/theory.pdf). Dies
       ist lediglich als Hinweis gedacht, nicht als genaue mathematische
       Schranke, auf den möglichen Verlust an Genauigkeit, unter der Annahme,
       dass der Hamilton-Operator selbst gut konditioniert ist. Gleichwohl ist
       es nützlich für die Wahl der Ordnung der B-Splines, k, und der Anzahl der
       Knoten ohne Vielfachheit, N. Als Faustregel verliert man etwa
       log10(kappa) Stellen an numerischer Genauigkeit, wobei kappa die
       Konditionszahl der Massenmatrix M ist (siehe theory/theory.pdf). Für
       64-Bit Gleitkommazahlen sollte man versuchen, die Konditionszahl kappa
       unter ~1e6 zu halten, sodass etwa 10 Stellen an Genauigkeit übrig
       bleiben. Man beachte, dass all dies voraussetzt, dass H gut konditioniert
       ist. Diese Berechnung ist daher mehr eine grobe Abschätzung als eine
       genaue Bestimmung des numerischen Fehlers. Im Allgemeinen sollten alle
       Zahlen in einem vernünftigen Rahmen gewählt werden, sodass den
       Ergebnissen zu trauen ist.

       WICHTIG: Je nachdem, wie die höchste gewünschte Hauptquantenzahl gewählt
                ist, muss rmax groß genug gewählt werden. Man beachte dabei,
                dass hohe effektive Kernladungszahlen zu radialen Eigenzuständen
                mit Trägern bei kleineren Abständen führen --- das heißt, rmax
                muss für Atome in der Regel größer sein als für Ionen.

    4. Halten Sie von hier an die Parameter species, k, N, nmax und rmax in
       interface/settings.c FEST, und ändern Sie nur die
       Bahndrehimpulsquantenzahl, l, und die Gesamtdrehimpulsquantenzahl, j.
       Erzeugen Sie für jedes gewünschte Paar (l, j) die Eigenenergien und
       radialen Eigenzustände, indem Sie die Makefile mit dem Argument "solve"
       ausführen. Beachten Sie, dass es typischerweise nötig sein wird, den
       Parameter offset für jedes Paar (l, j) einzeln zu wählen.

       WICHTIG: Die Knotendaten (Knotenvektor und Schrittweiten; siehe
                theory/theory.pdf) werden für jede Art X genau EINMAL im
                Datenverzeichnis von AlkCalc abgelegt. Löschen Sie daher, nach
                dem Erproben verschiedener Einstellungen der zugehörigen
                Parameter, die Datei data/knotdata-X.dat von Hand. Sie wird dann
                für die schließlich gewählten Parameter genau EINMAL neu
                erzeugt. Dies mag etwas kompliziert klingen, lässt sich aber wie
                folgt zusammenfassen:

                - Prüfen Sie, mit welcher Genauigkeit die Eigenenergien
                  berechnet werden, indem Sie verschiedene Werte der Parameter
                  k, N und rmax ausprobieren.

                - Sind Sie einmal zufrieden, halten Sie die Parameter k, N und
                  rmax für alle Paare (l, j) fest, und löschen Sie, bevor Sie
                  den ersten finalen Satz an Eigenenergien und radialen
                  Eigenzuständen erzeugen, die Datei mit den Knotendaten in
                  AlkCalc/data. Wird der Eigenwertlöser hiernach für das nächste
                  Paar (l, j) ausgeführt, wird die Datei mit den Knotendaten
                  NICHT überschrieben --- sie wird nur dann neu erzeugt, wenn
                  AlkCalc die Datei mit den Knotendaten nicht auffinden kann,
                  das heißt, wenn sie gelöscht wurde.


Wichtige weitere Hinweise.

    - In der Voreinstellung ist der Massenfaktor C aus theory/theory.pdf im
      Quellkode in src/eigensolver.c auf eins gesetzt, das heißt, die reduzierte
      Masse wird durch die Elektronenmasse angenähert. Der Grund hierfür ist,
      dass die Parameter für das parametrische Modellpotential (siehe Verwweise
      [Mar1994,Aym1996]) mit C = 1 berechnet wurden: Die berechneten
      Grundzustandsenergien passen besser zu den idealen Ionisationsenergien,
      wenn C = 1 ist. Werden stattdessen Modellparameter verwendet, die MIT dem
      Massenfaktor C gewonnen wurden, muss der Massenfaktor auch im Quellkode
      verwendet werden. Dies läuft darauf hinaus, eine einzige Zeile in der
      Datei src/potential.c einzukommentieren (und eine andere
      auszukommentieren). Bitte beachten Sie dazu auch die im Quellkode gegebene
      Erklärung.

    - AlkCalc enthält das Wasserstoffatom (1H) und das Heliumion (4HE+), beide
      mit Russell-Saunders-Kopplung (LS) und, in der Voreinstellung, mit
      Massenfaktor C = 1. Diese beiden Arten dienen dazu, Ergebnisse sowohl für
      Atome als auch für Ionen zu prüfen. Man beachte, dass die LS-Kopplung
      schwach ist, besonders für kleine Bahndrehimpulsquantenzahlen und große
      Hauptquantenzahlen n. Um zu prüfen, ob die Einstellungen des
      Eigenwertlösers in interface/settings.c vernünftig sind, genügt es daher
      meist, die Eigenenergien mit dem analytischen Ergebnis für 1H und 4HE+ zu
      vergleichen, selbst ohne die LS-Kopplung zu berücksichtigen.

    - In der Voreinstellung verwendet AlkCalc für die Knoten der B-Splines
      Schrittweiten, das heißt Abstände zwischen aufeinanderfolgenden Knoten,
      die über das Intervall [0, rmax] linear zunehmen. Dies lässt sich durch
      Änderung der Funktion step in src/eigensolver.c anpassen. Weitere Angaben
      finden sich unmittelbar im entsprechenden Quellkode.

    - Alle Ionenmassen in der Datei interface/species.dat sind die Masse des
      vollständigen Atoms ABZÜGLICH der Masse des entfernten Elektrons bzw. der
      entfernten Elektronen.

    - In interface/settings.c wird die Gesamtdrehimpulsquantenzahl, j,
      festgelegt. Diese Quantenzahl nimmt halbzahlige Werte an. Beim
      numerischen Umgang mit halbzahligen Quantenzahlen ist dabei immer eine
      Entscheidung zu treffen. Hier wird die folgende Übereinkunft getroffen: um
      den numerischen Wert zu "säubern", wird die Ganzzahl 2 * j gleich
      floor(2.0 * j + 0.5) gesetzt. Die Funktion floor ist als floor(x) = k in
      den nichtnegativen Ganzzahlen definiert (x reelle Zahl), wobei k die
      eindeutige Ganzzahl ist, für die x im halboffenen Intervall [k, k + 1)
      liegt. In der Praxis bedeutet dies, dass AlkCalc für jedes j im Interval
      [0.25, 0.75) den Wert j = 1 / 2 verwendet, für jedes j in [1.25, 1.75) den
      Wert j = 3 / 2, und so fort.


Verweise.

    Beachte, dass auch die Verweise aus den übrigen Dateien von AlkCalc hier
    aufgeführt sind.

    [Mar1994] M. Marinescu, H. R. Sadeghpour, and A. Dalgarno, "Dispersion
              coefficients for alkali-metal dimers", Phys. Rev. A 49, 982 (1994)

    [Aym1996] M. Aymar, C. H. Greene, E. Luc-Koenig, "Multichannel Rydberg
              spectroscopy of complex atoms", Rev. Mod. Phys. 68, 1015 (1996)

    [Sib2017] N. Šibalić, J. D. Pritchard, C. S. Adams, and K. J. Weatherill,
              "ARC: An open-source library for calculating properties of alkali
              Rydberg atoms", Comput. Phys. Commun. 220, 319–-331 (2017)

    [Web2017] S. Weber, C. Tresp, H. Menke, A. Urvoy, O. Firstenberg,
              H. P. Büchler, and S. Hofferberth, "Calculation of Rydberg
              interaction potentials", J. Phys. B: At. Mol. Opt. Phys. 50 133001
              (2017)

    [Don1987] J. J. Dongarra and E. Grosse, "Distribution of Mathematical
              Software via Electronic Mail", Commun. ACM 30, 403--407 (1987)

    [GLC]     "The GNU C Library (glibc)", Free Software Foundation,
              URL: https://www.gnu.org/software/libc/

    [GCC]     "GCC, The GNU Compiler Collection", Free Software Foundation,
              URL: https://gcc.gnu.org/software/gcc/

    [C99]     "INTERNATIONAL STANDARD ISO/IEC 9899:1999(E)" (American National
              Standard Institute, New York, 1999) 2. Aufl.

    [Gol1969] G. H. Golub and J. H. Welsch, "Calculation of Gauss Quadrature
              Rules", Math. Comp. 23, 221--230 (1969)

    [dBo2001] C. de Boor, "A Practical Guide to Splines" (Springer, New York,
              2001) 1. Aufl., ISBN: 978-0-387-95366-3

    [Amo1993] D. E. Amos, "Implementation of de Boor's algorithm",
              URL: http://www.netlib.org/slatec/src/dbspvd.f

    [And1999] E. Anderson, Z. Bai, C. Bischof, S. Blackford, J. Demmel,
              J. Dongarra, J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney,
              and D. Sorensen, "LAPACK User's Guide" (Society for Industrial and
              Applied Mathematics, Philadelphia, PA, 1999) 3. Aufl.,
              ISBN: 0-89871-447-8

    Verweise die ausschließlich im Kode, den restlichen READMEs und den
    Dateien mit Daten von AlkCalc vorkommen.

    [NISTcuu] "The NIST Reference on Constants, Units, and Uncertainty",
              URL: https://physics.nist.gov/cuu/Constants/

    [NISTaw]  "Atomic Weights and Isotopic Compositions with Relative Atomic
              Masses", URL: https://www.nist.gov/pml/atomic-weights-and-isotopic
              -compositions-relative-atomic-masses

    [NISTie]  "NIST Atomic Spectra Database Ionization Energies Data",
              URL: https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
