/* -------------------------------------------------------------------------- *
 * AlkCalc: Alkali-Metall-Atom und Erdalkali-Metall-Ionen-Rechner             *
 *                                                                            *
 * Autor dieser Datei: Simon Euchner                                          *
 * -------------------------------------------------------------------------- */

(English version: 'READMEeng.txt')


Kontakt.

    Bei etwaigen Zweifeln, Fragen und/order Anregungen zögern Sie bitte nicht,
    den Autor zu kontaktieren.

        Simon Euchner
        Elektronische Post: <euchner.se@gmail.com>


Einführung.

        Physik mit gefangenen Atomen und Ionen ist von großem Interesse in allen
    möglichen Naturwissenschaften. Insbesondere Quantensomputer und
    Quantensimulatoren, basierend auf ultra kalten Atomen und Ionen, gewannen
    große Beliebtheit über die letzten Jahrzehnte, da heutzutage die notwendige
    Kontrolle über diese komplexen Quantensysteme Wirklichkeit ist.

        Jedoch ist es für das Konstruieren solcher Systeme Zugang zu
    atomphysikalischen Größen unabdingbar. Insbesondere sind Atome und Ionen von
    Interesse, welche mittels eines Einelektronenproblems behandelbar sind. Dazu
    wird angenommen, dass sich nur ein einzelnes Elektron nicht in dessen
    Grundzustandskonfiguration befindet, und die Restlichen den Atomkern zu
    einem effektiven Kern abschirmen, welcher als ein einzelnes Teilchen
    beschrieben wird. Schon seit den 1990er Jahren sind Modellpotentiale
    bekannt, siehe Refn. [6,8], um diese Einelektronenprobleme numerisch zu
    simulieren. Weiter gibt es Programme die diese Modellpotentiale in
    Kombination mit experimentell bestimmen Quantendefekten nutzen, wie z. B. in
    Ref. [10]. Typischerweise sind die existenten Softwarepakete recht
    umfangreich, was dazu führt, dass es kompliziert und anstrengend wird ein
    weitreichendes Verständnis für deren Funktionsweise zu erlangen. Außerdem
    sind solche Bibliotheken oftmals für eine bestimmte, stark von der Hardware
    abstrahierte, Programmiersprache entwickelt, wie z. B. Python, was
    Wissenschaftler, zumindest Projektbezogen, dazu zwingt diese
    Programmiersprache selbst zu verwenden.

        Zur Lösung der oben genannten Probleme präsentieren wie die Bibliothek
    AlkCalc (Alkali-Metall-Atom und Erdalkali-Metall-Ionen Rechner), welche mit
    einer anderen Philosophie entwickelt wurde: Wir wollen ein einfaches
    Computerprogram zur Verfügung stellen, welches alle Bedürfnisse deckt und
    gleichzeitig so schlank wie möglich gehalten ist. AlkCalc ist in nativem C
    und FORTRAN programmiert und besteht aus zwei Teilen: Den
    Bibliotheksfunktionen, welche das Arbeiten mit Eigenenergien,
    Eigenfunktionen, Clebsch-Gordan Koeffizienten und (nahezu) beliebigen
    Radialmatrixelementen erlauben, und einem Teil der dazu dient, die
    gewünschten Eigenenergien und radialen Eigenzustände numerisch zu berechnen.
    Die dabei erzeugten Daten werden dann von den Bibliotheksfunktionen
    benutzt um eine einfache Schnittstelle für den Benutzer zu kreieren.
        AlkCalc ist ein Programm zum "Anfassen", in dem Sinne, dass es einfach
    genug ist um direkt mit dem Quellkode interagieren zu können. Weil AlkCalc
    in einer hardwarenahen Sprache geschrieben ist, normgerecht bzgl. des C99
    Standards, ist es praktisch plattformunabhängig. Wichtig zu beachten ist,
    dass abstrahierte Sprachen, wie z. B. Python, typischerweise mit einfachen
    Verfahren ausgestattet sind um C code auszuführen (z. B. Ctypes oder Cython,
    um bei dem Beispiel "Python" zu bleiben), was es zu einem Leichten macht,
    eine Schnittstelle zu einer anderen Sprache zu programmieren, z. B. Python,
    Julia, Matlab, etc. --- was auch immer Mode sein mag. Daher sehen wir
    AlkCalc nicht ausschließlich als eine Bibliothek zur unmittelbaren
    Benutzung, sondern genauso als eine Basis für andere Bibliotheken, gezielt
    entwickelt für eine der stark abstrahierten Programmiersprachen.
        Des Weiteren erlaubt die einfach gehaltene Natur von AlkCalc das
    einfache Erweitern des Hamilton operators, z. B. durch das Hinzufügen von
    Termen welche die Wechselwirkung mit externen Feldern beschreiben; das Ganze
    ohne darauf angewiesen zu sein sich störungstheoretischer Methoden zu
    bedienen.
        Abschließend sei angemerkt, dass die Parameter für das Modellpotential
    in einer einzigen Textdatei enthalten sind, d. h. es ist leicht auf diese
    zuzugreifen und neue hinzuzufügen. Im Zusammenhang damit sei erwähnt, dass
    alle Daten in Form von Text gespeichert sind, sodass keine weiteren
    Programme notwendig sind um binäre oder ähnliche Datenformate
    aufzuschlüsseln.

        Insgesamt ist AlkCalc eine Software für all die Wissenschaftler die
    durchsichtige Computerprogramme sowie volle Kontrolle über diese und damit
    die Methoden, die dazu dienen Daten für Ihre Forschung zu erzeugen,
    schätzen. Es ist alles völlig offen gelegt, keine Details versteckt, keine
    unnötig komplizierten Datenformate verwendet und die gesamte Bibliothek
    schlank gehalten und dennoch vollständig; wenn die Eigenenergien und
    Eigenzustände bekannt sind, ist das gesamte Problem gelöst! Alle
    Bestandteile sind vollständig Quelloffen und der Benutzer ist dazu
    angehalten den Quellkode zu lesen und mit diesem auch zu interagieren.


Struktur dieser README Datei.

        Diese README Datei ist auf die folgende Weise Strukturiert: Im Abschnitt
    "Software Voraussetzungen" ist aufgelistet welche (wenigen) Programme
    dritter auf der Maschine vorhanden sein müssen. Das Meiste an Software
    dritter ist schon in Form von deren Quellkode (natürlich nur der tatsächlich
    benutzte) zur Verfügung gestellt. Im Abschnitt "Installation" beschreiben
    wir den Installationsprozess von AlkCalc. Im Abschnitt "Datenerzeugung"
    beschreiben wir wie die Eigenenergien und radialen Eigenzustände mit Hilfe
    von AlkCalc numerisch berechnet werden können. Abschließend sind im
    Abschnitt "Wichtige zusätzliche Informationen" (technische) Aspekte
    diskutiert, welche unbedingt vor der Benutzung von AlkCalc zu beachten sind.
        Zusätzlich zu dieser README Datei ist in "theory.d/theory.pdf" eine
    Einführung zu den Hamilton Operatoren die AlkCalc programmiert ist zu
    diagonalisieren dargelegt. Des Weiteren ist die Methode der Finiten
    Elemente, welche dazu dient, das Eigenproblem in ein Matrixproblem zu
    transformieren, beschrieben. Abschließend enthält "theory.d/theory.pdf" die
    volle Dokumentation der Bibliotheksfunktionen und alle Informationen die ein
    Benutzer bedarf.


Software Voraussetzungen.

    ----------------------------------------------------------------
    Software             Beispiel (getestet)
    ----------------------------------------------------------------
    ----------------------------------------------------------------
    BLAS                 blas-3.12.1-2 [1]

    LAPACK               lapack-3.12.1-2 [2]

    C Compiler           gcc-15.2.1+r22+gc4e96a094636-1 [12]

    FORTRAN Compiler     gcc-fortran-15.2.1+r22+gc4e96a094636-1 [12]

    C-Bibliothek         glibc-2.42+r17+gd7274d718e6f-1 [13]
    ----------------------------------------------------------------

    - AlkCalc folgt dem C99 standard, Ref. [11], mit der zusätzlichen Annahme,
      dass die Ganzzahltypen,  int8_t, int32_t und int64_t vorhanden sind. Der
      C99 Standard definiert diese Typen als optional (siehe Absch. 7.18.1.1 in
      Ref. [11]). Die meisten modernen C-Bibliotheken (z. B. die von GNU)
      implementieren diese Typen. Diese haben den Vorteil, dass der C99 Standard
      sicherstellt, dass falls diese Typen definiert sind, diese mit dem
      Zweierkomplementsystem arbeiten und keine zusätzlichen bits am Ende der
      Zahlen existieren (siehe Absch. 7.18.1.1 in Ref. [11]). Wir benutzen diese
      Eigenschaften um das Detektieren von arithmetischem Überlauf zu
      erleichtern und grundsätzlich fest definierte Zahlenbereiche zu erhalten.
      Es sei beachtet, dass falls diese Typen nicht definiert sind (und die C
      Bibliothek dem C99 standard folgt) der Compiler eine Fehlermeldung
      ausgibt.


Installation.

    Der Installationsprozess ist in sieben Schritte unterteilt.

    1. Stellen Sie sicher, dass auf der Maschine eine BLAS (z. B. Ref. [1]) und
       eine LAPACK (z. B. Ref. [2]) Installation vorhanden ist.

    2. Navigieren Sie in den Ordner LANCZOS.D und führen Sie die Makefile mit
       dem Befehl "make" aus. Das baut den relevanten Anteil von ARPACK (siehe
       Ref. [3]). Die Makefile ist für GNU make geschrieben. Gegebenenfalls
       müssen Sie diese Datei also für Ihr Makesystem anpassen, oder die Befehle
       manuell ausführen.

    3. Navigieren Sie in den Ordner LUFac.d und führen Sie auch hier die
       Makefile aus. Das baut den relevanten Anteil von SuperLU (siehe
       Ref. [7]).

    4. Navigieren Sie in den Ordner interface.d und setzen Sie die Variable
       PATH_TO_STATES, welche den Ort definiert, an dem die Daten für die
       radialen Eigenzustände später abgelegt werden.

--- Notiz (*)

    5. Im Ordner interface.d, setzen Sie die relevanten Dateipfade in der Datei
       alkcalc.h. Dabei muss PATH_TO_ALKCALC der absolute Dateipfad sein, der
       zum Ordner führt in dem AlkCalc abgelegt ist. Mit der Variable
       PATH_TO_STATES wird AlkCalc der Ort mitgeteilt, an dem die radialen
       Eigenzustände gespeichert sind. Diese nehmen einiges an Speicherplatz in
       Anspruch (~ 10 GB). Durch die Möglichkeit den Dateipfad zu den radialen
       Eigenzuständen als Benutzer selbst zu bestimmen, erhält man Flexibilität;
       z. B. kann somit für das Speichern der radialen Eigenzuständen eine
       externe Festplatte verwendet werden.

    6. Führen Sie die Makefile im Hauptordner AlkCalc mit dem Befehl "lib" aus.
       Das baut die Bibliotheksfunktionen, welche zur Endnutzerinteraktion
       bestimmt sind. Für eine ausführliche Dokumentation dieser, konsultieren
       Sie bitte "theory.d/theory.pdf".

    7. Um die Bibliotheksfunktionen benutzen zu können müssen die entsprechenden
       Daten, d. h. die Eigenenergien und die radialen Eigenzustände, zunächst
       erzeugt werden. Diesen Prozess beschreiben wir im nächsten Abschnitt.

(*) An dieser Stelle ist der Installationsprozess soweit fortgeschritten, dass
    die Eigenenergien und die radialen Eigenzustände numerisch berechnet werden
    können. Die folgenden Schritte werden nur benötigt, falls die
    Bibliotheksfunktionen (siehe "theory.d/theory.pdf" für eine vollständige
    Dokumentation) gewünscht sind.


Datenerzeugung.

        In diesem Abschnitt beschreiben wir wie AlkCalc dazu benutzt werden kann
    um (zumindest teilweise) den vollen Einelektronenhamiltonoperator, definiert
    in "theory.d/theory.pdf", zu diagonalisieren. Dieser Teil ist unverzichtbar
    um die für die Bibliotheksfunktionen benötigten Daten (Eigenenergien und
    radiale Eigenzustände) zu erzeugt. Im Folgenden ist das Vorgehen in vier
    Schritten anhand des Atoms oder Ions X erklärt.

    1. Navigieren Sie in den Ordner interface.d und öffnen Sie die Datei
       species.dat. Stellen Sie sicher, dass in dieser Datei die notwendigen
       Einträge für X vorhanden sind. Falls X neu hinzugefügt werden muss,
       folgen Sie der Formatierung der bestehenden Einträge.

    2. Öffnen Sie die Datei settings.c und setzen Sie die Parameter. Der
       Identifikator für X ist definiert in species.dat. Die zwei wichtigsten
       Parameter sind "offset" und "shift". Diese sollten auf die folgende Weise
       gewählt werden: Angenommen die Grundzustandsenergie (Ionisationsenergie)
       von X ist -0.5 Hartree. Man setze dann "offset" so, dass das Spektrum um
       den Wert 1.0 lokalisiert ist, d. h. in unserem Beispiel dass "offset"
       etwa den Wert 1.3 haben sollte. Nun kann "shift" auf 1.0 gesetzt werden,
       was bedeutet dass der Algorithmus Eigenwerte um 1.0 ausgibt. Am
       sichersten ist es mit "offset" alle Eigenenergien größer als "shift" zu
       machen. In unserem Beispiel könnte das mit dem Wert 1.6 für "offset"
       erreicht werden. Bitte beachten Sie, dass die numerisch berechneten
       Eigenenergien automatisch für die Werte von "offset" und "shift"
       korrigiert werden, sodass immer automatisch die wahren (mit der negativen
       Ionisationsenergie als Grundzustandsenergie) Eigenenergien in Einheiten
       von Hartree ausgegeben werden. Bitte beachten Sie außerdem, dass die
       Parameter "offset" und "shift" eine signifikante Auswirkung auf die
       Laufzeit des Algorithmus haben. Um zu testen welche Wahl letztendlich die
       kürzeste Laufzeit liefert, wiederholen Sie Schritt 3 für verschiedene
       Wahlen.

    3. Führen Sie die Makefile im Hauptordner AlkCalc mit dem Argument "solve"
       aus. Das erzeugt die Eigenenergien und radialen Eigenzustände. Die
       radialen Eigenzustände werden am Benutzerdefinierten (siehe
       Absch. Installation) Ort abgespeichert, und die Eigenenergien in data.d.
       Überprüfen Sie ob die Grundzustandsenergie richtig berechnet wurde und
       verändern Sie ggf. die Parameter "offset" und "shift".

       WICHTIG: Abhängig davon wie hoch die maximale Hauptquantenzahl sein soll,
                muss der maximale Radius "rmax" groß genug gewählt werden. Dabei
                ist zu beachten, dass für große Kernladungen der Träger der
                radialen Eigenzustände bei kleinen Radien liegt, d. h.,
                typischerweise muss "rmax" für Atome größer gewählt werden als
                für Ionen.

    4. Von hier an lässt man die PARAMETER species, N, nmax, rmax, in
       "interface.d/settings.c" UNVERÄNDERT und verändert nur noch die
       Bahndrehimpulsquantenzahl (l) und die Gesamtdrehimpulsquantenzahl (j).
       Für jede Wahl (l, j) erzeugt man nun die Eigenenergien und radialen
       Eigenzustände durch das Ausführen der Makefile mit dem Argument "solve".
       Wichtig ist zu beachten, dass es nötig sein kann die Parameter "offset"
       und "shift" zu verändern, während die Daten für die verschiedenen Paare
       (l, j) erzeugt werden. Beispielsweise kann es passieren, dass numerisch
       eine Eigenenergie ausgegeben wird, welche kleiner als die
       Grundzustandsenergie ist: Verändern der Parameter "offset" und "shift"
       löst dieses Problem und führt auf die wahren niedrigsten (größer als die
       Grundzustandsenergie) Energien.

       WICHTIG: Die Daten für die Diskretisierung (Orte und Schrittweiten),
                befindlich in data.d für jede Art von Atom und Ion, werden genau
                EINMAL abgespeichert. Daher muss nach dem Testen von Parametern
                die zugehörige Datei mit den Daten für die Diskretisierung
                manuell gelöscht werden. Diese wird dann für die finalen
                Parameter EINMAL erzeugt. Das mag sich etwas verwirrend anhören
                aber kann wie folgt zusammengefasst werden:

                - Überprüfen Sie durch das Testen von verschiedenen Parametern
                  ob die Eigenenergien effizient berechnet werden, wie
                  beschrieben unter Punkt 3.

                - Sobald die Wahl für die Parameter getroffen ist, hält man alle
                  Parameter fest, bis auf (l, j). Bevor die ersten Daten
                  erzeugt werden, löscht man die Datei mit den Daten für die
                  Diskretisierung im Ordner "data.d". Sobald der Algorithmus zur
                  Datenerzeugung ein weiteres Mal ausgeführt wird, für das
                  nächste Paar (l, j), wird die Datei mit den Daten für die
                  Diskretisierung nicht noch einmal erzeigt --- die Datei wird
                  nur erzeugt, falls diese nicht gefunden werden konnte, d. h.
                  gelöscht wurde.


Wichtige zusätzliche Informationen.

    - Per Standardeinstellung ist der Faktor C für die Massenkorrektur,
      definiert in "theory.d/theory.pdf", im Quellkode in "src.d/eigensolver.c"
      auf eins gesetzt. Das ist der Tatsache geschuldet, dass die
      Modellparameter (siehe Refn. [6,8]) die wir standardmäßig für das
      Modellpotential verwenden ohne die Massenkorrektur ermittelt wurden. Das
      schließen wir aus der Tatsache, dass die Grundzustandsenergien besser mit
      den wahren Ionisationsenergien übereinstimmen wenn wir die Massenkorrektur
      im Quellkode ignorieren. Natürlich, wenn Modellparameter verwendet werden
      die UNTER dem Einfluss der Massenkorrektur ermittelt wurden, muss diese
      auch im Quellkode berücksichtigt werden. Um die Massenkorrektur im
      Quellkode einzubinden, genügt es eine einzelne Zeile in
      "src.d/potential.c" einzukommentieren (und dafür eine andere
      auszukommentieren). Bitte beachten Sie dazu auch die Informationen in
      "src.d/potential.c".

    - AlkCalc enthält den Atomtyp Wasserstoff (1H) und den Ionentyp Helium
      (4He+). Beide sind inklusive Russell-Saunders (LS) Kopplung beschrieben,
      sowie, per Standardeinstellung, OHNE Massenkorrektur, d. h. wir nehmen an,
      dass die reduzierte Masse der Elektronenmasse entspricht. Wir stellen
      diese Typen bereit um Daten erzeugt zu können die auch durch
      analytische Methoden zugänglich sind; LS Kopplung ist schwach
      (insbesondere für kleine Bahndrehimpulsquantenzahlen l und große
      Hauptquantenzahlen n). Mit den bekannten Formeln kann somit bestimmt
      werden ob die in "interface.d/settings.c" gesetzten Parameter sinnvoll
      sind.

    - Für die Diskretisierung des Intervalls [0, rmax] benutzt AlkCalc linear
      größer werdende Schrittweiten. Das kann jedoch leicht verändert werden
      indem man die Funktion "step" in "src.d/eigensolver.d" umdefiniert.
      Weitere Informationen dazu sind direkt dem zugehörigen Quellkode zu
      entnehmen.

    - Alle Ionenmassen in der Datei "interface.d/species.dat" sind die Massen
      des vollständigen Atoms MINUS der Masse des Elektrons.


Referenzen.

    Beachten Sie, dass die Referenzen aus dem Quellkode hier auch aufgelistet
    sind.

    [1] "OpenBLAS: An optimised BLAS Library",
        URL: http://www.openmathlib.org/OpenBLAS/
    [2] "LAPACK---Linear Algebra PACKage", URL: https://netlib.org/lapack/
    [3] "Opencollab, ARPACK-NG", URL: https://github.com/opencollab/arpack-ng
    [4] J. W. P. Wilkinson, K. Bolsmann, T. L. M. Guedes, M. M\"uller, and
        I. Lesanovsky, New J. Phys. 27, 064502 (2025)
    [5] "NIST: National Institute of Standards and Technology (NIST)",
        URL: https://www.nist.gov
    [6] M. Aymar, C. H. Greene, E. Luc-Koenig, Rev. Mod. Phys. 68, 1015 (1996)
    [7] "SuperLU 7.0.1",
        URL: https://portal.nersc.gov/project/sparse/superlu/ug.pdf
    [8] M. Marinescu, H. R. Sadeghpour, und A. Dalgarno, Phys. Rev. A 49, 982
        (1994)
    [9] "Comission on isotopic abundances and atomic weights (CIAAW)",
        URL: https://www.ciaaw.org/lithium.htm
   [10] N. Šibalić, J. D. Pritchard, C. S. Adams, und K. J. Weatherill, Comput.
        Phys. Commun. 220, 319–331 (2017)
   [11] INTERNATIONAL STANDARD ISO/IEC 9899:1999(E) (American National Standard
        Institute, New York, 1999) 2nd ed.
   [12] "GCC, The GNU Compiler Collection", URL: https://gcc.gnu.org
   [13] "The GNU C Library", URL: https://www.gnu.org/software/libc
