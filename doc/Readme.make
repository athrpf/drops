=============
  Makefile:
=============

-----------------------------
  Hinweise fuer DROPS-User:
-----------------------------

Zunaechst muss in der Datei drops.conf im DROPS-Rootverzeichnis die verwendete
Rechnerarchitektur eingetragen werden (z.B. LINUX). Die Compilereinstellungen
erfolgen dann in der Datei arch/<Architekur>/mk.conf .

Im DROPS-Rootverzeichnis befindet sich das top-level-Makefile. Mit "make <rule>"
wird die entsprechende Regel ausgefuehrt, wobei <rule> fuer eine der folgenden
Regeln steht:

    all        generiert automatisch die Abhaengigkeiten und
               erzeugt dann alle ausfuehrbaren Programme in DROPS.
    
    clean      loescht alle Objektdateien sowie alle ausfuehrbaren Dateien.
    
    distclean  wie "clean", loescht zusaetzlich alle Dateien mit Endungen
               .off, .dat sowie geom/topo.cpp und das dependency-file.
	       
    dep        erzeugt automatisch die Abhaengigkeiten und legt ein
               entsprechendes dependency-file an.
	       
In den jeweiligen Unterverzeichnissen befinden sich die lokalen Makefiles. Diese
verstehen als Regeln "all", "clean", "distclean" sowie die jeweiligen
Namen der executables und Objektdateien.


---------------------------------------
  Hinweise fuer die DROPS-Entwickler:
---------------------------------------

Damit die automatische Generierung der Abhaengigkeiten funktioniert, muessen
im Quelltext alle eingebundenen DROPS-Header-Files *immer* mit Pfadangabe 
versehen sein (auch wenn das Header-File im selben Verzeichnis steht!). 
Also muss z.B. in geom/boundary.cpp stehen: 
#include "geom/boundary.h"    statt    #include "boundary.h"

Wenn sich die Abhaengigkeiten geaendert haben, koennen diese automatisch 
mit Hilfe des top-level-Makefiles neu erzeugt werden: Dazu muss im 
DROPS-Rootverzeichnis der Befehl
    make dep
ausgefuehrt werden.

Wenn neue executables hinzugekommen sind, muessen diese im jeweiligen lokalen
Makefile eingetragen werden, indem sie der Variablen EXEC hinzugefuegt werden
und eine neue Regel zum Linken des executable angelegt wird.

Kleiner Makefile-Chrash-Kurs: 
------------------------------
Eine Regel sieht so aus:

    target: depend1 depend2 ...
       <Tab>   command

target und depends koennen auch Patterns enthalten; das Zeichen % steht dabei
als Platzhalter. 
Weitere nuetzliche Anwendung von Patterns: Ersetzung, z.B.
    FILES = xxx.c yyy.c zzz.c
    OBJ   = $(FILES:%.c=%.o)     # -> xxx.o yyy.o zzz.o
       
Automatische Variablen, die in command verwendet werden koennen:

$@ = target
$< = erste Abhaengigkeit depend1
$^ = alle Abhaengigkeiten
$* = target ohne Suffix, bzw. = % bei Patterns

ACHTUNG: Automatische Variablen koennen nicht in conditionals
    ifeq "..." "..."
verwendet werden!

Besteht der Wert einer automatischen Variable aus Verzeichnis und Dateiname,
so kann mit $(@D) bzw. $(@F) der directory- bzw. file-Anteil zurueckgegeben 
werden (in diesem Beispiel fuer die automatische Variable $@).

Referenzen:
- How to write a Makefile:  http://vertigo.hsrl.rutgers.edu/ug/make_help.html
- GNU Make:                 http://www.gnu.org/manual/make/index.html

