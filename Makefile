# Riemann solver.

# No need for anything impressive.  This is a simple
# code with only one source file.

CXX=g++

SOURCE=main.C Riemann.C Cloption.C
OBJECTS=main.o Riemann.o Cloption.o
INCLUDES=riemann.h cloption.h eos.h
EXEC=rshock_riemann

INSTALLDIR="$(HOME)/gcc/bin"

CXXFLAGS=-DDEBUG -g

$(EXEC): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXEC)

# Compilation: just use the built-in rules (gnu make)

Riemann.o: Riemann.C Cloption.C $(INCLUDES)

# Useful targets

install: $(EXEC)
	mv $(EXEC) $(INSTALLDIR)

uninstall:
	rm $(INSTALLDIR)/$(EXEC)

clean:
	rm $(OBJECTS)
	make uninstall
	make emacsclean

rcs:
	ci -l $(SOURCE) $(INCLUDES) Makefile README

emacsclean:
	rm -f *~

dist=RiemannObject.1.0
distfile=$(dist).tar
distdir=$(dist)

dist:
	mkdir $(dist)
	cp $(SOURCE) $(INCLUDES) Makefile README $(distdir)
	tar -cvf $(distfile) ./$(distdir)
	gzip $(dist).tar
	rm -rf $(distdir)
