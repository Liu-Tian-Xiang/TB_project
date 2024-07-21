CC = mpiicpx -std=c++17 -O3
LIBS = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lboost_program_options -L/usr/local/BOOST_INSTALL/lib
OUT = ecode
OBJDIR=./objects
SRCDIR=./src

LIBOBJS= \
	$(OBJDIR)/BandStructure.o \
	$(OBJDIR)/Material.o \
	$(OBJDIR)/SymmetryPoints.o \
	$(OBJDIR)/Hamiltonian.o \
	$(OBJDIR)/HamiltonianNew.o \
	$(OBJDIR)/HamiltonianTest.o \
	$(OBJDIR)/HamiltonianWorking.o \
	$(OBJDIR)/Hamiltonian_pzheevx.o \
	$(OBJDIR)/Hamiltonian_fold.o \
	$(OBJDIR)/Cell.o \
	$(OBJDIR)/main.o \


$(OUT): objdir $(LIBOBJS)
	$(CC) $(LIBOBJS) -o $@ $(LIBS) $(GSL) $(LAPACK) 


objdir:
	mkdir -p $(OBJDIR)

$(OBJDIR)/Cell.o: $(SRCDIR)/Cell.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/BandStructure.o: $(SRCDIR)/BandStructure.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/Material.o: $(SRCDIR)/Material.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/SymmetryPoints.o: $(SRCDIR)/SymmetryPoints.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/HamiltonianWorking.o: $(SRCDIR)/HamiltonianWorking.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/HamiltonianTest.o: $(SRCDIR)/HamiltonianTest.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/HamiltonianNew.o: $(SRCDIR)/HamiltonianNew.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/Hamiltonian.o: $(SRCDIR)/Hamiltonian.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/Hamiltonian_fold.o: $(SRCDIR)/Hamiltonian_fold.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/Hamiltonian_pzheevx.o: $(SRCDIR)/Hamiltonian_pzheevx.cpp
	$(CC) -g -c $< -o $@
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp
	$(CC) -g -c $< -o $@



clean:
	rm -rf $(OBJDIR) $(OUT)
