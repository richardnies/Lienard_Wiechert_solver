IDIR = ./
CXX = g++
CXXFLAGS = -Wall -Werror -pedantic --std=c++11 $(debug_or_optimize)

ODIR = obj
EXEDIR = ..

_DEPS = particle.h simulation.h array_helpers.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = particle.o simulation.o array_helpers.o LWS_driver.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(EXEDIR)/simulation: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 