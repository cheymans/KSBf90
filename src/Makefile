FC = ifort
FFLAGS = -CB -O2

LIB = -L/star/lib -lpgplot -L/usr/X11 -lX11 -L/usr/lib64 -lcfitsio

OBJS = Calculate_Moments.o Permutations.o

gal_correct:$(OBJS)
	$(FC) $(FFLAGS) $@.f90 -o $@.a $(OBJS) $(LIB)

psffit: $(OBJS)
	$(FC) $(FFLAGS) $@.f90 -o $@.a $(OBJS) $(LIB)

.f90.o: $(OBJS)
	$(FC) -c $(FFLAGS) $< $(LIB)

clean:
	rm *.o