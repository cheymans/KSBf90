FC = ifort
FFLAGS = -CB -O2

LIB = -L/star/lib -lpgplot -L/usr/X11 -lX11 -L/usr/lib64 -lcfitsio

OBJS = Moments.o

psffit: $(OBJS)
	rm $@.a
	$(FC) $(FFLAGS) $@.f90 -o $@.a $(OBJS) $(LIB)

clean:
	rm *.o