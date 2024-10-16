include makefile.inc
all: res1 res2 res3 res4
Objs1= driver1.o src.o mpset.o epc_mpupk.o

driver1: $(Objs1)
	$(F77LINK) $(F77LINKOPTS) -o driver1 $(Objs1) 
res1: driver1  
	driver1  >res1


Objs2= driver2.o src.o mpset.o epc_mpupk.o

driver2: $(Objs2)
	$(F77LINK) $(F77LINKOPTS) -o driver2 $(Objs2) 
res2: driver2  
	driver2  >res2


Objs3= driver3.o src.o mpset.o epc_mpupk.o

driver3: $(Objs3)
	$(F77LINK) $(F77LINKOPTS) -o driver3 $(Objs3) 
res3: driver3  
	driver3  >res3


Libs4= $(PORT) $(TIMER)
Objs4= driver4.o src.o mpset.o epc_mpupk.o

driver4: $(Objs4)
	$(F77LINK) $(F77LINKOPTS) -o driver4 $(Objs4)  $(Libs4)
res4: driver4 data4 
	driver4 <data4 >res4


