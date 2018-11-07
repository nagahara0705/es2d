FC=gfortran

TARGET=./inflow.out

$(TARGET):film_inflow.f90
	$(FC) -o $(TARGET) film_inflow.f90 
