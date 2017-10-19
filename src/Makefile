guc : get_unit_cell
	gfortran guc.o read_data.o distf1.o supercell.o intersection.o distances.o check_angles.o matrix_operations.o check_cell.o write_output.o -o guc
	rm *.o
	rm *.mod

get_unit_cell : distf1 read_data supercell supercell intersection distances check_angles matrix_operations check_cell write_output
	gfortran -c get_unit_cell.f90 -o guc.o

distf1 : distf1.f90
	gfortran -c distf1.f90 -o distf1.o

read_data : read_data.f90
	gfortran -c read_data.f90 -o read_data.o

supercell : supercell.f90
	gfortran -c supercell.f90 -o supercell.o

intersection : intersection.f90
	gfortran -c intersection.f90 -o intersection.o

distances : distances.f90
	gfortran -c distances.f90 -o distances.o

check_angles : check_angles.f90 matrix_operations
	gfortran -c check_angles.f90  -o check_angles.o

check_cell : check_cell.f90 matrix_operations
	gfortran -c check_cell.f90  -o check_cell.o

matrix_operations : matrix_operations.f90
	gfortran -c matrix_operations.f90 -o matrix_operations.o

write_output : write_output.f90
	gfortran -c write_output.f90 -o write_output.o
