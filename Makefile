c_test : src/C/jgs.h src/C/matrice.h src/C/gauss_seidel.c src/C/jacobi.c src/C/matrice.c tests/C/test.c
	gcc -o c_test tests/C/test.c src/C/gauss_seidel.c src/C/jacobi.c src/C/matrice.c -lm
	./c_test

p_test :
	python3 tests/Python/test.py

p_resol_ex1 :
	python3 tests/Python/resol_ex1.py

clean :
	rm c_test