.PHONY:	main
.PHONY:	genKeys
.PHONY: test

all: main genKeys

main:
	mpicc -lflint -lcrypto -lgmp -g -o main main.c util.c algo.c

test:
	mpicc -lflint -lcrypto -lgmp -g -o test test.c util.c algo.c
	./test

genKeys:
	mpicc -lflint -lcrypto -lgmp -o genKeys genKeys.c util.c algo.c
