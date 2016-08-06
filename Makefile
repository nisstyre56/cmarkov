default: markov.c markov.h
	$(MAKE) lib;
	$(CC) -g -DTOK_LIB -Wall -Wextra -std=gnu99 -Wpointer-arith -Wmissing-prototypes -Werror -lmaa -lm -L. -ltokenize -O3 ./markov.c -o markov_test -Wl,-rpath,/home/wes/markov;

unsafe: markov.c  markov.h
	$(MAKE) lib;
	$(CC) -DNDEBUG -DTOK_LIB -Wall -std=gnu99 -Wextra -Wpointer-arith -Wmissing-prototypes -lmaa -lm -L. -ltokenize -O3 ./markov.c -o markov -Wl,-rpath,/home/wes/markov;

lib: markov.c markov.h tokenize.c tokenize.h
	$(CC) -g -DTOK_LIB -c -fpic -Wall -Wextra -pedantic -Wpointer-arith -Werror -std=c99 -O3 ./tokenize.c
	$(CC) -shared -o libtokenize.so tokenize.o -lmaa;
	$(CC) -g -DTOK_LIB -c -fpic -Wall -Wextra -pedantic -Wpointer-arith -Werror -std=gnu99 -O3 ./markov.c -lmaa -lm
	$(CC) -shared -o markov.so markov.o -lmaa;
