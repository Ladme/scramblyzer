scramblyzer: main.c src/general.c src/composition.c src/rate.c src/flipflops.c
	gcc src/general.c src/composition.c src/rate.c src/flipflops.c main.c -I$(groan) -L$(groan) -D_POSIX_C_SOURCE=200809L -o scramblyzer -lgroan -lm -std=c99 -pedantic -Wall -Wextra -O3 -march=native

install: scramblyzer
	cp scramblyzer ${HOME}/.local/bin
