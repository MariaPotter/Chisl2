#!/usr/bin/env bash

# Найти все файлы .c, .cpp, .cc и т.п. и подать их на вход компилятору
find . -maxdepth 1 -name '*.c*' | xargs -r g++ -Wall -Wextra -lm -g -O -fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize=null -fsanitize=bounds-strict -fstack-protector-all
