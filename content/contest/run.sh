#!/bin/bash
g++ -std=c++20 -DBRUTE -O2 -Wall -Wextra -Wconversion -Wfatal-errors -fsanitize=address,undefined $1 && ./a.out
