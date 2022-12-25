#!/usr/bin/env bash
ERRORS=0
function countError {
    echo "Проверка завершилась с ошибкой."
    (( ERRORS++ ))
}
trap countError ERR 

echo "# cppcheck"
cppcheck --enable=warning,performance --error-exitcode=2 . && echo "OK" || false

echo "# clang-tidy"
# "--" should be last clang-tidy argument,
# so we concatenate output of "find" with "--"
find . -maxdepth 1 -name "*.c*" | xargs -r $(dirname $0)/clang-tidy.sh  && echo "OK" || false

echo "# Предупреждения GCC"
find . -maxdepth 1 -name "*.c*" | xargs -r g++ -Wall -Werror -Wextra -lm -g -O -fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize=null -fsanitize=bounds-strict -fstack-protector-all && echo "OK" || false

exit $ERRORS
