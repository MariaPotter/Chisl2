#!/usr/bin/env bash
for val in $(ls | grep "\.[ch]p*")
do
./.github/clang-format.sh ./"$val" > ./h_"$val" && mv ./h_"$val" ./"$val"
done

