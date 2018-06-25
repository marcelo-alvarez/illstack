#!/bin/bash

git pull
comment="$@"
if [ "$#" -gt 0 ]; then
    git commit -a -m "$comment"
else
    git commit -a
fi
git pull
git push -u origin
