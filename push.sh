#!/bin/bash

comment="$@"
if [ "$#" -gt 0 ]; then
    git commit -a -m "$comment"
else
    git commit -a
fi
git push -u origin
