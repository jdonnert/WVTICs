#!/bin/bash
pwd
astyle --style=kr --attach-inlines --indent=spaces=4 --suffix=none --pad-paren --pad-oper --pad-header --align-pointer=name --align-reference=name --add-brackets *.h *.c
