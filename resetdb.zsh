#!/usr/bin/env zsh
# remove all tables from database and start from scratch 
# using the SOUTH migration system
#
# Note:
#    manage.py reset is deprecated and not available from django 1.5
#
# run with the shell that is configured for the django project
# ($DJANGO_SETTINGS_MODULE and labrack in $PYTHONPATH)

python manage.py migrate labrack zero
python manage.py migrate