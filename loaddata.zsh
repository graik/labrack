#!/usr/bin/env zsh
# load some useful initial data (units and perhaps some examples)
# into a labrack database
#
# run with the shell that is configured for the django project
# ($DJANGO_SETTINGS_MODULE and labrack in $PYTHONPATH)

python manage.py loaddata initial_units.json
python manage.py loaddata initial_types.json
