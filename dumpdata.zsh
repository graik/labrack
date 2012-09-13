#!/usr/bin/env zsh
# dump current labrack database content (w/o user infos) into
# initial.json
# The dump can be further restricted to certain tables with,
# for example, dumpdata -i 4 labrack.sample > example_samples.json
#
# run with the shell that is configured for the django project
# ($DJANGO_SETTINGS_MODULE and labrack in $PYTHONPATH)

python manage.py dumpdata -i 4 labrack > initial.json