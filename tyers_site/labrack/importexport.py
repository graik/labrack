import csv
from django.http import HttpResponse


def generate_csv(self, request, queryset, fields, tableName):
    """
    Export selected objects as CSV file
    """

    response = HttpResponse(mimetype='text/csv')
    response['Content-Disposition'] = 'attachment; filename=' + tableName + '.csv'

    writer = csv.writer(response)
    writer.writerow(fields.keys())

    for obj in queryset:
        columns = []

        for name,value in fields.items():
            try:
                #columns.append( obj.__getattribute__(value))  
                ## can't get something like container.location
                columns.append( eval('obj.%s'%value) )
            except:
                columns.append("")  ## capture 'None' fields

        columns = [ c.encode('utf-8') if type(c) is unicode else c \
                    for c in columns]

        writer.writerow( columns )

    return response