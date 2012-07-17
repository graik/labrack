from django import template

register = template.Library()

@register.simple_tag
def hasWritePermission(original, user):
    ## Ref: http://thedjangoforum.com/board/thread/592/function-call-in-templates/
    original.writePermission(user)
    return ''
