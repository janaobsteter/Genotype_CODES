import sys

oldwriter = sys.stdout

class NewWriter:
    def write(self, data):
        # pyqtwindow.textarea.text += data
        oldwriter.write(data.upper())

sys.stdout = NewWriter()

print 'Test'

def baz(num, location, place):
    print 'number', num
    print 'location', location
    print 'place', place

def foo(name, *args, **kwargs):
    print name
    # baz(*args, **kwargs)
    print 'args', args
    print 'kwargs', kwargs
    kwargs['place']

foo('bar', 123, 'abc', place='x')

