
ruamel.std.zipfile
==================

`package ruamel.std.zipfile <https://bitbucket.org/ruamel/std.zipfile>`_ is a drop-in
improvements over the standard zipfile package

You can just replace::

  import zipfile

with::

  import ruamel.std.zipfile as zipfile


The package includes InMemoryZipFile, which allows easy creation of ZIP files in memory,
and allows for streaming or writing out to disc after full creation::

  with InMemoryZipFile() as imz:
      imz.append("test.txt", "Another test").append("test2.txt", "Still another")
      with open('some_file.zip', 'wb') as fp:
          fp.write(imz.data)

which will write a two file ZIP file, the first file of which is named
``test.txt`` with content ``Another test``. The ``.data`` content can
also be dynamically returned to a web browser, etc.


File deletion from ZIP
----------------------

Taking advantage of the delayed writing of ``InMemoryZipFile``, the
function ``delete_from_zip_file(file_name, pattern, file_names)``,
takes a ``string`` or ``pathlib.Path`` as file_name.

Any files matching the pattern, if provided, are deleted from the
file, as well as are any files matching ``file_names`` (a list of
string/Path, single non-list instances are allowed).

The following deletes any files ending in `.pth` and the file
`tmp/README.rst` from a ZIP file ``test.zip``::

  delete_from_zip_file('test.zip', pattern='.*.pth', file_names=['README.rst'])

or::

  delete_from_zip_file('test.zip', pattern='.*.pth', file_names='README.rst')

or::

  delete_from_zip_file('test.zip', pattern=re.compile('.*.pth'), file_names='README.rst')


*Please note that this a ``re`` pattern not a ``glob`` pattern.*
You can, but don't have to provide a pattern compiled with
``re.compile()``


The ZIP file is recreated and recompressed, take this into account
when deleting files (restrict the size of files you handle, combine
patterns instead of doing multiple calls).
