import numpy as np

class Savez(object):
    def __init__(self, zipfilename):
        import zipfile, tempfile, os, sys

        if isinstance(zipfilename, basestring):
            if not zipfilename.endswith('.npz'):
                zipfilename += '.npz'

        # original _savez has no compression
        compression = zipfile.ZIP_STORED

        if sys.version_info >= (2, 5):
            self.zip = zipfile.ZipFile(
                zipfilename, mode="w", allowZip64=True,
                compression=compression)

        # Stage arrays in a temporary file on disk,
        # before writing to zip.
        fd, tmpfile = tempfile.mkstemp(suffix='-numpy.npy')
        os.close(fd)
        self.tmpfile = tmpfile
        self.i = 0  # array counter

    def savez(self, *args, **kwds):
        import os
        import numpy.lib.format as format

        namedict = kwds
        for val in args:
            key = 'arr_%d' % self.i
            if key in namedict.keys():
                raise ValueError(
              "Cannot use un-named variables and keyword %s" % key)
            namedict[key] = val
            self.i += 1

        try:
            for key, val in namedict.iteritems():
                fname = key + '.npy'
                fid = open(self.tmpfile, 'wb')
                try:
                    format.write_array(fid, np.asanyarray(val))
                    fid.close()
                    fid = None
                    self.zip.write(self.tmpfile, arcname=fname)
                finally:
                    if fid:
                        fid.close()
        finally:
            os.remove(self.tmpfile)

    def close(self):
        self.zip.close()

def test_Savez():
    import tempfile, os
    tmp = 'tmp_testarchive'
    database = Savez(tmp)
    for i in range(4):
        array = np.linspace(0, 5+i, 3)
        kwargs = {'myarray_%02d' % i: array}
        database.savez(**kwargs)
    database.close()

    database = np.load(tmp+'.npz')

    expected = {
        'myarray_00': np.array([ 0. ,  2.5,  5. ]),
        'myarray_01': np.array([ 0.,  3.,  6.])
        'myarray_02': np.array([ 0. ,  3.5,  7. ]),
        'myarray_03': np.array([ 0.,  4.,  8.]),
        }
    for name in database:
        computed = database[name]
        diff = np.abs(expected[name] - computed).max()
        assert diff < 1E-13
    database.close
    os.remove(tmp+'.npz')

test_Savez()
