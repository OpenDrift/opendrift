class fakeproj:
    # For readers with unprojected domain, we emulate a
    # pyproj class with needed functions
    class _crs:
        is_geographic = False

    crs = _crs()

    def __call__(self, x, y, inverse=False):
        # Simply return x and y since these are also row/column indices
        return x, y

    def transform_bounds(self, left, bottom, right, top, *args, **kwargs):
        """
        Fake transform_bounds
        """
        return (left, bottom, right, top)

    srs = 'None'


