class NotCoveredError(Exception):
    pass

class OutsideSpatialCoverageError(NotCoveredError):
    pass

class OutsideTemporalCoverageError(NotCoveredError):
    pass

class VariableNotCoveredError(NotCoveredError):
    pass
